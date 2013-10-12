#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import os
import types
import pysam

import pypeline.common.makefile
from pypeline.common.makefile import \
     MakefileError, \
     REQUIRED_VALUE, \
     IsStr, \
     IsDictOf, \
     IsListOf, \
     StringIn, \
     IsFloat, \
     IsUnsignedInt, \
     IsBoolean, \
     IsNone, \
     ValueIn, \
     ValuesSubsetOf, \
     StringStartsWith, \
     StringEndsWith, \
     CLI_PARAMETERS, \
     And, \
     Or, \
     Not
from pypeline.common.fileutils import \
     swap_ext
from pypeline.common.utilities import \
     fill_dict



def read_makefiles(options, filenames):
    makefiles = []
    for filename in filenames:
        makefile = pypeline.common.makefile.read_makefile(filename, _VALIDATION)
        makefile = _mangle_makefile(options, makefile["Makefile"])
        makefiles.append(makefile)
    return makefiles


def _mangle_makefile(options, mkfile):
    _collapse_samples(mkfile)
    _update_regions(options, mkfile)
    _update_subsets(options, mkfile)
    _update_filtering(mkfile)
    _update_sample_sets(mkfile)
    _update_genotyping(mkfile)
    _update_msa(mkfile)
    _check_genders(mkfile)
    _check_max_read_depth(mkfile)
    _check_indels_and_msa(mkfile)
    mkfile["Nodes"] = ()

    return mkfile


def _collapse_samples(mkfile):
    groups, samples = {}, set()
    def _collect_samples(samples_dict, path = ()):
        current_samples = {}
        for (key, subdd) in samples_dict.iteritems():
            if key.startswith("<") and key.endswith(">"):
                key = key.lstrip("<").rstrip(">")
                current_samples.update(_collect_samples(subdd, path + (key,)))
            elif key not in samples:
                samples.add(key)
                subdd["Name"] = key
                current_samples[key] = subdd
            else:
                raise MakefileError("Duplicate sample-name: %r" % (key,))

        groups[path] = current_samples
        return current_samples

    _collect_samples(mkfile["Project"]["Samples"])
    mkfile["Project"]["Samples"] = groups.pop(())
    mkfile["Project"]["Groups"] = groups


def _select_samples(select, groups, samples, path):
    selection = set()
    for group in select:
        if group.startswith("<") and group.endswith(">"):
            key = tuple(group[1:-1].split("/"))
            if key not in groups:
                raise MakefileError("Unknown group specifed for filtering %r: %r" % (path, key))
            selection.update(groups[key])
        elif group in samples:
            selection.add(group)
        else:
            raise MakefileError("Unknown/Invalid group specifed for filtering %r: %r" % (path, group))
    return selection


def _update_regions(options, mkfile):
    mkfile["Project"]["Regions"] = mkfile["Project"].pop("RegionsOfInterest")

    for (name, subdd) in mkfile["Project"]["Regions"].iteritems():
        if "Prefix" not in subdd:
            raise MakefileError("No genome specified for regions %r" % (name,))

        subdd["Name"]   = name
        subdd["Desc"]   = "{Prefix}.{Name}".format(**subdd)
        subdd["BED"]    = os.path.join(options.regions_root, subdd["Desc"] + ".bed")
        subdd["FASTA"]  = os.path.join(options.prefix_root, subdd["Prefix"] + ".fasta")

        required_files = (
            ("Regions file", subdd["BED"], None),
            ("Reference sequence", subdd["FASTA"], None),
            ("Reference sequence index", subdd["FASTA"] + ".fai",
             "Please index using 'samtools faidx %s'" % (subdd["FASTA"],)))

        for (desc, path, instructions) in required_files:
            if not os.path.isfile(path):
                message = "%s does not exist for %r:\n  Path = %r" \
                                % (desc, name, path)
                if instructions:
                    message = "%s\n%s" % (message, instructions)
                raise MakefileError(message)

        # Collects seq. names / validate regions
        subdd["Sequences"] = {None : _collect_and_validate_sequences_and_subsets(subdd)}
        subdd["SubsetFiles"] = {None : ()}

        sampledd = subdd["Genotypes"] = {}
        for sample_name in mkfile["Project"]["Samples"]:
            fasta_file = ".".join((sample_name, subdd["Desc"], "fasta"))
            sampledd[sample_name] = os.path.join(options.destination,
                                              mkfile["Project"]["Title"],
                                              "genotypes",
                                              fasta_file)


def _collect_and_validate_sequences_and_subsets(regions):
    contigs = {}
    with open(regions["FASTA"] + ".fai") as faihandle:
        for line in faihandle:
            name, length, _ = line.split(None, 2)
            if name in contigs:
                raise MakefileError(("Reference contains multiple identically named sequences:\n"
                                     "  Path = %r\n  Name = %r\n"
                                     "Please ensure that all sequences have a unique name!")
                                     % (regions["FASTA"], name))
            contigs[name] = int(length)

    parser = pysam.asBed()
    sequences = set()
    with open(regions["BED"]) as bedhandle:
        for (line_num, line) in enumerate(bedhandle):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            try:
                bed = parser(line, len(line))
                # Force evaluation of (lazily parsed) properties
                bed_start = bed.start
                bed_end   = bed.end
            except ValueError, error:
                raise MakefileError(("Error parsing line in regions file:\n"
                                     "  File = %r\n  Line = %i\n %s")
                                     % (regions["BED"], line_num, error))

            if len(bed) < 6:
                raise MakefileError(("Region at line #%i (%s) does not contain enough fields;\n"
                                     "at least the first 6 fields are required. C.f. defination at\n"
                                     "  http://genome.ucsc.edu/FAQ/FAQformat.html#format1")
                                     % (line_num, repr(bed.name) if len(bed) > 3 else "unnamed record"))

            contig_len = contigs.get(bed.contig)
            if contig_len is None:
                raise MakefileError(("Regions file contains contigs not found in reference:\n"
                                     "  Path = %r\n  Name = %r\n"
                                     "Please ensure that all contig names match the reference names!")
                                     % (regions["BED"], bed.contig))
            elif not (0 <= int(bed_start) < int(bed_end) <= contig_len):
                raise MakefileError(("Regions file contains invalid region:\n"
                                     "  Path   = %r\n  Contig = %r\n"
                                     "  Start  = %s\n  End    = %s\n"
                                     "Start must be >= 0 and < End, and End must be <= %i!")
                                     % (regions["BED"], bed.contig, bed.start, bed.end, contig_len))
            elif bed.strand not in "+-":
                raise MakefileError(("Regions file contains invalid region: "
                                     "  Path   = %r\n  Line = %i\n  Name = %r\n"
                                     "Strand is %r, expected either '+' or '-'")
                                     % (regions["BED"], line_num, bed.name, bed.strand))

            sequences.add(bed.name)
    return frozenset(sequences)


def _update_subsets(_options, mkfile):
    subsets_by_regions = mkfile["Project"]["Regions"]
    def _collect_subsets(roi, subset, path):
        if roi not in subsets_by_regions:
            raise MakefileError("Subset of unknown region (%r) requested at %r" % (roi, path))

        roi_fname = swap_ext(subsets_by_regions[roi]["BED"], subset + ".names")
        if not os.path.isfile(roi_fname):
            raise MakefileError(("Subset file does not exist for Regions Of Interest:\n"
                                 "  Region = %r\n  Subset = %r\n  Path   = %r")
                                 % (roi, subset, roi_fname))

        sequences = set()
        with open(roi_fname) as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    sequences.add(line)

        known_seqs   = subsets_by_regions[roi]["Sequences"][None]
        unknown_seqs = sequences - known_seqs
        if unknown_seqs:
            message = ("Unknown sequences in subset file:\n"
                       "  File   = %r\n  Region = %r\n  Subset = %r\n"
                       "  Unknown sequence names =") \
                       % (roi_fname, roi, subset)
            unknown_seqs = list(sorted(unknown_seqs))
            if len(unknown_seqs) > 5:
                unknown_seqs = unknown_seqs[:5] + ["..."]
            message = "\n    - ".join([message] + unknown_seqs)
            raise MakefileError(message)

        subsets_by_regions[roi]["SubsetFiles"][subset] = (roi_fname,)
        subsets_by_regions[roi]["Sequences"][subset] = frozenset(sequences)


    for (key, subdd) in mkfile["PhylogeneticInference"].iteritems():
        for (subkey, roidd) in subdd["RegionsOfInterest"].iteritems():
            roidd["Name"] = subkey

            if roidd.get("SubsetRegions") is not None:
                path = "PhylogeneticInference:%s:RegionsOfInterest:%s" % (key, subkey)
                _collect_subsets(subkey, roidd["SubsetRegions"], path)

    for (roi, subset) in mkfile["PAML"]["codeml"]["SubsetRegions"].iteritems():
        _collect_subsets(roi, subset, "PAML:codeml:SubsetRegions")


def _update_filtering(mkfile):
    samples = mkfile["Project"]["Samples"]
    groups  = mkfile["Project"]["Groups"]

    filtering = {}
    for (target, filter_by) in mkfile["Project"]["FilterSingletons"].iteritems():
        if target not in samples:
            raise MakefileError("Unknown/Invalid group specifed for singleton filtering: %r" % (target,))

        path = "Project:FilterSingletons:%s" % (target,)
        filtering[target] = _select_samples(filter_by, groups, samples, path)

    mkfile["Project"]["FilterSingletons"] = filtering


def _check_genders(mkfile):
    regions_genders = set()
    for regions in mkfile["Project"]["Regions"].itervalues():
        current_genders = set(regions["HomozygousContigs"])
        if not regions_genders:
            regions_genders = current_genders
        elif regions_genders != current_genders:
            raise MakefileError("List of genders for regions %r does not match other regions" \
                                % (regions["Name"],))

    for sample in mkfile["Project"]["Samples"].itervalues():
        if sample["Gender"] not in regions_genders:
            raise MakefileError("Sample %r has unknown gender %r; known genders are %s" \
                                % (sample["Name"], sample["Gender"],
                                   ", ".join(map(repr, regions_genders))))


def _check_max_read_depth(mkfile):
    for (key, settings) in mkfile["Genotyping"].iteritems():
        max_depths = settings["VCF_Filter"]["MaxReadDepth"]
        if isinstance(max_depths, types.DictType):
            required_keys = set()
            for sample in mkfile["Project"]["Samples"].itervalues():
                if sample["GenotypingMethod"].lower() == "samtools":
                    required_keys.add(sample["Name"])

            # Extra keys are allowed, to make it easier to temporarily disable a sample
            missing_keys = required_keys - set(max_depths)
            if missing_keys:
                raise MakefileError("MaxReadDepth not specified for the following samples for %r:\n    - %s" \
                                    % (key, "\n    - ".join(sorted(missing_keys))))


def _check_indels_and_msa(mkfile):
    msa     = mkfile["MultipleSequenceAlignment"]
    regions = mkfile["Project"]["Regions"]
    for (name, subdd) in regions.iteritems():
        msa_enabled = msa[name]["Enabled"]

        if subdd["IncludeIndels"] and not msa_enabled:
            raise MakefileError("Regions %r includes indels, but MSA is disabled!" % (name,))


def _update_sample_sets(mkfile):
    samples = mkfile["Project"]["Samples"]
    groups  = mkfile["Project"]["Groups"]

    for (key, subdd) in mkfile["PhylogeneticInference"].iteritems():
        subdd["ExcludeSamples"] = \
          _select_samples(subdd["ExcludeSamples"], groups, samples, "PhylogeneticInference:%s:ExcludeSamples" % (key,))

        # Replace None with an empty list, to simplify code using this value
        root_trees_on = subdd["RootTreesOn"] or ()
        subdd["RootTreesOn"] = \
          _select_samples(root_trees_on, groups, samples, "PhylogeneticInference:%s:RootTreesOn" % (key,))

    mkfile["PAML"]["codeml"]["ExcludeSamples"] = \
      _select_samples(mkfile["PAML"]["codeml"]["ExcludeSamples"], groups, samples, "PAML:codeml:ExcludeSamples")


def _update_genotyping(mkfile):
    genotyping = mkfile["Genotyping"]
    defaults   = genotyping.pop("Defaults")
    defaults.setdefault("Padding", 5)
    defaults["VCF_Filter"].setdefault("MaxReadDepth", 0)

    for key in mkfile["Project"]["Regions"]:
        subdd = fill_dict(genotyping.get(key, {}), defaults)
        subdd["Random"]["--padding"] = subdd["Padding"]
        genotyping[key] = subdd

    regions = set(genotyping)
    unknown_regions = regions - set(mkfile["Project"]["Regions"])
    if unknown_regions:
        raise MakefileError("Unknown Regions of Interest in Genotyping: %s" \
                            % (", ".join(unknown_regions),))


def _update_msa(mkfile):
    msa      = mkfile["MultipleSequenceAlignment"]
    defaults = msa.pop("Defaults")
    defaults.setdefault("Program", "MAFFT")
    defaults["MAFFT"].setdefault("Algorithm", "MAFFT")

    for key in mkfile["Project"]["Regions"]:
        msa[key] = fill_dict(msa.get(key, {}), defaults)

    unknown_regions = set(msa) - set(mkfile["Project"]["Regions"])
    if unknown_regions:
        raise MakefileError("Unknown Regions of Interest in Genotyping: %s" \
                            % (", ".join(unknown_regions),))


# Recursive definition of sample tree
_VALIDATION_SUBSAMPLE_KEY = And(StringStartsWith("<"),
                                StringEndsWith(">"))
_VALIDATION_SAMPLES_KEY    = And(IsStr, Not(_VALIDATION_SUBSAMPLE_KEY))
_VALIDATION_SAMPLES = {
    _VALIDATION_SAMPLES_KEY : {
        "GenotypingMethod" : StringIn(("reference sequence", "random sampling", "samtools"),
                                       default = "samtools"),
        "SpeciesName"      : IsStr,
        "CommonName"       : IsStr,
        "Gender"           : IsStr(default = REQUIRED_VALUE),
    }
}
_VALIDATION_SAMPLES[_VALIDATION_SUBSAMPLE_KEY] = _VALIDATION_SAMPLES

_VALIDATION_GENOTYPES = {
    "Padding"  : IsUnsignedInt,
    "MPileup"  : {
        StringStartsWith("-") : CLI_PARAMETERS,
    },
    "BCFTools" : {
        StringStartsWith("-") : CLI_PARAMETERS,
    },
    "Random"   : {
        "--min-distance-to-indels" : IsUnsignedInt,
    },
    "VCF_Filter" : {
        "MaxReadDepth"  : Or(IsUnsignedInt, IsDictOf(IsStr, IsUnsignedInt)),

        "-k" : IsNone,        "--keep-ambigious-genotypes"    : IsNone,
        "-q" : IsUnsignedInt, "--min-quality"                 : IsUnsignedInt,
        "-f" : IsFloat,       "--min-allele-frequency"        : IsFloat,
        "-Q" : IsUnsignedInt, "--min-mapping-quality"         : IsUnsignedInt,
        "-d" : IsUnsignedInt, "--min-read-depth"              : IsUnsignedInt,
        "-D" : IsUnsignedInt, "--max-read-depth"              : IsUnsignedInt,
        "-a" : IsUnsignedInt, "--min-num-alt-bases"           : IsUnsignedInt,
        "-w" : IsUnsignedInt, "--min-distance-to-indels"      : IsUnsignedInt,
        "-W" : IsUnsignedInt, "--min-distance-between-indels" : IsUnsignedInt,
        "-1" : IsFloat,       "--min-strand-bias"             : IsFloat,
        "-2" : IsFloat,       "--min-baseq-bias"              : IsFloat,
        "-3" : IsFloat,       "--min-mapq-bias"               : IsFloat,
        "-4" : IsFloat,       "--min-end-distance-bias"       : IsFloat,
    },
}

_VALIDATION_MSA = {
    "Enabled"   : IsBoolean(default = True),
    "Program"   : StringIn(("mafft",)), # TODO: Add support for other programs

    "MAFFT" : {
        "Algorithm" : StringIn(("mafft", "auto",
                                "FFT-NS-1", "FFT-NS-2", "FFT-NS-i",
                                "NW-INS-i", "L-INS-i", "E-INS-i", "G-INS-i")),
        StringStartsWith("-") : CLI_PARAMETERS,
    },
}


_VALIDATION = {
    "Project" : {
        "Title" : IsStr(default = "Untitled"),
        "Samples" : _VALIDATION_SAMPLES,
        "RegionsOfInterest" : {
            IsStr : {
                "Prefix"        : IsStr(default = REQUIRED_VALUE),
                "Realigned"     : IsBoolean(default = False),
                "ProteinCoding" : IsBoolean(default = False),
                "IncludeIndels" : IsBoolean(default = True),
                "HomozygousContigs" : {
                    IsStr : IsListOf(IsStr),
                    },
                },
            },
        "FilterSingletons" : {
            IsStr : IsListOf(IsStr),
            },
        },
    "Genotyping" : {
        "Defaults" : _VALIDATION_GENOTYPES,
        IsStr :      _VALIDATION_GENOTYPES,
    },
    "MultipleSequenceAlignment" : {
        "Defaults" : _VALIDATION_MSA,
        IsStr :      _VALIDATION_MSA,
        },
    "PhylogeneticInference" : {
        IsStr : {
            # Which program to use; TODO: Add support for other programs
            "Program"        : StringIn(("examl",), default = "examl"),
            # Exclude one or more samples from the phylogeny
            "ExcludeSamples" : IsListOf(IsStr, default = []),
            # Which samples to root the final trees on / or midpoint rooting
            "RootTreesOn"    : Or(IsListOf(IsStr), IsNone, default = []),
            # Create a tree per gene, for each region of interest,
            # or create a supermatrix tree from all regions specified.
            "PerGeneTrees"   : IsBoolean(default = False),
            # Selection of regions of interest / settings per region
            "RegionsOfInterest" : {
                IsStr : {
                    "Partitions"    : Or(And(IsStr, ValuesSubsetOf("123456789X")), ValueIn([False]),
                                         default = REQUIRED_VALUE),
                    "SubsetRegions" : Or(IsStr, IsNone, default = None),
                },
            },
            "SubsetRegions"  : IsDictOf(IsStr, IsStr, default = {}),
            "ExaML" : {
                "Bootstraps" : IsUnsignedInt(default = 100),
                "Replicates" : IsUnsignedInt(default = 1),
                "Model"      : StringIn(("GAMMA", "PSR"),
                                        default = "gamma"),
            }
        }
    },
    "PAML" : {
        "codeml" : {
            "ExcludeSamples" : IsListOf(IsStr, default = []),
            "SubsetRegions"  : IsDictOf(IsStr, IsStr, default = {}),
            IsStr : {
                "ControlFile" : IsStr(default = REQUIRED_VALUE),
                "TreeFile"    : IsStr(default = REQUIRED_VALUE),
            },
        },
    },
}
