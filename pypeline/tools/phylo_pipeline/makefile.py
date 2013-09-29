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
     IsUnsignedInt, \
     IsBoolean, \
     StringStartsWith, \
     StringEndsWith, \
     CLI_PARAMETERS, \
     And, \
     Or, \
     Not
from pypeline.common.fileutils import \
     swap_ext



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
    _update_exclusions(mkfile)
    _check_genders(mkfile)
    _check_max_read_depth(mkfile)
    _check_indels_and_msa(mkfile)
    mkfile["Nodes"] = ()

    padding = mkfile["Genotyping"]["Padding"]
    mkfile["Genotyping"]["Random"]["--padding"] = padding

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
        subdd["FASTA"]  = os.path.join(options.genomes_root, subdd["Prefix"] + ".fasta")

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
    def _collect_subsets(tree, path):
        for (roi, subset) in tree.get("SubsetRegions", {}).iteritems():
            if roi not in subsets_by_regions:
                raise MakefileError("Subset of unknown region (%r) requested at %r" % (roi, path))

            roi_fname = swap_ext(subsets_by_regions[roi]["BED"], subset + ".names")
            if not os.path.isfile(roi_fname):
                raise MakefileError(("Subset file does not exist Regions Of Interest:\n"
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

            subsets_by_regions[roi]["Sequences"][subset] = frozenset(sequences)

    _collect_subsets(mkfile["PhylogeneticInference"], "PhylogeneticInference")
    _collect_subsets(mkfile["PAML"]["codeml"], "PAML:codeml")


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
    max_depths = mkfile["Genotyping"]["VCF_Filter"]["MaxReadDepth"]
    if isinstance(max_depths, types.DictType):
        required_keys = set()
        for sample in mkfile["Project"]["Samples"].itervalues():
            if sample["GenotypingMethod"].lower() == "samtools":
                required_keys.add(sample["Name"])

        # Extra keys are allowed, to make it easier to temporarily disable a sample
        missing_keys = required_keys - set(max_depths)
        if missing_keys:
            raise MakefileError("MaxReadDepth not specified for the following samples:\n    - %s" \
                                % ("\n    - ".join(sorted(missing_keys)),))


def _check_indels_and_msa(mkfile):
    if mkfile["MSAlignment"]["Enabled"]:
        return

    regions = mkfile["Project"]["Regions"]
    for (name, subdd) in regions.iteritems():
        if subdd["IncludeIndels"]:
            raise MakefileError("Regions %r includes indels, but MSA is disabled!" % (name,))


def _update_exclusions(mkfile):
    samples = mkfile["Project"]["Samples"]
    groups  = mkfile["Project"]["Groups"]

    mkfile["PhylogeneticInference"]["ExcludeSamples"] = \
      _select_samples(mkfile["PhylogeneticInference"]["ExcludeSamples"], groups, samples, "PhylogeneticInference:ExcludeSamples")

    mkfile["PAML"]["codeml"]["ExcludeSamples"] = \
      _select_samples(mkfile["PAML"]["codeml"]["ExcludeSamples"], groups, samples, "PAML:codeml:ExcludeSamples")


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
        "Padding"  : IsUnsignedInt(default = 5),
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
            "MaxReadDepth"  : Or(IsUnsignedInt, IsDictOf(IsStr, IsUnsignedInt),
                                 default = 0),
            StringStartsWith("-") : CLI_PARAMETERS,
            },
        },
    "MSAlignment" : {
        "Enabled"   : IsBoolean(default = True),
        "Default"   : StringIn(("mafft",), # TODO: Add support for other programs
                               default = "mafft"),
        "MAFFT" : {
            "Algorithm" : StringIn(("auto", "FFT-NS-1", "FFT-NS-2", "FFT-NS-i", "NW-INS-i", "L-INS-i", "E-INS-i", "G-INS-i"),
                                   default = "auto")
            },
        },
    "PhylogeneticInference" : {
        "ExcludeSamples" : IsListOf(IsStr, default = []),
        "SubsetRegions"  : IsDictOf(IsStr, IsStr, default = {}),
        "Default" : StringIn(("examl",), # TODO: Add support for other programs
                             default = "examl"),
        "ExaML" : {
            "Bootstraps" : IsUnsignedInt(default = 100),
            "Replicates" : IsUnsignedInt(default = 1),
            "Model"      : StringIn(("GAMMA", "PSR"),
                                    default = "gamma"),
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
