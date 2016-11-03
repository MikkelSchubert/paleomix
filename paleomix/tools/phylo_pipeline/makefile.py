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
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
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
import pysam
import types

import paleomix.common.makefile
from paleomix.common.makefile import \
    MakefileError, \
    REQUIRED_VALUE, \
    IsDictOf, \
    IsListOf, \
    IsInt, \
    IsStr, \
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

from paleomix.common.fileutils import \
    swap_ext, \
    add_postfix
from paleomix.common.utilities import \
    fill_dict
from paleomix.common.console import \
    print_info, \
    print_warn
from paleomix.common.text import \
    parse_padded_table
from paleomix.common.bedtools import \
    read_bed_file, \
    BEDError
from paleomix.common.formats.fasta import \
    FASTA


def read_makefiles(options, filenames, commands):
    print_info("Reading makefile(s):")
    steps = frozenset(key for (key, _) in commands)

    makefiles = []
    for filename in filenames:
        makefile = paleomix.common.makefile.read_makefile(filename, _VALIDATION)
        makefile = _mangle_makefile(options, makefile["Makefile"], steps)
        makefiles.append(makefile)
    return makefiles


def _mangle_makefile(options, mkfile, steps):
    _collapse_samples(mkfile)
    _update_regions(options, mkfile)
    _update_subsets(mkfile, steps)
    _update_filtering(mkfile)
    _update_sample_sets(mkfile)
    _update_genotyping(mkfile)
    _update_msa(mkfile)
    _update_homozygous_contigs(mkfile)
    _check_bam_sequences(options, mkfile, steps)
    _check_sexes(mkfile)
    _update_and_check_max_read_depth(options, mkfile)
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
    print_info("    - Validating regions of interest ...")
    mkfile["Project"]["Regions"] = mkfile["Project"].pop("RegionsOfInterest")

    if not mkfile["Project"]["Regions"]:
        raise MakefileError('No regions of interest have been specified; '
                            'no analyses will be performed.')

    for (name, subdd) in mkfile["Project"]["Regions"].iteritems():
        if "Prefix" not in subdd:
            raise MakefileError("No genome specified for regions %r" % (name,))

        subdd["Name"]   = name
        subdd["Desc"]   = "{Prefix}.{Name}".format(**subdd)
        subdd["BED"]    = os.path.join(options.regions_root, subdd["Desc"] + ".bed")
        subdd["FASTA"]  = os.path.join(options.prefix_root, subdd["Prefix"] + ".fasta")

        required_files = (
            ("Regions file", subdd["BED"]),
            ("Reference sequence", subdd["FASTA"]),
        )

        for (desc, path) in required_files:
            if not os.path.isfile(path):
                raise MakefileError("%s does not exist for %r:\n  Path = %r"
                                    % (desc, name, path))

        # Collects seq. names / validate regions
        try:
            sequences = _collect_sequence_names(bed_file=subdd["BED"],
                                                fasta_file=subdd["FASTA"])
        except (IOError, BEDError), error:
            raise MakefileError("Error reading regions-of-interest %r:\n%s"
                                % (name, error))

        subdd["Sequences"] = {None: sequences}
        subdd["SubsetFiles"] = {None: ()}
        sampledd = subdd["Genotypes"] = {}
        for sample_name in mkfile["Project"]["Samples"]:
            fasta_file = ".".join((sample_name, subdd["Desc"], "fasta"))
            sampledd[sample_name] = os.path.join(options.destination,
                                                 mkfile["Project"]["Title"],
                                                 "genotypes",
                                                 fasta_file)


def _collect_fasta_contigs(filename, cache={}):
    if filename in cache:
        return cache[filename]

    if not os.path.exists(filename + ".fai"):
        print_info("      - Index does not exist for %r; this may "
                   "take a while ..." % (filename,))

    cache[filename] = contigs = dict(FASTA.index_and_collect_contigs(filename))
    return contigs


def _collect_sequence_names(bed_file, fasta_file, min_columns=6):
    contigs = _collect_fasta_contigs(fasta_file)
    sequences = {}

    for record in read_bed_file(bed_file, min_columns=6, contigs=contigs):
        current = (record.contig, record.strand)
        reference = sequences.setdefault(record.name, current)

        if current[0] != reference[0]:
            raise MakefileError("Regions in %r with the same name (%r) "
                                "are located on different contigs (%r and "
                                "%r); note that PALEOMIX assumes that "
                                "regions with the same name constitute "
                                "parts of a single consecutive sequence, "
                                "which must therefore be located on one "
                                "strand of a single sequence. Please "
                                "rename one or more of these regions to"
                                "continue." % (bed_file, record.name,
                                               current[0], reference[0]))
        elif current[1] != reference[1]:
            raise MakefileError("Regions in %r with the same name (%r) "
                                "are located on different strands; note "
                                "that PALEOMIX assumes that regions with "
                                "the same name constitute parts of a "
                                "single consecutive sequence, and that "
                                "these must therefore be located on the "
                                "same strand." % (bed_file, record.name,))

    return frozenset(sequences)


def _update_subsets(mkfile, steps):
    subsets_by_regions = mkfile["Project"]["Regions"]

    def _collect_subsets(roi, subset, path):
        if roi not in subsets_by_regions:
            raise MakefileError("Subset of unknown region (%r) requested at %r"
                                % (roi, path))

        roi_fname = swap_ext(subsets_by_regions[roi]["BED"], subset + ".names")
        if not os.path.isfile(roi_fname):
            raise MakefileError("Subset file does not exist for Regions Of "
                                "Interest:\n  Region = %r\n  Subset = %r\n"
                                "  Path   = %r"
                                % (roi, subset, roi_fname))

        sequences = set()
        with open(roi_fname) as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    sequences.add(line)

        known_seqs = subsets_by_regions[roi]["Sequences"][None]
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

    if "phylogeny:examl" in steps:
        for (key, subdd) in mkfile["PhylogeneticInference"].iteritems():
            for (subkey, roidd) in subdd["RegionsOfInterest"].iteritems():
                if subkey not in subsets_by_regions:
                    message = \
                        "Unknown regions name in phylogenetic inference:\n" \
                        "\tPath = PhylogeneticInference:%s:RegionsOfInterest" \
                        "\n\tName = %s"
                    raise MakefileError(message % (key, subkey))

                roidd["Name"] = subkey

                if roidd.get("SubsetRegions") is not None:
                    path = "PhylogeneticInference:%s:RegionsOfInterest:%s" % (key, subkey)
                    _collect_subsets(subkey, roidd["SubsetRegions"], path)

    if "paml:codeml" in steps:
        for (roi, subset) in mkfile["PAML"]["codeml"]["SubsetRegions"].iteritems():
            _collect_subsets(roi, subset, "PAML:codeml:SubsetRegions")


def _update_filtering(mkfile):
    samples = mkfile["Project"]["Samples"]
    groups  = mkfile["Project"]["Groups"]

    filtering = {}
    for (target, filter_by) in mkfile["Project"]["FilterSingletons"].iteritems():
        if target.startswith("<") and target.endswith(">"):
            raise MakefileError("Singleton-filtering must be specified per "
                                "sample, not by groups: %r" % (target,))
        elif target not in samples:
            raise MakefileError("Unknown/Invalid sample specifed for singleton filtering: %r" % (target,))
        elif target in filter_by:
            raise MakefileError("Attempting to filter singleton in sample using itself as comparison: %r" % (target,))

        path = "Project:FilterSingletons:%s" % (target,)
        filtering[target] = _select_samples(filter_by, groups, samples, path)

        # Implicit inclusion is allowed, since that is useful in some cases,
        # where we want to filter a sample based on the group it is a member of
        if target in filtering[target]:
            # The target itself must be excluded, as including it is invalid
            filtering[target] = filtering[target] - set((target,))
            print_warn("Warning: Sample %r is singleton-filtered using a "
                       "group it is also a member of; this may be by mistake."
                       % (target,))

        if not filtering[target]:
            raise MakefileError("No samples specified by which to "
                                "singleton-filter by for %r" % (target,))

    mkfile["Project"]["FilterSingletons"] = filtering


def _update_homozygous_contigs(mkfile):
    """Treat unspecified values for HomozygousContigs as an empty list, in
    order that the user does not need to specify "[]" for empty lists.
    """
    for regions in mkfile["Project"]["Regions"].itervalues():
        hcontigs = regions["HomozygousContigs"]

        for key, contigs in hcontigs.items():
            if contigs is None:
                hcontigs[key] = []


def _check_bam_sequences(options, mkfile, steps):
    """Check that the BAM files contains the reference sequences found in the
    FASTA file, matched by name and length; extra sequences are permitted. This
    check is only done if genotyping is to be carried out, to reduce the
    overhead of reading the BAM file headers.

    """
    if ("genotype" not in steps) and ("genotyping" not in steps):
        return

    print_info("    - Validating BAM files ...")
    bam_files = {}
    for regions in mkfile["Project"]["Regions"].itervalues():
        for sample in mkfile["Project"]["Samples"].itervalues():
            filename = os.path.join(options.samples_root, "%s.%s.bam"
                                    % (sample["Name"], regions["Prefix"]))
            if regions["Realigned"]:
                filename = add_postfix(filename, ".realigned")

            if os.path.exists(filename):
                bam_files[filename] = _collect_fasta_contigs(regions["FASTA"])

    for (filename, contigs) in bam_files.iteritems():
        with pysam.Samfile(filename) as handle:
            bam_contigs = dict(zip(handle.references, handle.lengths))

            for (contig, length) in contigs.iteritems():
                bam_length = bam_contigs.get(contig)

                if bam_length is None:
                    message = ("Reference sequence missing from BAM file; "
                               "BAM file aligned against different prefix?\n"
                               "    BAM file = %s\n    Sequence name = %s") \
                               % (filename, contig)
                    raise MakefileError(message)
                elif bam_length != length:
                    message = ("Length of reference sequence in FASTA differs "
                               "from length of sequence in BAM file; BAM file "
                               "aligned against different prefix?\n"
                               "    BAM file = %s\n"
                               "    Length in FASTA = %s\n"
                               "    Length in BAM = %s") \
                               % (filename, length, bam_length)
                    raise MakefileError(message)


def _check_sexes(mkfile):
    all_contigs = set()
    contigs_sexes = set()
    regions_sexes = set()
    for regions in mkfile["Project"]["Regions"].itervalues():
        all_contigs.update(_collect_fasta_contigs(regions["FASTA"]))

        for contigs in regions["HomozygousContigs"].itervalues():
            contigs_sexes.update(contigs)

        current_sexes = set(regions["HomozygousContigs"])
        if not regions_sexes:
            regions_sexes = current_sexes
        elif regions_sexes != current_sexes:
            raise MakefileError("List of sexes for regions %r does not "
                                "match other regions" % (regions["Name"],))

    if not regions_sexes:
        raise MakefileError("No sexes have been specified in makefile; "
                            "please list all sample sexes and assosiated "
                            "homozygous contigs (if any).")

    for sample in mkfile["Project"]["Samples"].itervalues():
        if sample.get("Sex") is None:
            if sample.get("Gender") is None:
                raise MakefileError("Please specify a sex for sample %r, or "
                                    "'NA' if not applicable."
                                    % (sample["Name"]))

            sample["Sex"] = sample.pop("Gender")
        elif sample.get("Gender") is not None:
            raise MakefileError("Both a Sex and a Gender has been specified "
                                "sample %r; the Gender field is deprecated, "
                                "please only use the Sex field."
                                % (sample["Name"]))

        if sample["Sex"] not in regions_sexes:
            sexes = ", ".join(map(repr, regions_sexes))
            message = "Sample %r has unknown sex %r; known sexes are %s" \
                % (sample["Name"], sample["Sex"], sexes)
            raise MakefileError(message)

    unknown_contigs = contigs_sexes - all_contigs
    if unknown_contigs:
        print_warn("WARNING: Unknown contig(s) in 'HomozygousContigs':\n"
                   "    - " + "\n    - ".join(unknown_contigs))
        print_warn("Please verify that the list(s) of contigs is correct!")


def _update_and_check_max_read_depth(options, mkfile):
    if any(subdd["VCF_Filter"]["MaxReadDepth"] == "auto"
           for subdd in mkfile["Genotyping"].itervalues()):
        print_info("    - Determinining max-depth from depth-histograms ...")

    for (key, settings) in mkfile["Genotyping"].iteritems():
        required_keys = set()
        for sample in mkfile["Project"]["Samples"].itervalues():
            if sample["GenotypingMethod"].lower() == "samtools":
                required_keys.add(sample["Name"])

        max_depths = settings["VCF_Filter"]["MaxReadDepth"]
        if isinstance(max_depths, types.DictType):
            # Extra keys are allowed, to make it easier
            # to temporarily disable a sample
            missing_keys = required_keys - set(max_depths)
            if missing_keys:
                missing_keys = "\n    - ".join(sorted(missing_keys))
                message = "MaxReadDepth not specified for the following " \
                          "samples for %r:\n    - %s" % (key, missing_keys)
                raise MakefileError(message)

        elif isinstance(max_depths, types.StringTypes):
            assert max_depths.lower() == "auto", max_depths
            prefix = mkfile["Project"]["Regions"][key]["Prefix"]

            settings["VCF_Filter"]["MaxReadDepth"] \
                = _read_max_depths(options, prefix, required_keys)
        else:
            max_depths = dict.fromkeys(required_keys, max_depths)
            settings["VCF_Filter"]["MaxReadDepth"] = max_depths


def _read_max_depths(options, prefix, required_keys):
    missing = []
    max_depths = {}
    for sample in required_keys:
        fname = "%s.%s.depths" % (sample, prefix)
        fpath = os.path.join(options.samples_root, fname)
        max_depths[sample] = fpath

        if not os.path.exists(fpath):
            missing.append((sample, fpath))

    if missing:
        raise MakefileError("Could not determine 'MaxReadDepth' values "
                            "automatically; .depth files are missing for one "
                            "or more samples: \n  - " +
                            "\n  - ".join("%s: %s" % item for item in missing) +
                            "\n\nEnsure that the .depth files are available, "
                            "or specify a value for 'MaxReadDepth' manually.")

    for sample, fpath in max_depths.iteritems():
        max_depths[sample] = _read_max_depth(fpath, prefix, sample)

    return max_depths


def _read_max_depth(filename, prefix, sample):
    if filename in _DEPTHS_CACHE:
        return _DEPTHS_CACHE[filename]

    max_depth = None
    max_depths = {}
    try:
        with open(filename) as handle:
            for row in parse_padded_table(handle):
                if row["Name"] != "*" and \
                        row["Sample"] == "*" and \
                        row["Library"] == "*" and \
                        row["Contig"] == "*":

                    if row["Name"] in max_depths:
                        raise MakefileError("Depth histogram %r contains "
                                            "multiple 'MaxDepth' records for "
                                            "sample %r; please rebuild!"
                                            % (filename, row["Name"]))

                    max_depths[row["Name"]] = row["MaxDepth"]
    except (OSError, IOError), error:
        raise MakefileError("Error reading depth-histogram (%s): %s"
                            % (filename, error))

    if sample in max_depths:
        max_depth = max_depths[sample]
    else:
        name_counts = {}
        name_mapping = {}
        for cand_sample, cand_max in max_depths.iteritems():
            name = cand_sample.split('.', 1)[0]
            name_mapping[name] = cand_sample
            name_counts[name] = name_counts.get(name, 0) + 1

        if name_mapping.get(sample) == 1:
            # Sample name (with some extensions) found
            # This is typical if 'paleomix depths' has been run manually.
            max_depth = max_depths[name_mapping[sample]]
        elif len(max_depths) == 1:
            # Just one sampel in the depth histogram; even though it does not
            # match, we assuem that this is the correct table. This is because
            # manually generating files / renaming files would otherwise cause
            # failure when using 'MaxDepth: auto'.
            (cand_sample, max_depth), = max_depths.items()
            print_warn("        - Name in depths file not as expected; "
                       "found %r, not %r:"
                       % (cand_sample, sample))

    if max_depth is None:
        raise MakefileError("MaxDepth for %r not found in depth-histogram: %r"
                            % (sample, filename))
    elif max_depth == "NA":
        raise MakefileError("MaxDepth is not calculated for sample %r; "
                            "cannot determine MaxDepth values automatically."
                            % (filename,))
    elif not max_depth.isdigit():
        raise MakefileError("MaxDepth is not a valid for sample %r in %r; "
                            "expected integer, found %r."
                            % (sample, filename, max_depth))

    max_depth = int(max_depth)

    print_info("        - %s.%s = %i" % (sample, prefix, max_depth))
    _DEPTHS_CACHE[filename] = max_depth
    return max_depth


_DEPTHS_CACHE = {}


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

    for (key, subdd) in genotyping.iteritems():
        if subdd.get("GenotypeEntirePrefix"):
            message = "GenotypeEntirePrefix is only allowed for prefixes " \
                      "using default parameters, but is set for %r" % (key,)
            raise MakefileError(message)

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
_VALIDATION_SAMPLES_KEY = And(IsStr, Not(_VALIDATION_SUBSAMPLE_KEY))
_VALIDATION_SAMPLES = {
    _VALIDATION_SAMPLES_KEY: {
        "GenotypingMethod": StringIn(("reference sequence",
                                      "random sampling",
                                      "samtools"),
                                     default="samtools"),
        "SpeciesName": IsStr,  # Not used; left for backwards compatibility
        "CommonName": IsStr,   # Not used; left for backwards compatibility
        "Sex": IsStr(),
        "Gender": IsStr(),
    }
}
_VALIDATION_SAMPLES[_VALIDATION_SUBSAMPLE_KEY] = _VALIDATION_SAMPLES

# Genotyping settings; note that explicit lists must not be used here, to allow
# proper inheritance of default values. Use IsListOf instead.
_VALIDATION_GENOTYPES = {
    "Padding": IsUnsignedInt,
    "GenotypeEntirePrefix": IsBoolean(default=False),
    "MPileup": {
        StringStartsWith("-"): Or(IsInt, IsStr, IsNone),
    },
    "BCFTools": {
        StringStartsWith("-"): Or(IsInt, IsStr, IsNone),
    },
    "Random": {
        "--min-distance-to-indels": IsUnsignedInt,
    },
    "VCF_Filter": {
        "MaxReadDepth": Or(IsUnsignedInt, IsDictOf(IsStr, IsUnsignedInt),
                           StringIn(("auto",))),

        "--keep-ambigious-genotypes": IsNone,
        "--min-quality": IsUnsignedInt,
        "--min-allele-frequency": IsFloat,
        "--min-mapping-quality": IsUnsignedInt,
        "--min-read-depth": IsUnsignedInt,
        "--max-read-depth": IsUnsignedInt,
        "--min-num-alt-bases": IsUnsignedInt,
        "--min-distance-to-indels": IsUnsignedInt,
        "--min-distance-between-indels": IsUnsignedInt,
        "--min-strand-bias": IsFloat,
        "--min-baseq-bias": IsFloat,
        "--min-mapq-bias": IsFloat,
        "--min-end-distance-bias": IsFloat,
    },
}

_VALIDATION_MSA = {
    "Enabled": IsBoolean(default=True),
    "Program": StringIn(("mafft",)),  # TODO: Add support for other programs

    "MAFFT": {
        "Algorithm": StringIn(("mafft", "auto",
                               "FFT-NS-1", "FFT-NS-2", "FFT-NS-i",
                               "NW-INS-i", "L-INS-i", "E-INS-i", "G-INS-i")),
        StringStartsWith("-"): CLI_PARAMETERS,
    },
}


_VALIDATION = {
    "Project": {
        "Title": IsStr(default="Untitled"),
        "Samples": _VALIDATION_SAMPLES,
        "RegionsOfInterest": {
            IsStr: {
                "Prefix": IsStr(default=REQUIRED_VALUE),
                "Realigned": IsBoolean(default=False),
                "ProteinCoding": IsBoolean(default=False),
                "IncludeIndels": IsBoolean(default=True),
                "HomozygousContigs": {
                    IsStr: Or(IsNone, IsListOf(IsStr)),

                    # The sex 'NA' defaults to no homozygous chromosomes
                    "NA": Or(IsNone, IsListOf(IsStr),
                             default=[]),
                },
            },
        },
        "FilterSingletons": {
            IsStr: [IsStr],
        },
    },
    "Genotyping": {
        "Defaults": _VALIDATION_GENOTYPES,
        IsStr: _VALIDATION_GENOTYPES,
    },
    "MultipleSequenceAlignment": {
        "Defaults": _VALIDATION_MSA,
        IsStr: _VALIDATION_MSA,
    },
    "PhylogeneticInference": {
        IsStr: {
            # Which program to use; TODO: Add support for other programs
            "Program": StringIn(("examl",), default="examl"),
            # Exclude one or more samples from the phylogeny
            "ExcludeSamples": [IsStr],
            # Which samples to root the final trees on / or midpoint rooting
            "RootTreesOn": [IsStr],
            # Create a tree per gene, for each region of interest,
            # or create a supermatrix tree from all regions specified.
            "PerGeneTrees": IsBoolean(default=False),
            # Selection of regions of interest / settings per region
            "RegionsOfInterest": {
                IsStr: {
                    "Partitions": Or(And(IsStr,
                                         ValuesSubsetOf("123456789X")),
                                     ValueIn([False]),
                                     default=REQUIRED_VALUE),
                    "SubsetRegions": Or(IsStr, IsNone, default=None),
                },
            },
            "SubsetRegions": {
                IsStr: IsStr,
            },
            "ExaML": {
                "Bootstraps": IsUnsignedInt(default=100),
                "Replicates": IsUnsignedInt(default=1),
                "Model": StringIn(("GAMMA", "PSR"),
                                  default="gamma"),
            }
        }
    },
    "PAML": {
        "codeml": {
            "ExcludeSamples": [IsStr],
            "SubsetRegions": {
                IsStr: IsStr,
            },
            IsStr: {
                "ControlFile": IsStr(default=REQUIRED_VALUE),
                "TreeFile": IsStr(default=REQUIRED_VALUE),
            },
        },
    },
}
