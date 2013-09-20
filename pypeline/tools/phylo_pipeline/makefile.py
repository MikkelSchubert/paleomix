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
import types

import pypeline.common.makefile
from pypeline.common.makefile import \
     MakefileError, \
     REQUIRED_VALUE, \
     IsInt, \
     IsStr, \
     IsDictOf, \
     IsListOf, \
     StringIn, \
     IsUnsignedInt, \
     IsBoolean, \
     ValueGE, \
     StringStartsWith, \
     StringEndsWith, \
     CLI_PARAMETERS, \
     And, \
     Or, \
     Not


def read_makefiles(filenames):
    makefiles = []
    for filename in filenames:
        makefile = pypeline.common.makefile.read_makefile(filename, _VALIDATION)
        makefile = _mangle_makefile(makefile["Makefile"])
        makefiles.append(makefile)
    return makefiles


def _mangle_makefile(mkfile):
    _collapse_taxa(mkfile)
    _update_intervals(mkfile)
    _update_filtering(mkfile)
    _update_exclusions(mkfile)
    _check_genders(mkfile)
    _check_max_read_depth(mkfile)
    mkfile["Nodes"] = ()

    padding = mkfile["Genotyping"]["Padding"]
    mkfile["Genotyping"]["Random"]["--padding"] = padding

    return mkfile


def _collapse_taxa(mkfile):
    groups, taxa = {}, set()
    def _collect_taxa(taxa_dict, path = ()):
        current_taxa = {}
        for (key, subdd) in taxa_dict.iteritems():
            if key.startswith("<") and key.endswith(">"):
                key = key.lstrip("<").rstrip(">")
                current_taxa.update(_collect_taxa(subdd, path + (key,)))
            elif key not in taxa:
                taxa.add(key)
                subdd["Name"] = key
                current_taxa[key] = subdd
            else:
                raise MakefileError("Duplicate taxa-name: %r" % (key,))

        groups[path] = current_taxa
        return current_taxa

    _collect_taxa(mkfile["Project"]["Taxa"])
    mkfile["Project"]["Taxa"] = groups.pop(())
    mkfile["Project"]["Groups"] = groups


def _select_taxa(select, groups, taxa, path):
    selection = set()
    for group in select:
        if group.startswith("<") and group.endswith(">"):
            key = tuple(group[1:-1].split("/"))
            if key not in groups:
                raise MakefileError("Unknown group specifed for filtering %r: %r" % (path, key))
            selection.update(groups[key])
        elif group in taxa:
            selection.add(group)
        else:
            raise MakefileError("Unknown/Invalid group specifed for filtering %r: %r" % (path, group))
    return selection


def _update_intervals(mkfile):
    intervals = mkfile["Project"]["Intervals"]
    for (interval, subdd) in intervals.iteritems():
        if "Genome" not in subdd:
            raise MakefileError("No genome specified for interval %r" % (interval,))
        subdd["Name"] = interval


def _update_filtering(mkfile):
    taxa   = mkfile["Project"]["Taxa"]
    groups = mkfile["Project"]["Groups"]

    filtering = {}
    for (target, filter_by) in mkfile["Project"]["Filter Singletons"].iteritems():
        if target not in taxa:
            raise MakefileError("Unknown/Invalid group specifed for singleton filtering: %r" % (target,))

        path = "Project:Filter Singletons:%s" % (target,)
        filtering[target] = _select_taxa(filter_by, groups, taxa, path)

    mkfile["Project"]["Filter Singletons"] = filtering


def _check_genders(mkfile):
    interval_genders = set()
    for interval in mkfile["Project"]["Intervals"].itervalues():
        current_genders = set(interval["Homozygous Contigs"])
        if not interval_genders:
            interval_genders = current_genders
        elif interval_genders != current_genders:
            raise MakefileError("List of genders for interval %r does not match other intervals" \
                                % (interval["Name"],))

    for taxon in mkfile["Project"]["Taxa"].itervalues():
        if taxon["Gender"] not in interval_genders:
            raise MakefileError("Taxon %r has unknown gender %r; known genders are %s" \
                                % (taxon["Name"], taxon["Gender"],
                                   ", ".join(map(repr, interval_genders))))


def _check_max_read_depth(mkfile):
    max_depths = mkfile["Genotyping"]["VCF_Filter"]["MaxReadDepth"]
    if isinstance(max_depths, types.DictType):
        required_keys = set()
        for taxon in mkfile["Project"]["Taxa"].itervalues():
            if taxon["Genotyping Method"].lower() == "samtools":
                required_keys.add(taxon["Name"])

        # Extra keys are allowed, to make it easier to temporarily disable a taxon
        missing_keys = required_keys - set(max_depths)
        if missing_keys:
            raise MakefileError("MaxReadDepth not specified for the following taxa:\n    - %s" \
                                % ("\n    - ".join(sorted(missing_keys)),))


def _update_exclusions(mkfile):
    taxa   = mkfile["Project"]["Taxa"]
    groups = mkfile["Project"]["Groups"]

    mkfile["Phylogenetic Inference"]["ExcludeGroups"] = \
      _select_taxa(mkfile["Phylogenetic Inference"]["ExcludeGroups"], groups, taxa, "Phylogenetic Inference:ExcludeGroups")

    mkfile["PAML"]["codeml"]["ExcludeGroups"] = \
      _select_taxa(mkfile["PAML"]["codeml"]["ExcludeGroups"], groups, taxa, "PAML:codeml:ExcludeGroups")


# Recursive definition of taxa tree
_VALIDATION_SUBTAXA_KEY = And(StringStartsWith("<"),
                              StringEndsWith(">"))
_VALIDATION_TAXA_KEY    = And(IsStr, Not(_VALIDATION_SUBTAXA_KEY))
_VALIDATION_TAXA = {
    _VALIDATION_TAXA_KEY : {
        "Genotyping Method" : StringIn(("reference sequence", "random sampling", "samtools"),
                                       default = "samtools"),
        "Species Name"      : IsStr,
        "Common Name"       : IsStr,
        "Gender"            : IsStr(default = REQUIRED_VALUE),
    }
}
_VALIDATION_TAXA[_VALIDATION_SUBTAXA_KEY] = _VALIDATION_TAXA


_VALIDATION = {
    "Project" : {
        "Title" : IsStr(default = "Untitled"),
        "Taxa" : _VALIDATION_TAXA,
        "Intervals" : {
            IsStr : {
                "Genome"         : IsStr(default = REQUIRED_VALUE),
                "Protein coding" : IsBoolean(default = False),
                "Include indels" : IsBoolean(default = True),
                "Homozygous Contigs" : {
                    IsStr : IsListOf(IsStr),
                    },
                },
            },
        "Filter Singletons" : {
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
            StringStartsWith("-") : CLI_PARAMETERS,
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
    "Phylogenetic Inference" : {
        "ExcludeGroups" : IsListOf(IsStr, default = []),
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
            "ExcludeGroups" : IsListOf(IsStr, default = []),
            IsStr : {
                "Control File" : IsStr(default = REQUIRED_VALUE),
                "Tree File"    : IsStr(default = REQUIRED_VALUE),
            },
        },
    },
}
