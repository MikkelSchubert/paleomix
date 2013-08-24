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
import pypeline.common.makefile
from pypeline.common.makefile import \
     MakefileError, \
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


class MAKEFileError(RuntimeError):
    pass


def read_makefile(filename):
    makefile = pypeline.common.makefile.read_makefile(filename, _VALIDATION)
    return _mangle_makefile(makefile["Makefile"])


def _mangle_makefile(mkfile):
    _collapse_taxa(mkfile)
    _update_intervals(mkfile)
    _update_filtering(mkfile)
    mkfile["Nodes"] = ()

    padding = mkfile["Genotyping"]["Padding"]
    mkfile["Genotyping"]["Random"]["--padding"] = padding

    excluded_groups = set(mkfile["Phylogenetic Inference"]["ExcludeGroups"])
    unknown_groups  = excluded_groups - set(mkfile["Project"]["Taxa"])
    if unknown_groups:
        raise MakefileError("Unknown taxa in 'Phylogenetic Inference:ExcludeGroups': '%s'" \
                            % ("', '".join(unknown_groups)))

    for (ctl_name, ctl_files) in mkfile["PAML"]["codeml"].iteritems():
        if ctl_name not in ("ExcludeGroups",):
            if not (ctl_files.get("Control File") and ctl_files.get("Tree File")):
                raise MakefileError("Both tree and control file must be specified for PAML:codeml:%s" \
                                    % (ctl_name,))

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


def _update_intervals(mkfile):
    intervals = mkfile["Project"]["Intervals"]
    for (interval, subdd) in intervals.iteritems():
        if "Genome" not in subdd:
            raise MAKEFileError("No genome specified for interval %r" % (interval,))
        subdd["Name"] = interval


def _update_filtering(mkfile):
    taxa   = mkfile["Project"]["Taxa"]
    groups = mkfile["Project"]["Groups"]

    filtering = {}
    for (target, filter_by) in mkfile["Project"]["Filter Singletons"].iteritems():
        filtering[target] = set()
        for group in filter_by:
            if group.startswith("<") and group.endswith(">"):
                key = tuple(group[1:-1].split("/"))
                if key not in groups:
                    raise MAKEFileError("Unknown group specifed for filtering '%s': %s" % (target, key))
                filtering[target].update(groups[key])
            elif group in taxa:
                filtering[target].add(group)
            else:
                raise MAKEFileError("Unknown/Invalid group specifed for filtering '%s': '%s'" % (target, group))
    mkfile["Project"]["Filter Singletons"] = filtering


# Recursive definition of taxa tree
_VALIDATION_SUBTAXA_KEY = And(StringStartsWith("<"),
                              StringEndsWith(">"))
_VALIDATION_TAXA_KEY    = And(IsStr, Not(_VALIDATION_SUBTAXA_KEY))
_VALIDATION_TAXA = {
    _VALIDATION_TAXA_KEY : {
        "Genotyping Method" : StringIn(("reference sequence", "random sampling", "samtools")),
        "Species Name"      : IsStr,
        "Common Name"       : IsStr,
        "Gender"            : IsStr,
    }
}
_VALIDATION_TAXA[_VALIDATION_SUBTAXA_KEY] = _VALIDATION_TAXA


_VALIDATION = {
    "Project" : {
        "Title" : IsStr(default = "Untitled"),
        "Taxa" : _VALIDATION_TAXA,
        "Intervals" : {
            IsStr : {
                "Genome"         : IsStr,
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
        "Default"  : StringIn(("random", "samtools"),
                              default = "samtools"),
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
        "Default"   : StringIn(("mafft",),
                               default = "mafft"),
        "MAFFT" : {
            "Algorithm" : StringIn(("auto", "FFT-NS-1", "FFT-NS-2", "FFT-NS-i", "NW-INS-i", "L-INS-i", "E-INS-i", "G-INS-i"),
                                   default = "auto")
            },
        },
    "Phylogenetic Inference" : {
        "ExcludeGroups" : IsListOf(IsStr, default = []),
        "Default" : StringIn(("raxml", "raxml-light", "examl"),
                             default = "examl"),
        "ExaML" : {
            "Threads"    : IsUnsignedInt(default = 1),
            "Bootstraps" : IsUnsignedInt(default = 100),
            "Replicates" : IsUnsignedInt(default = 1),
            "Model"      : StringIn(("GAMMA", "PSR"),
                                    default = "gamma"),
            "Outgroup"   : IsStr,
        }
    },
    "PAML" : {
        "codeml" : {
            "ExcludeGroups" : IsListOf(IsStr),
            IsStr : {
                "Control File" : IsStr,
                "Tree File"    : IsStr,
            },
        },
    },
}
