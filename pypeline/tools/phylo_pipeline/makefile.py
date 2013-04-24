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

import yaml

from pypeline.common.makefile import *

class MAKEFileError(RuntimeError):
    pass


def read_makefiles(filenames):
    makefiles = []
    for filename in filenames:
        makefile = read_makefile(filename, {}, _VALIDATION)
        makefile = makefile["Makefile"] # Not using extra stats
        makefile = _mangle_makefile(makefile)

        makefiles.append(makefile)

    return makefiles


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

    return mkfile




def _collapse_taxa(mkfile):
    groups = {}
    def _collect_taxa(dd, path = ()):
        current_taxa = {}
        for (key, subdd) in dd.iteritems():
            if key.startswith("<") and key.endswith(">"):
                key = key.lstrip("<").rstrip(">")
                current_taxa.update(_collect_taxa(subdd, path + (key,)))
            else:
                subdd["Name"] = key
                current_taxa[key] = subdd

        groups[path] = current_taxa
        return current_taxa

    _collect_taxa(mkfile["Project"]["Taxa"])
    mkfile["Project"]["Taxa"] = groups.pop(())
    mkfile["Project"]["Groups"] = groups


def _update_intervals(mkfile):
    intervals = mkfile["Project"]["Intervals"]
    for (interval, subdd) in intervals.iteritems():
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


def _validate_taxa(path, dd, taxa = set()):
    if not isinstance(dd, types.DictType):
        raise MAKEFileError("Expected dicts in Taxa tree for '%s', found %s: %s" \
                                % (filename, dd.__class__.__name__, dd))

    for (key, subdd) in dd.iteritems():
        if key.startswith("<") and key.endswith(">"):
            _validate_taxa(path + (key,), subdd, taxa)
        elif key.lower() not in taxa:
            taxa.add(key.lower())
            validate_makefile(subdd, _TAXA_VALIDATION, path + (key,))
        else:
            raise MAKEFileError("Taxa specified multiple times: %s" % key)


_DEFAULTS = {
    "Project" : {
        "Title" : "Untitled",
        "Taxa"  : {
            },
        "Intervals" : {
            },
        "Filter Singletons" : {
            },
        },
    "Genotyping" : {
        "Default" : "SAMTools",
        "Random" : {
            "Padding"    : 5,
            },
        "SAMTools" : {
            "Padding"    : 5,
            },
        },
    "MSAlignment" : {
        "Default" : "mafft",
        "MAFFT" : {
            "Algorithm" : "auto",
            },
        },
    "Phylogenetic Inference" : {
        "Default" : "ExaML",
        "ExaML" : {
            "Bootstraps" : 100,
            "Replicates" : 1,
            "Model"      : "gamma",
        }
    },
    "PAML" : {
        "codeml" : {
        },
    },
}


_VALIDATION = {
    "Project" : {
        "Title" : IsStr,
        "Taxa" : _validate_taxa,
        "Intervals" : {
            IsStr : {
                "Genome"         : IsStr,
                "Protein coding" : IsBoolean,
                "Orthology map"  : IsStr,
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
        "Padding"    : IsInt,
        "Default" : OneOf("random", "samtools", case_sensitive = False),
        AnyOf("MPileup", "BCFTools", "Random") : {
            IsStrWithPrefix("-") : CLI_PARAMETERS,
            },
        "VCF_Filter" : {
            "Mappability"   : IsStr,
            "MaxReadDepth"  : Or(IsUnsignedInt, IsDictOf(IsStr, IsInt)),
            IsStrWithPrefix("-") : CLI_PARAMETERS,
            },
        },
    "MSAlignment" : {
        "Default"   : AnyOf("mafft", case_sensitive = False),
        "MAFFT" : {
            "Algorithm" : AnyOf("auto", "g-ins-i", case_sensitive = False), # TODO
            },
        },
    "Phylogenetic Inference" : {
        "ExcludeGroups" : IsListOf(IsStr),
        "Default" : AnyOf("raxml", "raxml-light", "examl", case_sensitive = False),
        "ExaML" : {
            "Threads"    : IsInt,
            "Bootstraps" : IsUnsignedInt,
            "Replicates" : IsUnsignedInt,
            "Model"      : AnyOf("gamma", case_sensitive = False), # TODO
            "Outgroup"   : IsStr,
        }
    },
    "PAML" : {
        "codeml" : {
            "ExcludeGroups" : IsListOf(IsStr),
            "Control Files" : IsDictOf(IsStr, IsStr),
            "Tree File"     : IsStr,
        },
    },
}


_TAXA_VALIDATION = {
    "Genotyping Method" : OneOf("reference sequence", "random sampling", "samtools", case_sensitive = False),
    "Species Name"      : IsStr,
    "Common Name"       : IsStr,
    "Gender"            : IsStr,
    "Genomes"           : {
        IsStr             : IsStr,
        }
    }

