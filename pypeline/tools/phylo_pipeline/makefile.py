import os
import types

import yaml

class MAKEFileError(RuntimeError):
    pass


def read_makefiles(filenames):
    return [read_makefile(filename) for filename in filenames]


def read_makefile(filename):
    try:
        with open(filename) as handle:
            mkfile = yaml.safe_load(handle)
    except (OSError, IOError, yaml.scanner.ScannerError), e:
        raise MAKEFileError("Error reading '%s':\n%s" % (filename, e))

    _validate_against(None, mkfile, _MKFILE_REFERENCE)
    _collapse_taxa(mkfile)
    _update_intervals(mkfile)
    _update_filtering(mkfile)
    mkfile["Nodes"] = ()
    
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
        

def _validate_against(key, obs, reference):
    if callable(reference):
        reference(key, obs)
    elif isinstance(reference, set):
        if isinstance(obs, types.StringTypes):
            obs = obs.lower()
        
        if obs not in reference:
            raise MAKEFileError("'%s' found, expected one of %s" \
                                % (obs, ", ".join(map(str, reference))))
    elif isinstance(obs, types.DictType):
        default = reference.get("*")
        for key in obs:
            if key not in reference and not default:
                raise MAKEFileError("Unexpected key '%s' found, expected one of '%s'" \
                                        % (key, "', '".join(reference)))
            _validate_against(key, obs[key], reference.get(key, default))
    else:
        raise MAKEFileError("'%s' found, expected dict" % (obs.__class__.__name__,))
    
    
def _validate_taxa(key, dd, taxa = set()):
    if not isinstance(dd, types.DictType):
        raise MAKEFileError("Expected dicts in Taxa tree for '%s', found %s: %s" \
                                % (filename, dd.__class__.__name__, dd))
    for (key, subdd) in dd.iteritems():
        if key.startswith("<") and key.endswith(">"):
            _validate_taxa(key, subdd)
        elif key.lower() not in taxa:
            taxa.add(key.lower())
            _validate_against(key, subdd, _TAXA_REFERENCE)
        else:
            raise MAKEFileError("Taxa specified multiple times: %s" % key)
                                
                            
def is_str_list(key, value):
    if not isinstance(value, types.ListType):
        raise MAKEFileError("Expected list for key '%s', found %s: %s" \
                                % (key, value.__class__.__name__, value))
    elif not all(isinstance(field, types.StringTypes) for field in value):
        raise MAKEFileError("Expected list of strings for key '%s', found %s" \
                                % (key, value))
    

def is_positive_int(key, value):
    if not isinstance(value, types.IntType):
        raise MAKEFileError("Expected positive integer for key '%s', found %s." \
                                % (key, value.__class__.__name__))
    elif isinstance(value, types.BooleanType):
        raise MAKEFileError("Expected positive integer for key '%s', found boolean." \
                                % (key,))
    elif value < 0:    
        raise MAKEFileError("Expected positive integer for key '%s', found %i." \
                                % (key, value))

    
def is_type_of(cls):
    def _is_type_of(key, value):
        if not isinstance(value, cls):
            raise MAKEFileError("Expected %s for key '%s', found %s: %s" \
                                    % (cls.__name__, key, value.__class__.__name__, value))
    return _is_type_of

   
    
    
_MKFILE_REFERENCE = {
    "Project" : {
        "Title" : is_type_of(str),
        "Taxa" : _validate_taxa,
        "Intervals" : {
            "*" : {
                "Genome" : is_type_of(str),
                "Protein coding" : is_type_of(bool),
                "Orthology Map" : is_type_of(str),
                "Homozygous Contigs" : {
                    "*" : is_str_list,
                    },
                },
            },
        "Filter Singletons" : {
            "*" : is_str_list,
            },
        },
    "Genotyping" : {
        "Default" : set(("random", "samtools")),
        "Random" : {
            "Padding" : is_positive_int,
            "MinDistanceToIndels" : is_positive_int,
            },
        "SAMTools" : {
            "Padding"    :   is_positive_int,
            "MinDepth"   :   is_positive_int,
            "MaxDepth"   :   is_positive_int,
            "MinQuality" :   is_positive_int,
            "MinDistanceToIndels" : is_positive_int,
            },
        },
    "MSAlignment" : {
        "Default"   : set(("mafft",)),
        "MAFFT" : {
            "Algorithm" : set(("auto", "g-ins-i")), # TODO
            },
        },
    "Phylogenetic inference" : {
        "Default" : set(("raxml", "raxml-light", "examl")),
        "RAxML-Family" : {
            "Bootstraps" : is_positive_int,
            "Replicates" : is_positive_int,
            "Model": set(("gtrgamma", )), # TODO
            "Outgroup" : is_type_of(str),
        }
    }
}


_TAXA_REFERENCE = {
    "Genotyping Method" : set(("reference sequence", "random sampling", "samtools")),
    "Species Name"      : is_type_of(str),
    "Common Name"       : is_type_of(types.StringTypes),
    "Gender"            : is_type_of(str),
    "Genomes"           : {
        "*"             : is_type_of(str),
        }
    }
