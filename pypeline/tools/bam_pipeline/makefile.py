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
import glob
import yaml
import types
import string
import hashlib
import datetime

import pypeline.tools.bam_pipeline.paths as paths


class MakefileError(RuntimeError):
    pass


DEFAULT_OPTIONS = {
    # Sequencing platform, used to tag read-groups.
    "Platform"       : "Illumina",
    # Offset for quality scores in FASTQ files.
    "QualityOffset" : 33,

    # BWA: Use seed during mapping
    "BWA_UseSeed"    : True,
    # BWA: Max edit-distance (int), or missing prob under 0.02 err. rate (float)
    "BWA_MaxEdit"    : 0.04,
    # BWA: Minimum quality of reads retained after mapping
    "BWA_MinQuality" : 0,
    # Contains PCR duplicates, filter if true
    "PCRDuplicates"  : True,
    # Exclude READ_TYPES from alignment/analysis
    "ExcludeReads"   : [],
}

READ_TYPES = ("Paired", "Single", "Collapsed")




def read_makefiles(filenames):
    makefiles = []
    for filename in filenames:
        makefiles.append(_read_makefile(filename))
    return makefiles


def _read_makefile(filename):
    try:
        with open(filename) as makefile:
            string = makefile.read()
            data = yaml.safe_load(string)
    except Exception, e:
        raise MakefileError(e)

    options  = _read_options(data, dict(DEFAULT_OPTIONS))
    prefixes = _read_prefixes(data)
    
    mtime = os.path.getmtime(os.path.realpath(filename))
    return {"Options"  : options,
            "Prefixes" : prefixes,
            "Targets"  : _read_structure(data, options, prefixes, 3),
            "Makefile" : {"Filename" : filename,
                          "Hash"     : hashlib.sha1(string).hexdigest(),
                          "MTime"    : datetime.datetime.fromtimestamp(mtime).strftime("%F %T.%f %z")}}


def _read_options(data, options):
    opt_obj = data.pop("Options", {})
    if not isinstance(opt_obj, dict):
        raise MakefileError("Options must be a dictionary, not a %s" % opt_obj.__class__.__name__)
    
    for (key, value) in opt_obj.iteritems():
        options[key] = value
            
    return options


def _read_prefixes(data, inherit = None):
    if "Prefixes" not in data:
        return inherit
    elif not isinstance(data["Prefixes"], list):
        raise MakefileError("Prefixes must be a list, not a %s" % data["Prefixes"].__class__.__name__)

    prefixes = {}
    for prefix in data.pop("Prefixes", []):
        if not isinstance(prefix, dict):
            raise MakefileError("Each prefix must be a dictionary, not a %s" % prefix.__class__.__name__)

        label = prefix.pop("Label", "")
        if len(prefix) > 1:
            raise MakefileError("Each prefix must contain a Label, and a path, found %s" % (prefix,))
        name, filename = prefix.items()[0]

        if name != "*":
            items = [(name, filename)]
        else:
            items = []
            for fname in glob.glob(filename):
                name = os.path.basename(fname).split(".")[0]
                items.append((name, fname))

            if not items:
                raise MakefileError("Did not find any matches for glob '%s'" % filename)                

        for (name, path) in items:
            if (name in READ_TYPES) or (name in [(s + "Reads") for s in READ_TYPES]) or (name == "Options"):
                raise MakefileError("Prefixes cannot be named '%s', please use another name." % name)
            elif (label == "*") or (set(name) & set(string.whitespace)):
                raise MakefileError("The label must not contain white-space, and cannot be '*': %s" % name)
            elif (set(name) & set(string.whitespace)):
                raise MakefileError("Prefix name must not contain whitespace:\n\t- Prefix: %s" % name)
            elif name in prefixes:
                raise MakefileError("Duplicate prefix name found: '%s'" % name)
        
            reference = paths.reference_sequence(path)
            if not reference:
                raise MakefileError("""Could not find reference sequence for prefix '%s':
       Reference sequences MUST have the extensions .fasta or .fa, and
       be located at '${prefix}', '${prefix}.fa' or '${prefix}.fasta'.""" % name)

            prefixes[name] = {"Name"      : name,
                              "Label"     : label,
                              "Path"      : path,
                              "Reference" : reference}
   
    if not prefixes:
        raise MakefileError("At least one BWA prefix must be specified in the makefile!")
    return prefixes


def _read_structure(data, options, prefixes, levels):
    options  = _read_options(data, options)

    results = {}
    for (cur_key, cur_data) in data.iteritems():
        function = _read_structure if levels else _read_barcode
        results[cur_key] = function(cur_data, options, prefixes, levels - 1)
    return results


def _read_barcode(data, options, prefixes, _levels):
    if isinstance(data, types.StringTypes):
        return _read_barcode_raw(data, options, prefixes)
    elif isinstance(data, types.DictType):
        if all((key in prefixes) for key in data):
            return _read_barcode_trimmed(data, options, prefixes)
        elif all((key in READ_TYPES) for key in data):
            return _read_barcode_aligned(data, options, prefixes)
        else:
            raise MakefileError("Error at Barcode level; keys must either be prefix-names, OR 'Paired', 'Single' or 'Collapsed'. Found: %s" \
                                    % (", ".join(data),))
    else:
        raise MakefileError("Expected string or dictionary at Barcode level, found %s" \
                                % type(data).__name__)


def _read_barcode_raw(data, options, prefixes):
    files = paths.collect_files(data)
    if not any(files.values()):
        raise MakefileError("Could not find files using search-string '%s'." % data)
    elif ("SE" not in files) and (len(files["PE_1"]) != len(files["PE_2"])):
        raise MakefileError("""Number of mate 1/2 files are not the same:"
    - Search string: %s
    - Found %i mate 1 files
    - Found %i mate 2 files""" % (data, len(files["PE_1"]), len(files["PE_2"])))
        
    return { "Reads"    : { "Raw" : files },
             "Options"  : options,
             "Prefixes" : prefixes }


def _read_barcode_trimmed(data, options, prefixes):
    bams = {}
    for (key, subdata) in data.iteritems():
        if not isinstance(subdata, types.StringTypes):
            raise MakfileError("Expected string after prefix-name at Barcode level, found %s: %s" \
                                   % (subdata.__class__.__name__, subdata))
        # External bams are treated under the assumptin that they may contain paired reads
        bams[key] = { "Paired" : subdata }
        
    return { "Reads"    : { "BAM" : bams },
             "Options"  : options,
             "Prefixes" : prefixes }


def _read_barcode_aligned(data, options, prefixes):
    result = { "Reads"    : { "Trimmed" : {} },
               "Options"  : options,
               "Prefixes" : prefixes }
    for (key, subdata) in data.iteritems():
        if not isinstance(subdata, types.StringTypes):
            raise MakfileError("Expected string after read-type at Barcode level, found %s'" \
                                   % (subdata.__class__.__name__))
        result["Reads"]["Trimmed"][key] = subdata
    return result


