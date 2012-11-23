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

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd
from pypeline.atomicparams import *
from pypeline.atomicset import ParallelCmds


class BWAIndexNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file, prefix = None, dependencies = ()):
        params = AtomicParams("bwa")
        params.set_parameter("index")
        params.set_parameter("%(IN_FILE)s")

        # Destination prefix, in temp folder
        params.set_parameter("-p", "%(TEMP_OUT_PREFIX)s")

        prefix = prefix if prefix else input_file
        params.set_paths(IN_FILE = input_file,
                         TEMP_OUT_PREFIX = os.path.basename(prefix),
                         **_prefix_files(prefix, iotype = "OUT"))
        
        return {"prefix":  prefix,
                "command": params}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = parameters.command.finalize()
        description =  "<BWA Index '%s' -> '%s.*'>" % (parameters.input_file, 
                                                       parameters.prefix)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             dependencies = parameters.dependencies)


class SE_BWANode(CommandNode):
    @create_customizable_cli_parameters
    def customize(self, input_file, output_file, reference, prefix, threads = 1, dependencies = ()):
        threads = _get_max_threads(reference, threads)

        aln   = AtomicParams("bwa")
        aln.set_parameter("aln")
        aln.set_parameter(prefix)
        aln.set_parameter("%(IN_FILE)s")
        aln.set_parameter("-t", threads)
        aln.set_paths(IN_FILE = input_file,
                      OUT_STDOUT = AtomicCmd.PIPE,
                      **_prefix_files(prefix))
        
        samse = AtomicParams("bwa")
        samse.set_parameter("samse")
        samse.set_parameter(prefix)
        samse.set_parameter("-")
        samse.set_parameter("%(IN_FILE)s")
        samse.set_paths(IN_STDIN = aln,
                        IN_FILE  = input_file,
                        OUT_STDOUT = AtomicCmd.PIPE)

        order, commands = _process_output(samse, output_file, reference)
        commands["samse"] = samse
        commands["aln"]   = aln

        return {"commands" : commands,
                "order"    : ["aln", "samse"] + order,
                "threads"  : threads}

        
    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = ParallelCmds([parameters.commands[key].finalize() for key in parameters.order])
        description =  "<SE_BWA (%i threads): '%s'>" % (parameters.threads, parameters.input_file)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             threads      = parameters.threads,
                             dependencies = parameters.dependencies)



class PE_BWANode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file_1, input_file_2, output_file, reference, prefix, threads = 2, dependencies = ()):
        threads = _get_max_threads(reference, threads)

        alns = []
        for (iindex, filename) in enumerate((input_file_1, input_file_2), start = 1):
            aln = AtomicParams("bwa")
            aln.set_parameter("aln")
            aln.set_parameter(prefix)
            aln.set_parameter("%(IN_FILE)s")
            aln.set_parameter("-f", "%(TEMP_OUT_SAI)s")
            aln.set_parameter("-t", max(1, threads / 2))
            aln.set_paths(IN_FILE = filename,
                          OUT_STDOUT = AtomicCmd.PIPE,
                          TEMP_OUT_SAI = "pair_%i.sai" % iindex,
                          **_prefix_files(prefix))
            alns.append(aln)
        aln_1, aln_2 = alns
        
        sampe = AtomicParams("bwa")
        sampe.set_parameter("sampe")
        sampe.set_parameter(prefix)
        sampe.set_parameter("%(TEMP_IN_SAI_1)s")
        sampe.set_parameter("%(TEMP_IN_SAI_2)s")
        sampe.set_parameter("%(IN_FILE_1)s")
        sampe.set_parameter("%(IN_FILE_2)s")
        sampe.set_parameter("-P", fixed = False)
        sampe.set_paths(IN_FILE_1     = input_file_1,
                        IN_FILE_2     = input_file_2,
                        TEMP_IN_SAI_1 = "pair_1.sai",
                        TEMP_IN_SAI_2 = "pair_2.sai",
                        OUT_STDOUT    = AtomicCmd.PIPE)

        order, commands = _process_output(sampe, output_file, reference)
        commands["sampe"] = sampe
        commands["aln_1"] = aln_1
        commands["aln_2"] = aln_2

        return {"commands" : commands,
                "order"    : ["aln_1", "aln_2", "sampe"] + order,
                # At least one thread per 'aln' process
                "threads"  : max(2, threads)}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        command = ParallelCmds([parameters.commands[key].finalize() for key in parameters.order])
        description =  "<PE_BWA (%i threads): '%s'>" % (parameters.threads, parameters.input_file_1)
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             threads      = parameters.threads,
                             dependencies = parameters.dependencies)


    def _setup(self, _config, temp):
        os.mkfifo(os.path.join(temp, "pair_1.sai"))
        os.mkfifo(os.path.join(temp, "pair_2.sai"))



class BWASWNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file_1, output_file, reference, prefix, input_file_2 = None, threads = 1, dependencies = ()):
        threads = _get_max_threads(reference, threads)

        aln = AtomicParams("bwa")
        aln.set_parameter("bwasw")
        aln.set_parameter(prefix)
        aln.set_parameter("%(IN_FILE_1)s")

        files   = {"IN_FILE_1"  : input_file_1,
                   "OUT_STDOUT" : AtomicCmd.PIPE}
        if input_file_2:
            aln.set_parameter("%(IN_FILE_2)s")
            files["IN_FILE_2"] = input_file_2
        files.update(_prefix_files(prefix))
        aln.set_paths(**files)

        aln.set_parameter("-t", threads)


        order, commands = _process_output(aln, output_file, reference)
        commands["aln"] = aln
        
        return {"commands" : commands,
                "order"    : ["aln"] + order,
                "threads"  : threads}


    @use_customizable_cli_parameters
    def __init__(self, parameters):
        if parameters.input_file_2:
            description =  "<PE_BWASW (%i threads): '%s', '%s' -> '%s'>" \
                % (parameters.threads, parameters.input_file_1, parameters.input_file_2, parameters.output_file)
        else:
            description =  "<BWASW (%i threads): '%s' -> '%s'>" \
                % (parameters.threads, parameters.input_file_1, parameters.output_file)

        command = ParallelCmds([parameters.commands[key].finalize() for key in parameters.order])
        CommandNode.__init__(self, 
                             command      = command,
                             description  = description,
                             threads      = parameters.threads,
                             dependencies = parameters.dependencies)


def _process_output(stdin, output_file, reference):
    convert = AtomicParams("safeSAM2BAM")
    convert.set_parameter("--flag-as-sorted")
    convert.set_paths(IN_STDIN   = stdin,
                      OUT_STDOUT = AtomicCmd.PIPE)

    flt = AtomicParams("samtools")
    flt.set_parameter("view")
    flt.set_parameter("-")
    flt.set_parameter("-b") # Output BAM
    flt.set_parameter("-u") # Output uncompressed BAM
    flt.set_parameter("-F", "0x4", sep = "", fixed = False) # Remove misses
    flt.set_paths(IN_STDIN  = convert,
                  OUT_STDOUT = AtomicCmd.PIPE)

    sort = AtomicParams("samtools")
    sort.set_parameter("sort")
    sort.set_parameter("-o") # Output to STDOUT on completion
    sort.set_parameter("-")
    sort.set_parameter("%(TEMP_OUT_BAM)s")
    sort.set_paths(IN_STDIN     = flt,
                   OUT_STDOUT   = AtomicCmd.PIPE,
                   TEMP_OUT_BAM = "sorted")

    calmd = AtomicParams("samtools")
    calmd.set_parameter("calmd")
    calmd.set_parameter("-b") # Output BAM
    calmd.set_parameter("-")
    calmd.set_parameter("%(IN_REF)s")
    calmd.set_paths(IN_REF   = reference,
                    IN_STDIN = sort,
                    OUT_STDOUT = output_file)

    order = ["convert", "filter", "sort", "calmd"]
    dd = {"convert" : convert,
          "filter"  : flt,
          "sort"    : sort,
          "calmd"   : calmd}

    return order, dd



def _prefix_files(prefix, iotype = "IN"):
    files = {}
    for postfix in ("amb", "ann", "bwt", "pac", "rbwt", "rpac", "rsa", "sa"):
        files["%s_PREFIX_%s" % (iotype, postfix.upper())] = prefix + "." + postfix
    return files


def _get_max_threads(reference, threads):
    if not os.path.exists(reference):
        return threads
    elif os.path.getsize(reference) < 2 ** 20: # 1MB
        return 1

    return threads
