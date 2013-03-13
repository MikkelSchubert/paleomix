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
from pypeline.nodes.samtools import SAMTOOLS_VERSION

import pypeline.common.versions as versions


BWA_VERSION = versions.Requirement(call   = ("bwa",),
                                   search = r"Version: (\d+)\.(\d+)\.(\d+)",
                                   checks = versions.And(versions.GE(0, 5, 9),
                                                         versions.LT(0, 6, 0)))

# Required by safeSam2Bam for 'PG' tagging support / known good version
# Cannot be a lambda due to need to be able to pickle function
def _get_pysam_version():
    return __import__("pysam").__version__
PYSAM_VERSION = versions.Requirement(name   = "module 'pysam'",
                                     call   = _get_pysam_version,
                                     search = b"(\d+)\.(\d+)\.(\d+)",
                                     checks = versions.GE(0, 7, 4))

class BWAIndexNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_file, prefix = None, dependencies = ()):
        prefix = prefix if prefix else input_file
        params = _BWAParams(("bwa", "index"), prefix, iotype = "OUT",
                            IN_FILE = input_file,
                            TEMP_OUT_PREFIX = os.path.basename(prefix),
                            CHECK_BWA = BWA_VERSION)

        # Input fasta sequence
        params.push_positional("%(IN_FILE)s")
        # Destination prefix, in temp folder
        params.set_parameter("-p", "%(TEMP_OUT_PREFIX)s")

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


def BWANode(input_file_1, input_file_2, **kwargs):
    if input_file_1 and input_file_2:
        return PE_BWANode(input_file_1 = input_file_1,
                          input_file_2 = input_file_2,
                          **kwargs)
    elif input_file_1:
        return SE_BWANode(input_file = input_file_1,
                          **kwargs)


class SE_BWANode(CommandNode):
    @create_customizable_cli_parameters
    def customize(self, input_file, output_file, reference, prefix, threads = 1, dependencies = ()):
        threads = _get_max_threads(reference, threads)

        aln = _BWAParams(("bwa", "aln"), prefix,
                         IN_FILE = input_file,
                         OUT_STDOUT = AtomicCmd.PIPE,
                         CHECK_BWA = BWA_VERSION)
        aln.push_positional(prefix)
        aln.push_positional("%(IN_FILE)s")
        aln.set_parameter("-t", threads)

        samse = _BWAParams(("bwa", "samse"), prefix,
                           IN_STDIN = aln,
                           IN_FILE  = input_file,
                           OUT_STDOUT = AtomicCmd.PIPE,
                           CHECK_BWA = BWA_VERSION)
        samse.push_positional(prefix)
        samse.push_positional("-")
        samse.push_positional("%(IN_FILE)s")

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
            aln = _BWAParams(("bwa", "aln"), prefix,
                             IN_FILE = filename,
                             OUT_STDOUT = AtomicCmd.PIPE,
                             TEMP_OUT_SAI = "pair_%i.sai" % iindex,
                             CHECK_BWA = BWA_VERSION)
            aln.push_positional(prefix)
            aln.push_positional("%(IN_FILE)s")
            aln.set_parameter("-f", "%(TEMP_OUT_SAI)s")
            aln.set_parameter("-t", max(1, threads / 2))
            alns.append(aln)
        aln_1, aln_2 = alns

        sampe = _BWAParams(("bwa", "sampe"), prefix,
                           IN_FILE_1     = input_file_1,
                           IN_FILE_2     = input_file_2,
                           TEMP_IN_SAI_1 = "pair_1.sai",
                           TEMP_IN_SAI_2 = "pair_2.sai",
                           OUT_STDOUT    = AtomicCmd.PIPE,
                           CHECK_BWA = BWA_VERSION)
        sampe.push_positional(prefix)
        sampe.push_positional("%(TEMP_IN_SAI_1)s")
        sampe.push_positional("%(TEMP_IN_SAI_2)s")
        sampe.push_positional("%(IN_FILE_1)s")
        sampe.push_positional("%(IN_FILE_2)s")
        sampe.set_parameter("-P", fixed = False)

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

        aln = _BWAParams(("bwa", "bwasw"), prefix,
                         IN_FILE_1  = input_file_1,
                         OUT_STDOUT = AtomicCmd.PIPE)
        aln.set_parameter(prefix)
        aln.set_parameter("%(IN_FILE_1)s")

        if input_file_2:
            aln.set_parameter("%(IN_FILE_2)s")
            aln.set_paths("IN_FILE_2", input_file_2)

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
    convert.set_paths(IN_STDIN    = stdin,
                      OUT_STDOUT  = AtomicCmd.PIPE,
                      CHECK_PYSAM = PYSAM_VERSION)

    flt = AtomicParams(["samtools", "view"])
    flt.push_positional("-")
    flt.set_parameter("-b") # Output BAM
    flt.set_parameter("-u") # Output uncompressed BAM
    flt.set_parameter("-F", "0x4", sep = "", fixed = False) # Remove misses
    flt.set_paths(IN_STDIN  = convert,
                  OUT_STDOUT = AtomicCmd.PIPE,
                  CHECK_SAM = SAMTOOLS_VERSION)

    sort = AtomicParams(("samtools", "sort"))
    sort.set_parameter("-o") # Output to STDOUT on completion
    sort.push_positional("-")
    sort.push_positional("%(TEMP_OUT_BAM)s")
    sort.set_paths(IN_STDIN     = flt,
                   OUT_STDOUT   = AtomicCmd.PIPE,
                   TEMP_OUT_BAM = "sorted",
                   CHECK_SAM = SAMTOOLS_VERSION)

    calmd = AtomicParams(("samtools", "calmd"))
    calmd.push_positional("-")
    calmd.push_positional("%(IN_REF)s")
    calmd.set_parameter("-b") # Output BAM
    calmd.set_paths(IN_REF   = reference,
                    IN_STDIN = sort,
                    OUT_STDOUT = output_file,
                    CHECK_SAM = SAMTOOLS_VERSION)

    order = ["convert", "filter", "sort", "calmd"]
    dd = {"convert" : convert,
          "filter"  : flt,
          "sort"    : sort,
          "calmd"   : calmd}

    return order, dd



def _BWAParams(call, prefix, iotype = "IN", **kwargs):
    params = AtomicParams(call, **kwargs)
    for postfix in ("amb", "ann", "bwt", "pac", "rbwt", "rpac", "rsa", "sa"):
        key = "%s_PREFIX_%s" % (iotype, postfix.upper())
        params.set_paths(key, prefix + "." + postfix)

    return params


def _get_max_threads(reference, threads):
    if not os.path.exists(reference):
        return threads
    elif os.path.getsize(reference) < 2 ** 20: # 1MB
        return 1

    return threads
