import os
import random

import pypeline.fileutils

from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd




class RAxMLReduceNode(CommandNode):
    def __init__(self, inalignment, outalignment, inpartitions, outpartitions, dependencies = ()):
        call = ["raxmlHPC",
                "-f", "c",
                "-m", "GTRGAMMA", # Model required, but not used
                "-n", "GTRGAMMA",
                "-w", "%(TEMP_DIR)s",
                "-s", "%(TEMP_IN_ALIGNMENT)s",
                "-q", "%(TEMP_IN_PARTITIONS)s"]

        self._kwargs = {"IN_ALIGNMENT"        : inalignment,
                        "IN_PARTITIONS"       : inpartitions,
                            
                        "TEMP_IN_ALIGNMENT"   : "RAXML_alignment",
                        "TEMP_IN_PARTITIONS"  : "RAXML_partitions",
                        "TEMP_OUT_INFO"       : "RAxML_info.GTRGAMMA",

                        "OUT_ALIGNMENT"       : outalignment,
                        "OUT_PARTITIONS"      : outpartitions}

        CommandNode.__init__(self,
                             command      = AtomicCmd(call, **self._kwargs),
                             description  = "<RAxMLReduce: '%s' -> '%s'>" \
                                     % (inalignment, outalignment),
                             dependencies = dependencies)


    def _setup(self, config, temp):
        for key in ("IN_ALIGNMENT", "IN_PARTITIONS"):
            source      = self._kwargs[key]
            destination = os.path.join(temp, self._kwargs["TEMP_" + key])

            pypeline.fileutils.copy_file(source, destination)

        CommandNode._setup(self, config, temp)


    def _teardown(self, config, temp):
        for postfix in ("ALIGNMENT", "PARTITIONS"):
            filenames = [self._kwargs["TEMP_IN_" + postfix],
                         self._kwargs["TEMP_IN_" + postfix] + ".reduced",
                         os.path.basename(self._kwargs["OUT_" + postfix])]

            for (source, destination) in zip(filenames, filenames[1:]):
                source      = os.path.join(temp, source)
                destination = os.path.join(temp, destination)

                if not os.path.exists(destination):
                    os.rename(source, destination)
                elif source != destination:
                    os.remove(source)
        
        CommandNode._teardown(self, config, temp)




class RAxMLRapidBSNode(CommandNode):
    def __init__(self, infile, partitions, destination, model = "GTRGAMMAI", threads = 1, dependencies = ()):
        call = ["raxmlHPC-PTHREADS" if (threads > 1) else "raxmlHPC",
                "-f", "a",
                "-m", model,
                "-n", model,
                "-s", "%(IN_ALIGNMENT)s",
                "-q", "%(IN_PARTITIONS)s",
                "-w", "%(TEMP_DIR)s",
                "-x", int(random.random() * 2**32),
                "-p", int(random.random() * 2**32),
                "-N", "autoMRE"]
        if threads > 1:
            call.extend(("-T", threads))


        out_tmpl = os.path.join(destination, "RAxML_%s." + model)
        command = AtomicCmd(call,
                            IN_ALIGNMENT    = infile,
                            IN_PARTITIONS   = partitions,

                            OUT_INFO        = out_tmpl % "info",
                            OUT_BESTTREE    = out_tmpl % "bestTree",
                            OUT_BOOTSTRAP   = out_tmpl % "bootstrap",
                            OUT_BIPART      = out_tmpl % "bipartitions",
                            OUT_BIPARTLABEL = out_tmpl % "bipartitionsBranchLabels")

        CommandNode.__init__(self,
                             command      = command,
                             description  = "<RAxMLRapidBS: '%s' -> '%s'>" \
                                     % (infile, destination),
                             dependencies = dependencies)

