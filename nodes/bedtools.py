import os

import node

from atomiccmd import AtomicCmd



class SlopBedNode(node.Node):
    def __init__(self, config, destination, infile, outfile, genome, amount = 0, dependencies = ()):
        command = AtomicCmd(destination,
                            ["slopBed", "-b", str(amount), 
                             "-i", "%(IN_FILE)s",
                             "-g", "%(IN_GENOME)s"],
                            IN_FILE = infile,
                            IN_GENOME = genome,
                            stdout  = outfile)

        description = "<SlopBed: '%s' -> '%s'>" \
            % (infile, os.path.join(destination, outfile))

        node.Node.__init__(self, 
                           description  = description,
                           command      = command,
                           dependencies = dependencies)
