import os

import node

from atomiccmd import AtomicCmd



class SlopBedNode(node.SimpleNode):
    def __init__(self, config, destination, infile, outfile, genome, from_start = 0, from_end = 0, strand_relative = False, dependencies = ()):
        if type(from_start) != type(from_end):
            raise ValueError("'from_start' and 'from_end' should be of same type!")

        call = ["slopBed", 
                "-i", "%(IN_FILE)s", 
                "-g", "%(IN_GENOME)s",
                "-l", str(from_start),
                "-r", str(from_end)]
        
        if strand_relative:
            call.append("-s")            
        if type(from_start) is float:
            call.append("-pct")
        
        command = AtomicCmd(destination,
                            call,
                            IN_FILE = infile,
                            IN_GENOME = genome,
                            stdout  = outfile)

        description = "<SlopBed: '%s' -> '%s'>" \
            % (infile, os.path.join(destination, outfile))

        node.SimpleNode.__init__(self, 
                           description  = description,
                           command      = command,
                           dependencies = dependencies)
