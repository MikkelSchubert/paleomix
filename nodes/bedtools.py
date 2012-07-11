from pypeline.node import CommandNode
from pypeline.atomiccmd import AtomicCmd



class SlopBedNode(CommandNode):
    def __init__(self, infile, outfile, genome, from_start = 0, from_end = 0, strand_relative = False, dependencies = ()):
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
        
        command = AtomicCmd(call,
                            IN_FILE = infile,
                            IN_GENOME = genome,
                            stdout  = outfile)

        description = "<SlopBed: '%s' -> '%s'>" % (infile, outfile)

        CommandNode.__init__(self, 
                             description  = description,
                             command      = command,
                             dependencies = dependencies)
