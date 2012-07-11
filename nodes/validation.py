from pypeline.node import Node
from pypeline.pylib.samwrap import SamfileReader
from pypeline.fileutils import reroot_path, move_file



class PairedStatisticsNode(Node):
    def __init__(self, infile, dependencies = ()):
        self._infile = infile
        self._outfile = infile + ".paired_stats"

        Node.__init__(self, 
                      description  = "<PairedStatistics: '%s' -> '%s'>" \
                          % (infile, self._outfile),
                      input_files  = [infile],
                      output_files = [infile + ".paired_stats"],
                      dependencies = dependencies)


    def _run(self, _config, temp):
        table = {}
        with SamfileReader(self._infile) as bamfile:
            for record in bamfile:
                key = (dict(record.tags)["RG"], record.tid)
    
                try:
                    subtable = table[key]
                except KeyError:
                    subtable = table[key] = [0, 0, 0]        

                if record.flag & 0x40: # first of pair
                    subtable[1] += 1
                elif record.flag & 0x80: # second of pair
                    subtable[2] += 1
                else: # Singleton
                    subtable[0] += 1

            references = bamfile.references
            lengths    = bamfile.lengths

            
        with open(reroot_path(temp, self._outfile), "w") as table_file:
            table_file.write("Group\tChr\tChrLen\tSingle\tFirst\tSecond\n")

            for ((group, tid), subtable) in sorted(table.items()):
                prefix = [group, references[tid], lengths[tid]]
                row    = [str(value) for value in (prefix + subtable)]

                table_file.write("\t".join(row) + "\n")

        move_file(reroot_path(temp, self._outfile), self._outfile)
