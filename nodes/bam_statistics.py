import os

import pysam

import node
import fileutils



class PairedStatisticsNode(node.Node):
    def __init__(self, config, infile, dependencies = ()):
        self._infile = infile
        self._outfile = os.path.split(infile)[-1] + ".paired_stats"
        self._destination = os.path.split(infile)[0]

        node.Node.__init__(self, 
                           description  = "<PairedStatistics: '%s' -> '%s'>" \
                               % (infile, os.path.join(self._destination, self._outfile)),
                           command      = None,
                           dependencies = dependencies)


    def output_exists(self):
        return fileutils.valid_file(os.path.join(self._destination, self._outfile))


    def run(self, config):
        table, count = {}, 0
        bamfile = pysam.Samfile(self._infile)
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

        temp = fileutils.create_temp_dir(config.temp_root)
        with open(os.path.join(temp, self._outfile), "w") as table_file:
            table_file.write("Group\tChr\tSingle\tFirst\tSecond\tRatio\tHits\n")

            references = dict(zip(bamfile.references, bamfile.lengths))
            for ((group, tid), subtable) in sorted(table.items()):
                chrom = bamfile.getrname(tid)
                total = sum(subtable)
    
                if total:
                    fractions  = ["%.2f" % (count / float(total)) for count in subtable]
                    fractions += ["%.4f" % (float(total) / references[chrom]), str(total)]
                else:
                    fractions  = ["NA"] * 3 + ["0.0000", "0"]
        
                table_file.write("\t".join([group, chrom] + fractions) + "\n")

        bamfile.close()
        os.rename(os.path.join(temp, self._outfile), 
                  os.path.join(self._destination, self._outfile))
        os.rmdir(temp)

        return True
