import os

from optparse import OptionParser
import pypeline.tools.phylo_pipeline.makefile

def parse_options(argv, parser = None):
    parser = OptionParser()
    parser.add_option("--run",                default = False, action="store_true")
    parser.add_option("--non-verbose",        default = False, action="store_true")
    parser.add_option("--max-threads",        default = 12, type = int)
    parser.add_option("--temp-root",          default = "./temp")
    parser.add_option("--samples-root",       default = "./data/samples")
    parser.add_option("--intervals-root",     default = "./data/intervals")
    parser.add_option("--genomes-root",       default = "./data/genomes")
    parser.add_option("--destination",        default = "./results")
    
    (options, args) = parser.parse_args(argv)

    makefiles = pypeline.tools.phylo_pipeline.makefile.read_makefiles(args)
    if not makefiles:
        return None, None
    
    return options, makefiles


def get_genome_for_interval(interval, taxon):
    default = interval.get("Genome")
    genomes = taxon.get("Genomes", {})

    return genomes.get(interval["Name"], default)

    
def collect_bed_files(options, interval, taxa):
    bedfiles = {}
    for taxon in taxa.itervalues():
        name      = taxon["Name"]
        prefix    = get_prefix(interval, taxon)
        bedfile   = os.path.join(options.intervals_root, prefix + ".bed")
        
        bedfiles[name] = bedfile
    return bedfiles


def collect_fasta_files(options, interval, taxa):
    fastafiles = {}
    for taxon in taxa.itervalues():
        name      = taxon["Name"]
        prefix    = get_prefix(interval, taxon)
        fastafile = os.path.join(options.destination, "genotypes", "%s.%s.fasta" % (name, prefix))

        fastafiles[name] = fastafile
    return fastafiles
        

def collect_sequences(options, interval, taxa):
    bedfiles  = collect_bed_files(options, interval, taxa)
    if len(set(bedfiles.itervalues())) > 1:
        raise RuntimeError("Support for combining different intervals files not implemented!")

    # Same set of sequences for all genomes
    sequences = set()
    with open(bedfiles.itervalues().next()) as bedhandle:
        for line in bedhandle:
            sequences.add(line.split()[3])
            if len(sequences) > 10: 
                print "WARNING: REMOVE THIS! " + __file__
                break
    seqmap = dict(zip(sequences, sequences))

    return dict((name, dict.fromkeys(taxa, name)) for name in seqmap)


        
def get_prefix(interval, taxon = None):
    if not taxon:
        return "{Genome}.{Name}"
    
    genome  = get_genome_for_interval(interval, taxon)
    return "%s.%s" % (genome, interval["Name"])
