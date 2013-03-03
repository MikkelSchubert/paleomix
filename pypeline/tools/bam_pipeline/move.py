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
import optparse

import pysam

import pypeline.ui as ui


def parse_args(arg):
    fields = arg.rstrip("/").split("/")
    if len(fields) != 5:
        return None

    return { "Target" : fields[0],
             "Sample" : fields[2],
             "Library" : fields[3],
             "Lane" : fields[4] }

def collect_genomes(record):
    genomes = []
    for filename in os.listdir(record["Target"]):
        fpath = os.path.join(record["Target"], filename)
        if (filename != "reads") and os.path.isdir(fpath):
            genomes.append(filename)
    return genomes


def collect_bams(root):
    filenames = []
    for filename in os.listdir(root):
        fpath = os.path.join(root, filename)
        if os.path.isfile(fpath) and fpath.lower().endswith(".bam"):
            filenames.append(filename)
    return filenames


def move_reads(source, destination):
    tmpl = "{Target}/reads/{Sample}/{Library}/{Lane}"
    src  = tmpl.format(**source)
    dst  = tmpl.format(**destination)

    print "# Moving reads"
    print "mkdir -p '%s'" % (os.path.dirname(dst,),)
    print "mv -v '%s' '%s'" % (src, dst)
    print "rmdir --ignore-fail-on-non-empty -vp '%s'" % (os.path.dirname(src),)


def move_bams(source, destination):
    for genome in collect_genomes(source):
        tmpl = "{Target}/%s/{Sample}/{Library}/{Lane}" % genome
        src  = tmpl.format(**source)
        dst  = tmpl.format(**destination)

        print
        print "# Moving BAMs for genome '%s'" % genome
        print "mkdir -p '%s'" % (os.path.dirname(dst,),)
        print "mv -v '%s' '%s'" % (src, dst)
        print "rmdir --ignore-fail-on-non-empty -vp '%s'" % (os.path.dirname(src),)


def retag_bams(options, source, destination):
    for genome in collect_genomes(source):
        tmpl = "{Target}/%s/{Sample}/{Library}/{Lane}" % genome
        src  = tmpl.format(**source)
        dst  = tmpl.format(**destination)

        print
        for filename in collect_bams(src):
            handle = pysam.Samfile(os.path.join(src, filename))
            try:
                readgroup = handle.header.get("RG", ())
                if len(readgroup) != 1:
                    raise RuntimeError("...")
                readgroup = readgroup[0]
            finally:
                handle.close()


            readgroup["ID"] = destination["Library"]
            readgroup["SM"] = destination["Sample"]
            readgroup["LB"] = destination["Library"]
            readgroup["PU"] = destination["Lane"]

            fpath = os.path.join(dst, filename)
            picard_tmpl = "java -jar '%s/AddOrReplaceReadGroups.jar' I=%s O=%s %s" \
                % (options.jar_root, fpath, fpath + ".retagged.bam",
                   " ".join("=".join(item) for item in readgroup.iteritems()))

            print picard_tmpl
            print "mv -v '%s' '%s'" % (fpath + ".retagged.bam", fpath)

        handle.close()


def rm_files(record):
    print
    for genome in collect_genomes(record):
        print "rm -vf {Target}/{0}/{Sample}/{Library}.*".format(genome, **record)
        print "rm -vf {Target}/{0}.*".format(genome, **record)
        print "rm -vfr {Target}.{0}.mapDamage/{Library}*".format(genome, **record)

        filenames = []
        for postfix in ("", ".realigned"):
            for ext in (".bai", ".bam", ".coverage", ".validated"):
                filenames.append("{Target}.{0}{1}{2}".format(genome, postfix, ext, **record))
        filenames.append("{Target}.summary".format(**record))

        print "rm -vf %s" % (" ".join(filenames), )



def main(argv):
    parser = optparse.OptionParser()
    parser.add_option("--picard-root", default = os.path.join(os.path.expanduser('~'), "install", "picard-tools"),
                      help = "Folder containing Picard JARs (http://picard.sf.net)")
    options, args = parser.parse_args(argv)

    if len(args) != 2:
        ui.print_err("Usage: bam_pipeline move SRC DST")
        ui.print_err("  where: SRC and DST are paths in the form TARGET/reads/SAMPLE/LIBRARY/LANE")
        ui.print_err("Note that the second folder of the path (here \"reads/\") is ignored.")

        return 1

    source      = parse_args(args[0])
    destination = parse_args(args[1])

    move_reads(source, destination)
    move_bams(source, destination)
    retag_bams(options, source, destination)
    rm_files(source)
    rm_files(destination)
    print

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

