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
import pysam


class SamfileReader(pysam.Samfile):
    def __init__(self, filename, **kws):
        kws["mode"] = None # Forces auto-detection of format, read-only
        pysam.Samfile.__init__(self, filename, **kws)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()


class FastafileReader:
    """Wrapper class for pysam.Fastafile. Exists due to a bug in pysam.Fastafile
       or samtools (?), which causes fetch to fail (apparently randomly) if
       sequences are read out of order. Additionally, adds the ability to fetch a
       list of sequences found in an indexed fasta file."""

    def __init__(self, filename):
        self._filename = filename
        
        # This will create the .fai file if missing ...
        # TODO: Catch IOError and check if this is due to lack of +w for the folder,
        #       which is required to create a .fai.
        self._fastafile = pysam.Fastafile(filename)

        fpath = self._filename + ".fai"
        with open(fpath, "r") as f:
            # Columns are <name>, <length>, <start-offset>, \
            #             <line-length>, <line-length-with-whitespace>
            self._sequences = [l.split("\t", 2)[0] for l in f]

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
        
    def get_sequences(self):
        """Given a path to a FASTA or .fai file, returns tuple for each sequence
           in the index, containing the sequence name and length."""
        return list(self._sequences)
    
    def fetch(self, contig):
        """."""
        if contig not in self._sequences:
            raise KeyError("Contig is not found in FASTA file: '%s'" % contig)
        sequence = self._fastafile.fetch(contig)
        assert sequence
        return sequence

    def close(self):
        return self._fastafile.close()


class Tabixfile(pysam.Tabixfile):
    def __new__(cls, filename, mode = "r"):
        return pysam.Tabixfile.__new__(cls, filename, mode)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()


def read_tabix(lines, parser):
    """Reads a tabix indexed file, and return a dictionary of contigs."""
    table = {}
    for line in lines:
        if line.startswith("#"):
            continue

        line = line.rstrip()
        record = parser(line, len(line))
        contig = record.contig
        
        try:
            table[contig].append(record)
        except KeyError:
            table[contig] = [record]

    return table

def read_tabix_VCF(lines):
    return read_tabix(lines, pysam.asVCF())

def read_tabix_BED(lines):
    return read_tabix(lines, pysam.asBed())

def read_tabix_GTF(lines):
    parser = pysam.asGTF()
    for line in lines:
        yield parser(line, len(line))
