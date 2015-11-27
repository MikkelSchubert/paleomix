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
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import re
import random

from paleomix.node import \
     Node, \
     NodeError
from paleomix.common.fileutils import \
     move_file, \
     reroot_path




class PHYLIPBootstrapNode(Node):
    """Generates a bootstrap alignment for a partition PHYLIP file;

    Note that only the PHYLIP / partitions format produced by the Node
    FastaToPartitionedInterleavedPhyNode is supported, in addition to the
    formats produced by RAxMLReduceNode.

    Parameters:
      -- input_alignment  - The input alignment file in PHYLIP format
      -- input_partition  - The input partition file in RAxML format
      -- output_alignment - The output alignment file in PHYLIP format
                            The simple (RAxML like) sequential format is used.
      -- seed             - RNG seed for selecting alignment columns."""

    def __init__(self, input_alignment, input_partition, output_alignment,
                 seed = None, dependencies = ()):
        self._input_phy  = input_alignment
        self._input_part = input_partition
        self._output_phy = output_alignment
        self._seed       = seed

        Node.__init__(self,
                      description  = "<PHYLIPBootstrap: %r -> %r>" \
                        % (input_alignment, output_alignment),
                      input_files  = (input_alignment, input_partition),
                      output_files = (output_alignment,),
                      dependencies = dependencies)


    def _run(self, _config, temp):
        if self._seed is not None:
            rng = random.Random(self._seed)
        partitions = _read_partitions(self._input_part)
        header, names, sequences = _read_sequences(self._input_phy)
        bootstraps = self._bootstrap_sequences(sequences, partitions, rng)

        temp_fpath = reroot_path(temp, self._output_phy)
        with open(temp_fpath, "w") as output_phy:
            output_phy.write(header)

            for (name, fragments) in zip(names, bootstraps):
                output_phy.write(name)
                output_phy.write(" ")
                for sequence in fragments:
                    output_phy.write(sequence)
                output_phy.write("\n")

        move_file(temp_fpath, self._output_phy)


    @classmethod
    def _bootstrap_sequences(cls, sequences, partitions, rng):
        final_partitions = [[] for _ in sequences]
        for (start, end) in partitions:
            # Convert alignment to columns, and randomly select among those
            columns = zip(*(sequence[start:end] for sequence in sequences))
            bootstrap_partition = (rng.choice(columns) for _ in columns)

            # Convert randomly selected columns back into sequences
            for (dest, partition) in zip(final_partitions,
                                         zip(*bootstrap_partition)):
                dest.append("".join(partition))

        return final_partitions



_RE_PARTITION = re.compile(r"^[A-Z]+, [^ ]+ = (\d+)-(\d+)$")
_RE_PARTITION_SINGLE = re.compile(r"^[A-Z]+, [^ ]+ = (\d+)$")

def _read_partitions(filename):
    """Read a partition file, as produced by the pipeline itself, and
    returns a list of tuples containing the (start, end) coordinates;
    each line is expected to follow the following format:

    DNA, Name = Start-End

    Multiple regions, or skips are not supported."""
    partitions = []
    with open(filename) as handle:
        for (line_num, line) in enumerate(handle):
            result = _RE_PARTITION.match(line.rstrip())
            if result:
                start, end = result.groups()
            else:
                result = _RE_PARTITION_SINGLE.match(line.rstrip())
                if not result:
                    message = ("Line %i in partitions file does not follow "
                               "expected format:\n"
                               "  Expected, either = 'DNA, Name = Start-End'\n"
                               "                or = 'DNA, Name = Start'\n"
                               "  Found = %r") % (line_num, line.rstrip())
                    raise NodeError(message)
                start, = result.groups()
                end   = start

            partitions.append((int(start) - 1, int(end)))
    return partitions


def _read_sequences(filename):
    """Collects the sequences from a PHYLIP file, and returns the header,
    the names of the sequences, and the sequences themselves. The parser
    supports interleaved sequences (as produced by the pipeline), or simple
    sequential (each paired name and sequence on a single line) as produced
    by RAxML's reduce functionality. PHYLIP files containing multiple entries
    are not supported."""
    line, header = " ", None
    with open(filename) as handle:
        # Find header
        num_sequences = num_bases = 0
        while line:
            line = handle.readline()
            if line.strip():
                header = line
                num_sequences, num_bases = map(int, line.split())
                break

        names     = [None for _ in xrange(num_sequences)]
        sequences = [[] for _ in xrange(num_sequences)]

        line_num = 0
        while line:
            line = handle.readline()
            line_strip = line.strip()
            if line_strip:
                # The first N sequences are expected to contain sample names
                index = line_num % num_sequences
                if line_num < num_sequences:
                    name, line_strip = line_strip.split(None, 1)
                    names[index] = name

                sequences[index].extend(line_strip.split())
                line_num += 1

    if len(sequences) != num_sequences:
        message = ("Expected %i sequences, but found %i in PHYLIP file:\n"
                   "    Filename = %r") % (num_sequences,
                                           len(sequences),
                                           filename)
        raise NodeError(message)

    for (index, fragments) in enumerate(sequences):
        sequences[index] = "".join(fragments)
        if len(sequences[index]) != num_bases:
            message = ("Expected %ibp sequences, found %ibp sequence for %r\n"
                       " Filename = %r") % (num_bases,
                                            len(sequences[index]),
                                            names[index],
                                            filename)
            raise NodeError(message)

    return header, names, sequences
