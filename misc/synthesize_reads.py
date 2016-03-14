#!/usr/bin/python -3
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
from __future__ import print_function

import sys
import math
import gzip
import random

from optparse import \
    OptionParser, \
    OptionGroup

from paleomix.common.sequences import \
    reverse_complement
from paleomix.common.formats.fasta import \
    FASTA
from paleomix.common.utilities import \
    fragment
from paleomix.common.sampling import \
    weighted_sampling


def _dexp(lambda_value, position):
    return lambda_value * math.exp(-lambda_value * position)


def _rexp(lambda_value, rng):
    return - math.log(rng.random()) / lambda_value


def toint(value):
    return int(round(value))


# Adapter added to the 5' end of the forward strand (read from 5' ...)
PCR1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC%sATCTCGTATGCCGTCTTCTGCTTG"
# Adapter added to the 5' end of the reverse strand (read from 3' ...):
# rev. compl of the forward
PCR2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"


def _get_indel_length(indel_lambda, rng):
    return 1 + toint(_rexp(indel_lambda, rng))


def _get_weighted_choices(rng, sub_rate, indel_rate):
    choices_by_nt = {}
    for src_nt in "ACGT":
        choices = "ACGTID"
        probs = [sub_rate / 4] * 4     # ACGT
        probs += [indel_rate / 2] * 2  # ID
        probs[choices.index(src_nt)] = 1 - sum(probs) + sub_rate / 4
        choices_by_nt[src_nt] = weighted_sampling(choices, probs, rng)
    return choices_by_nt


def _mutate_sequence(rng, choices, refseq, indel_lambda=0):
    position = 0
    sequence, positions = [], []
    while position < len(refseq):
        ref_nt = refseq[position]
        if ref_nt not in "ACGT":
            read_nt = rng.choice("ACGT")
        else:
            read_nt = choices[ref_nt].next()

        if read_nt == "D":
            for _ in xrange(_get_indel_length(indel_lambda, rng)):
                position += 1
        elif read_nt == "I":
            for _ in xrange(_get_indel_length(indel_lambda, rng)):
                sequence.append(rng.choice("ACGT"))
                positions.append(position)
        else:
            sequence.append(read_nt)
            positions.append(position)
            position += 1
    return "".join(sequence), positions


class Specimen(object):
    """Represents a specimen, from which samples are derived.

    These are mutated by the addition of changes to the sequence

    """
    def __init__(self, options, filename):
        genome = list(FASTA.from_file(filename))
        assert len(genome) == 1, len(genome)

        self._genome = genome[0].sequence.upper()
        self._sequence = None
        self._positions = None
        self._annotations = None

        self._mutate(options)

    def _mutate(self, options):
        rng = random.Random(options.specimen_seed)
        choices = _get_weighted_choices(rng, options.specimen_sub_rate,
                                        options.specimen_indel_rate)
        self._sequence, self._positions = \
            _mutate_sequence(rng, choices, self._genome,
                             options.specimen_indel_lambda)

    @property
    def sequence(self):
        return self._sequence

    @property
    def positions(self):
        return self._positions

    @property
    def annotations(self):
        return self._annotations


class Sample(object):
    def __init__(self, options, specimen):
        self._specimen = specimen
        self._random = random.Random(options.sample_seed)
        self._options = options

        frac_endog = self._random.gauss(options.sample_endog_mu,
                                        options.sample_endog_sigma)
        self._frac_endog = min(1, max(0.01, frac_endog))
        self._endog_id = 0
        self._contam_id = 0

    def get_fragment(self):
        """Returns either a DNA fragmnet, representing either a fragment of
        the sample genome, or a randomly generated DNA sequence representing
        contaminant DNA that is not related to the species."""
        if self._random.random() <= self._frac_endog:
            return self._get_endogenous_sequence()
        return self._get_contaminant_sequence()

    def _get_contaminant_sequence(self):
        length = self._get_frag_len()
        sequence = [self._random.choice("ACGT") for _ in xrange(length)]

        self._contam_id += 1
        name = "Seq_junk_%i" % (self._contam_id,)
        return (False, name, "".join(sequence))

    def _get_endogenous_sequence(self):
        length = self._get_frag_len()
        max_position = len(self._specimen.sequence) - length
        position = self._random.randint(0, max_position)
        strand = self._random.choice(("fw", "rv"))

        sequence = self._specimen.sequence[position:position + length]
        real_pos = self._specimen.positions[position]
        if strand == "rv":
            sequence = reverse_complement("".join(sequence))

        self._endog_id += 1
        name = "Seq_%i_%i_%i_%s" % (self._endog_id, real_pos, length, strand)
        return (True, name, sequence)

    def _get_frag_len(self):
        length = toint(self._random.gauss(self._options.sample_frag_len_mu,
                                          self._options.sample_frag_len_sigma))

        return max(self._options.sample_frag_len_min,
                   min(self._options.sample_frag_len_max, length))


class Damage(object):
    def __init__(self, options, sample):
        self._options = options
        self._sample = sample
        self._random = random.Random(options.damage_seed)
        self._rates = self._calc_damage_rates(options)

    def get_fragment(self):
        is_endogenous, name, sequence = self._sample.get_fragment()
        if is_endogenous and self._options.damage:
            sequence = self._damage_sequence(sequence)
        return (name, sequence)

    def _damage_sequence(self, sequence):
        result = []
        length = len(sequence)
        for (position, nucleotide) in enumerate(sequence):
            if nucleotide == "C":
                if self._random.random() < self._rates[position]:
                    nucleotide = "T"
            elif nucleotide == "G":
                rv_position = length - position - 1
                if self._random.random() < self._rates[rv_position]:
                    nucleotide = "A"
            result.append(nucleotide)
        return "".join(result)

    @classmethod
    def _calc_damage_rates(cls, options):
        rate = options.damage_lambda
        rates = [_dexp(rate, position)
                 for position in range(options.sample_frag_len_max)]
        return rates


class Library(object):
    def __init__(self, options, sample):
        self._options = options
        self._sample = sample
        self._cache = []
        self._rng = random.Random(options.library_seed)

        self.barcode = options.library_barcode
        if self.barcode is None:
            self.barcode = "".join(self._rng.choice("ACGT") for _ in range(6))
        assert len(self.barcode) == 6, options.barcode

        pcr1 = PCR1 % (self.barcode,)
        self.lanes = self._generate_lanes(options, self._rng, sample, pcr1)

    @classmethod
    def _generate_lanes(cls, options, rng, sample, pcr1):
        lane_counts = []
        for _ in xrange(options.lanes_num):
            lane_counts.append(toint(random.gauss(options.lanes_reads_mu,
                                                  options.lanes_reads_sigma)))
        reads = cls._generate_reads(options, rng, sample,
                                    sum(lane_counts), pcr1)

        lanes = []
        for count in lane_counts:
            lanes.append(Lane(options, reads[:count]))
            reads = reads[count:]
        return lanes

    @classmethod
    def _generate_reads(cls, options, rng, sample, minimum, pcr1):
        reads = []
        while len(reads) < minimum:
            name, sequence = sample.get_fragment()
            cur_forward = sequence + pcr1
            cur_reverse = reverse_complement(sequence) + PCR2
            # Number of PCR copies -- minimum 1
            num_dupes = toint(_rexp(options.library_pcr_lambda, rng)) + 1
            for dupe_id in xrange(num_dupes):
                cur_name = "%s_%s" % (name, dupe_id)
                reads.append((cur_name, cur_forward, cur_reverse))
        random.shuffle(reads)
        return reads


class Lane(object):
    def __init__(self, options, reads):
        rng = random.Random()
        choices = _get_weighted_choices(rng, options.reads_sub_rate,
                                        options.reads_indel_rate)

        self._sequences = []
        for (name, forward, reverse) in reads:
            forward, _ = _mutate_sequence(rng, choices, forward,
                                          options.reads_indel_lambda)

            if len(forward) < options.reads_len:
                forward += "A" * (options.reads_len - len(forward))
            elif len(forward) > options.reads_len:
                forward = forward[:options.reads_len]

            reverse, _ = _mutate_sequence(rng, choices, reverse,
                                          options.reads_indel_lambda)

            if len(reverse) < options.reads_len:
                reverse += "T" * (options.reads_len - len(reverse))
            elif len(reverse) > options.reads_len:
                reverse = reverse[:options.reads_len]

            self._sequences.append((name, "".join(forward), "".join(reverse)))

    @property
    def sequences(self):
        return self._sequences


def parse_args(argv):
    parser = OptionParser()

    group = OptionGroup(parser, "Specimen")
    group.add_option("--specimen-seed", default=None,
                     help="Seed used to initialize the 'speciment', for the "
                          "creation of a random genotype. Set to a specific "
                          "values if runs are to be done for the same "
                          "genotype.")
    group.add_option("--specimen-sub-rate", default=0.005,  type=float)
    group.add_option("--specimen-indel-rate", default=0.0005, type=float)
    group.add_option("--specimen-indel-lambda", default=0.9, type=float)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Samples from specimens")
    group.add_option("--sample-seed", default=None)
    group.add_option("--sample-frag-length-mu",
                     dest="sample_frag_len_mu", default=100, type=int)
    group.add_option("--sample-frag-length-sigma",
                     dest="sample_frag_len_sigma", default=30, type=int)
    group.add_option("--sample-frag-length-min",
                     dest="sample_frag_len_min", default=0, type=int)
    group.add_option("--sample-frag-length-max",
                     dest="sample_frag_len_max", default=500, type=int)
    group.add_option("--sample-endogenous_mu",
                     dest="sample_endog_mu", default=0.75, type=float)
    group.add_option("--sample-endogenous_sigma",
                     dest="sample_endog_sigma", default=0.10, type=float)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Post mortem damage of samples")
    group.add_option("--damage", dest="damage",
                     default=False, action="store_true")
    group.add_option("--damage-seed", dest="damage_seed", default=None)
    group.add_option("--damage-lambda", dest="damage_lambda",
                     default=0.25, type=float)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Libraries from samples")
    group.add_option("--library-seed", dest="library_seed", default=None)
    group.add_option("--library-pcr-lambda", dest="library_pcr_lambda",
                     default=3, type=float)
    group.add_option("--library-barcode", dest="library_barcode", default=None)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Lanes from libraries")
    group.add_option("--lanes", dest="lanes_num", default=3, type=int)
    group.add_option("--lanes-reads-mu", dest="lanes_reads_mu",
                     default=10000, type=int)
    group.add_option("--lanes-reads-sigma", dest="lanes_reads_sigma",
                     default=2500, type=int)
    group.add_option("--lanes-reads-per-file", dest="lanes_per_file",
                     default=2500, type=int)
    parser.add_option_group(group)

    group = OptionGroup(parser, "Reads from lanes")
    group.add_option("--reads-sub-rate", dest="reads_sub_rate",
                     default=0.005, type=float)
    group.add_option("--reads-indel-rate", dest="reads_indel_rate",
                     default=0.0005, type=float)
    group.add_option("--reads-indel-lambda",
                     dest="reads_indel_lambda", default=0.9, type=float)
    group.add_option("--reads-length", dest="reads_len", default=100, type=int)
    parser.add_option_group(group)

    options, args = parser.parse_args(argv)
    if len(args) != 2:
        sys.stderr.write("Usage: %s <genome> <prefix>\n" % sys.argv[0])
        return None, None
    return options, args


def main(argv):
    options, args = parse_args(argv)
    if not options:
        return 1

    print("Generating %i lane(s) of synthetic reads ...\nDISCLAIMER: For "
          "demonstration of PALEOMIX only; the synthetic data is not "
          "biologically meaningful!" % (options.lanes_num,))

    specimen = Specimen(options, args[0])
    sample = Sample(options, specimen)
    damage = Damage(options, sample)
    library = Library(options, damage)

    for (lnum, lane) in enumerate(library.lanes, start=1):
        fragments = fragment(options.lanes_per_file, lane.sequences)
        for (readsnum, reads) in enumerate(fragments, start=1):
            templ = "%s%s_L%i_R%%s_%02i.fastq.gz" % (args[1], library.barcode,
                                                     lnum, readsnum)

            print("  Writing %s" % (templ % "{Pair}",))
            with gzip.open(templ % 1, "w") as out_1:
                with gzip.open(templ % 2, "w") as out_2:
                    for (name, seq_1, seq_2) in reads:
                        out_1.write("@%s%s/1\n%s\n" % (library.barcode, name, seq_1))
                        out_1.write("+\n%s\n" % ("I" * len(seq_1),))
                        out_2.write("@%s%s/2\n%s\n" % (library.barcode, name, seq_2))
                        out_2.write("+\n%s\n" % ("H" * len(seq_2),))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
