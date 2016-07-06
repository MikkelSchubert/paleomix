#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2016 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import argparse
import datetime
import itertools
import os
import sys

import pysam

from paleomix.common.sequences import NT_CODES

import paleomix.common.fileutils as fileutils
import paleomix.tools.zonkey.common as common


_CHUNK_SIZE = 1000000


_SETTINGS_TEMPLATE = """
# Database format; is incremented when the format changes
Format: 1

# Revision number; is incremented when the database (but not format) changes
Revision: {Revision}

# Arguments passed to plink
Plink: "--horse"

# Number of autosomal chromosomes; required for e.g. PCA analyses.
# This includes autosomal chromosomes not included in the analyses.
NChroms: {NChroms}

# N bases of padding used for mitochondrial sequences; the last N bases of the
# alignments are expected to be the same as the first N bases, in order to
# allow alignments at this region of the genome, and are combined to generate
# final consensus.
MitoPadding: 30

# The minimum distance between SNPs, assuming an even distribution of SNPs
# across the genome. Used when --treemix-k is set to 'auto', which is the
# default behavior. Value from McCue 2012 (doi:10.1371/journal.pgen.1002451).
SNPDistance: 150000

"""

_BUILD_SH_TEMPLATE = """#!/bin/bash

set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes

MITO_FA="mitochondria.fasta"
if [ ! -e "${MITO_FA}" ];
then
    echo "WARNING: Mitochondrial FASTA ('${MITO_FA}') not found!"
    MITO_FA=""
fi

SIM_TXT="simulations.txt"
if [ ! -e "${SIM_TXT}" ];
then
    echo "WARNING: Simulations ('${SIM_TXT}') not found!"
    SIM_TXT=""
fi

EXAMPLES="examples"
if [ ! -d "${EXAMPLES}" ];
then
    echo "WARNING: Examples ('${EXAMPLES}') not found!"
    EXAMPLES=""
fi

FILENAME="zonkey{REVISION}.tar"
SOURCES="settings.yaml contigs.txt samples.txt ${MITO_FA} ${SIM_TXT} ${EXAMPLES} genotypes.txt build.sh"

rm -vf "${FILENAME}"

if ! tar cvf "${FILENAME}" ${SOURCES};
then
    echo "Removing partial files ..."
    rm -vf "${FILENAME}"
    exit 1
fi
"""


class ZonkeyError(RuntimeError):
    pass


def _try_cast_int(value):
    try:
        return int(value)
    except ValueError:
        return value


def _write_build_sh(args, filename):
    sys.stderr.write('Writing %r ...\n' % (filename,))
    if os.path.exists(filename) and not args.overwrite:
        sys.stderr.write('  File exists; skipping.\n')
        return

    tmpl = _BUILD_SH_TEMPLATE.replace("{REVISION}", args.revision)
    with open(filename, 'w') as handle:
        handle.write(tmpl)


def _write_genotypes(args, data, filename):
    sys.stderr.write('Writing %r ...\n' % (filename,))
    if os.path.exists(filename) and not args.overwrite:
        sys.stderr.write('  File exists; skipping.\n')
        return

    samples = data['samples']
    keys = tuple(sorted(samples))

    ref_handle = pysam.Fastafile(args.reference)

    with open(filename, 'w') as handle:
        header = ('Chrom', 'Pos', 'Ref', ';'.join(keys))
        handle.write('%s\n' % ('\t'.join(header)))

        for contig, size in sorted(data['contigs'].items()):
            # Skip non-autosomal contigs
            if not isinstance(contig, int):
                continue

            sys.stderr.write('  - %s:   0%%\r' % (contig,))
            for pos in xrange(0, size, _CHUNK_SIZE):
                sys.stderr.write('  - %s: % 3i%%\r'
                                 % (contig, (100 * pos) / size))

                chunks = []
                for key in keys:
                    real_name = samples[key]['contigs'][contig]
                    fasta_handle = samples[key]['handle']
                    chunk = fasta_handle.fetch(real_name,
                                               pos,
                                               pos + _CHUNK_SIZE)

                    chunks.append(chunk)

                ref_chunk = ref_handle.fetch(real_name, pos, pos + _CHUNK_SIZE)

                for idx, row in enumerate(itertools.izip(*chunks)):
                    if 'N' in row:
                        continue

                    nucleotides = set()
                    for nuc in row:
                        nucleotides.update(NT_CODES[nuc])

                    if len(nucleotides) == 2:
                        handle.write('%s\t%i\t%s\t%s\n'
                                     % (contig, pos + idx + 1,
                                        ref_chunk[idx], ''.join(row)))
            sys.stderr.write('  - %s: 100%%\n' % (contig,))


def _write_settings(args, contigs, filename):
    sys.stderr.write('Writing %r ...\n' % (filename,))
    if os.path.exists(filename) and not args.overwrite:
        sys.stderr.write('  File exists; skipping.\n')
        return

    # Determine the highest numbered chromosome; this is required by,
    # for example, SmartPCA.
    nchroms = max(name for name in contigs if isinstance(name, int))

    with open(filename, 'w') as handle:
        handle.write(_SETTINGS_TEMPLATE.format(Revision=args.revision,
                                               NChroms=nchroms))


def _write_contigs(args, filename):
    sys.stderr.write('Writing %r ...\n' % (filename,))
    if os.path.exists(filename) and not args.overwrite:
        sys.stderr.write('  File exists; skipping.\n')
        return

    fasta_handle = pysam.Fastafile(args.reference)
    contigs = _read_contigs(args.reference)
    lines = ['ID\tSize\tNs\tChecksum']

    for name, (real_name, size) in sorted(contigs.items()):
        sys.stderr.write('  - %s:   0%%\r' % (name,))
        n_uncalled = 0
        for pos in xrange(0, size, _CHUNK_SIZE):
            sys.stderr.write('  - %s: % 3i%%\r' % (name, (100 * pos) / size))
            chunk = fasta_handle.fetch(real_name, pos, pos + _CHUNK_SIZE)
            n_uncalled += chunk.count('n')
            n_uncalled += chunk.count('N')
            n_uncalled += chunk.count('-')

        sys.stderr.write('  - %s: 100%%\n' % (name,))
        lines.append('%s\t%i\t%i\t%s'
                     % (name, size, n_uncalled, 'NA'))
    lines.append('')

    with open(filename, 'w') as handle:
        handle.write('\n'.join(lines))


def _write_samples(args, samples, filename):
    sys.stderr.write('Writing %r ...\n' % (filename,))
    if os.path.exists(filename) and not args.overwrite:
        sys.stderr.write('  File exists; skipping.\n')
        return

    lines = ['ID\tGroup(2)\tGroup(3)\tSpecies\tSex\tSampleID\tPublication']

    for name in sorted(samples):
        lines.append('%s\t-\t-\tNA\tNA\t%s\tNA'
                     % (name, name))
    lines.append('')

    with open(filename, 'w') as handle:
        handle.write('\n'.join(lines))


def _process_contigs(reference, samples):
    ref_contigs = _read_contigs(reference)
    for name, (_, size) in ref_contigs.items():
        ref_contigs[name] = size

    for sample_name, obs_data in samples.items():
        obs_contigs = obs_data['contigs']
        for ref_name, ref_size in ref_contigs.iteritems():
            if ref_name not in obs_contigs:
                raise ZonkeyError('Contig missing for sample %r: %r'
                                  % (sample_name, ref_name))

            obs_name, obs_size = obs_contigs[ref_name]
            if obs_size != ref_size:
                raise ZonkeyError('Contig %r for sample %r has wrong size; '
                                  '%i observed vs %i expected'
                                  % (obs_name, sample_name,
                                     obs_size, ref_size))

            obs_contigs[ref_name] = obs_name

    return {'samples': samples,
            'contigs': ref_contigs}


def _read_contigs(filename):
    contigs = {}
    with open(filename + '.fai') as handle:
        for line in handle:
            name, size, _ = line.split('\t', 2)
            if name in contigs:
                raise ZonkeyError('FASTA file %r contains multiple contigs '
                                  'with the same name (%r); this is not '
                                  'supported.' % (filename, name))

            fixed_name = common.contig_name_to_plink_name(name)
            if fixed_name is not None:
                contigs[_try_cast_int(fixed_name)] = (name, int(size))

    return contigs


def _collect_samples(reference, filenames):
    samples = {}
    for filename in filenames:
        basename = os.path.basename(filename).split('.', 1)[0]
        if basename in samples:
            raise ZonkeyError('Duplicate sample name %r'
                              % (filename,))

        # Open first to insure that file is indexed
        handle = pysam.Fastafile(filename)
        contigs = _read_contigs(filename)
        if not contigs:
            raise ZonkeyError('No usable contigs found in %r.'
                              % (filename,))

        samples[basename] = {'handle': handle,
                             'contigs': contigs}

    return _process_contigs(reference, samples)


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    parser.add_argument('reference')
    parser.add_argument('FASTAs', nargs="+")
    parser.add_argument('--overwrite', default=False,
                        action='store_true')

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    args.revision = datetime.datetime.today().strftime('%Y%m%d')

    data = _collect_samples(args.reference, args.FASTAs)
    if not data:
        return 1

    fileutils.make_dirs(args.root)

    _write_contigs(args, os.path.join(args.root, 'contigs.txt'))
    _write_samples(args, data['samples'], os.path.join(args.root, 'samples.txt'))
    _write_settings(args, data['contigs'], os.path.join(args.root, 'settings.yaml'))
    _write_genotypes(args, data, os.path.join(args.root, 'genotypes.txt'))
    _write_build_sh(args, os.path.join(args.root, 'build.sh'))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


