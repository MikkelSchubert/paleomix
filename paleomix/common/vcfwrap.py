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
"""Wrapper and utility functions for VCF handling, using
the VCF data-structures from pysam."""

import re
import collections


_re_tmpl = "(^|;)%s=([^;]+)(;|$)"
_re_cache = {}


def get_info(vcf, field, default = None, type = str):
    """Returns the value of the specified field from the info column
    of a VCF record. The type is determined by 'type' parameter, which
    may be any function. If no matching key is found, or if the key is
    not associated with a value, the function returns None by default."""
    try:
        regexp = _re_cache[field]
    except KeyError:
        regexp = _re_cache[field] = re.compile(_re_tmpl % (field,))

    match = regexp.search(vcf.info)
    if not match:
        return default

    return type(match.groups()[1])


Indel = collections.namedtuple("Indel", ["in_reference", "pos", "prefix", "what", "postfix"])

def parse_indel(vcf):
    """Parses the VCF record of an indel, and returns a tuple containing the
    position (0-based) of the previous base, a boolean indicating whether or
    not the subsequent sequence is found in the reference sequence, and a 
    string containing the bases added to / removed from the reference.

    Thus (7, False, "ACGT", "AC") indicates that the sequence ACGT has been 
    inserted following the 8th nucleotide, compared with the reference, and
    that the insertion is followed by the bases "AC" on the reference."""
    if not is_indel(vcf):
        raise ValueError("SNP passed to 'parse_indel'!")
    elif "," in vcf.alt:
        raise ValueError("VCF records with multiple indels not supported!")
    elif vcf.ref[0] != vcf.alt[0]:
        raise ValueError("Sequences do not match VCF spec, first base differs: "
                         "%s:%s -- %s > %s" % (vcf.contig, vcf.pos + 1, vcf.ref, vcf.alt))

    ref_len, alt_len = len(vcf.ref), len(vcf.alt)
    # The length of the insertion / deletion
    len_diff = abs(alt_len - ref_len)

    # Wheter or not the sequence 'what' is found in the reference
    in_reference = (ref_len >= alt_len)

    # The sequence added or removed from the reference
    longest = max(vcf.ref, vcf.alt, key = len)
    shortest = min(vcf.ref, vcf.alt, key = len)
    what = longest[1:len_diff + 1]

    postfix = shortest[1:]
    if longest[len_diff + 1:] != postfix:
        raise ValueError("Sequence postfix does not match; malformed indel!")

    return Indel(in_reference, vcf.pos, vcf.ref[0], what, postfix)


def is_indel(vcf):
    """Returns true if the VCF entry represents an indel."""
    # FIXME: Is this a universal key for indels?
    return "INDEL" in vcf.info


def get_genotype(vcf, sample=0, _re=re.compile(r'[|/]')):
    """Returns the most likely genotype of a sample in a vcf record. If no
    single most likely genotype can be determined, the function returns 'N' for
    both bases."""
    nucleotides = []
    nucleotides.extend(vcf.ref.split(","))
    nucleotides.extend(vcf.alt.split(","))

    result = []
    for genotype in _re.split(get_format(vcf, sample)["GT"]):
        result.append(nucleotides[int(genotype)])

    return result


# The corresponding nucleotides for each value in the VCF PL field
_genotype_indices = [(jj, ii)
                     for ii in range(0, 10)
                     for jj in range(0, ii + 1)]


def get_ml_genotype(vcf, sample=0):
    """Returns the most likely genotype of a sample in a vcf record. If no
    single most likely genotype can be determined, the function returns 'N' for
    both bases."""
    genotypes = []
    genotypes.extend(vcf.ref.split(","))
    genotypes.extend(vcf.alt.split(","))

    PL = map(int, get_format(vcf, sample)["PL"].split(","))

    expected_length = (len(genotypes) * (len(genotypes) + 1)) // 2
    if len(PL) != expected_length:
        raise ValueError("Expected %i PL values, found %i"
                         % (expected_length, len(PL)))

    if PL.count(min(PL)) > 1:
        # No single most likely genotype
        return ("N", "N")

    most_likely = min(xrange(len(PL)), key=PL.__getitem__)
    prefix, postfix = _genotype_indices[most_likely]

    return (genotypes[prefix], genotypes[postfix])


def get_format(vcf, sample=0):
    return dict(zip(vcf.format.split(":"),
                    vcf[sample].split(":")))
