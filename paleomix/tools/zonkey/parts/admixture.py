#!/usr/bin/python
#
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
#
"""

Parsing and validation of admixture results.

"""
import collections


CUTOFF = 0.001


class AdmixtureError(RuntimeError):
    pass


def read_admixture_results(filename, data, k_groups, cutoff=CUTOFF):
    key = "Group(%i)" % (k_groups,)
    names = tuple(data.sample_order) + ("-",)
    table = _admixture_read_results(filename, names)
    _admixture_validate_ancestral_groups(data, table, k_groups, cutoff)

    ancestral_groups = [[set(), value] for value in table["-"]]
    for sample, row in table.iteritems():
        if sample == '-':
            continue

        group = data.samples[sample][key]
        for index, value in enumerate(row):
            if value >= cutoff:
                ancestral_groups[index][0].add(group)

    return ancestral_groups


def get_percentiles(data, sample1, sample2, nreads, k_groups, has_ts, value):
    results = {'Sample1': sample1,
               'Sample2': sample2}

    nreads_lower = set(row['NReads'] for row in data.simulations
                       if row['NReads'] <= nreads)
    nreads_upper = set(row['NReads'] for row in data.simulations
                       if row['NReads'] >= nreads)

    if nreads_lower:
        selection = _select_simulations(data=data,
                                        sample1=sample1,
                                        sample2=sample2,
                                        nreads=max(nreads_lower),
                                        k_groups=k_groups,
                                        has_ts=has_ts)
        lower_bound, upper_bound = _get_percentile_range(selection, value)
        results['Lower'] = {'NReads': max(nreads_lower),
                            'Lower': lower_bound,
                            'Upper': upper_bound}

    if nreads_upper:
        selection = _select_simulations(data=data,
                                        sample1=sample1,
                                        sample2=sample2,
                                        nreads=min(nreads_upper),
                                        k_groups=k_groups,
                                        has_ts=has_ts)
        lower_bound, upper_bound = _get_percentile_range(selection, value)
        results['Upper'] = {'NReads': min(nreads_upper),
                            'Lower': lower_bound,
                            'Upper': upper_bound}

    return results


def _select_simulations(data, sample1, sample2, nreads, k_groups, has_ts):
    selection = []
    samples = frozenset((sample1, sample2))
    for row in data.simulations:
        if row['K'] != k_groups or row['HasTS'] != has_ts:
            continue
        elif row['NReads'] != nreads:
            continue
        elif frozenset((row['Sample1'], row['Sample2'])) != samples:
            continue

        selection.append(row)

    return selection


def _get_percentile_range(selection, value):
    selection = [(row['Percentile'], row['Value'])
                 for row in selection]
    selection.sort()

    lower_bound = 0.0
    upper_bound = 1.0

    for cur_pct, cur_value in selection:
        if cur_value > value:
            break

        lower_bound = cur_pct

    for cur_pct, cur_value in reversed(selection):
        if cur_value < value:
            break

        upper_bound = cur_pct

    return lower_bound, upper_bound


def _admixture_read_results(filename, samples):
    with open(filename) as handle:
        lines = handle.readlines()

    if len(samples) != len(lines):
        raise AdmixtureError("unexpected number of lines in admixture file; "
                             "expected %i samples, found %i"
                             % (len(samples), len(lines)))

    result = {}
    for name, line in zip(samples, lines):
        result[name] = [float(value) for value in line.split()]

    return result


def _admixture_validate_ancestral_groups(data, table, k_groups, cutoff):
    key = "Group(%i)" % (k_groups,)
    groups = collections.defaultdict(dict)
    for sample, row in table.iteritems():
        if sample not in data.samples:
            continue

        group = data.samples[sample][key]
        for index, value in enumerate(row):
            if value >= cutoff:
                groups[group][index] = True

    mixed_groups = []
    for group, memberships in sorted(groups.iteritems()):
        count = len(memberships)

        if count > 1:
            mixed_groups.append("member(s) of reference group %s assigned to "
                                "%i ancestral populations" % (group, count))

    if mixed_groups:
        raise AdmixtureError("Inconsistent ADMIXTURE results: %s; "
                             "cannot determine ancestry!"
                             % ("; ".join(mixed_groups)))
