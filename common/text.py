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
import string
import operator

from utilities import try_cast



def padded_table(table, min_padding = 4):
    sizes = [0] * len(table[0])
    for row in table:
        for (ii, field) in enumerate(row):
            sizes[ii] = max(sizes[ii], len(str(field)))
    sizes = [(size + min_padding) for size in sizes]

    for row in table:
        padded = []
        for (ii, field) in enumerate(row):
            padded.append(string.ljust(str(field), sizes[ii]))
        yield "".join(padded).rstrip()


def parse_padded_table(lines, sep = None, header = None):
    def _do_split(line):
        fields = line.split(sep)
        fields[-1] = fields[-1].strip()
        return fields

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        elif header is None:
            header = _do_split(line)
            continue

        yield dict(zip(header, _do_split(line)))

    
def merge_table_records(tables, sum_columns = ()):
    if not any(tables):
        yield {}

    filtered = filter(None, tables)
    fields = filtered[0][0].keys()
    keys = [field for field in fields if field not in sum_columns]

    if set(sum_columns) - set(fields):
        raise RuntimeError("Unexpected column-names in 'sum_columns'")    

    merged = {}
    for table in filtered:
        for row in table:
            key = tuple(row[field] for field in keys)
            values = [int(row[field]) for field in sum_columns]
            
            if key not in merged:
                merged[key] = values
            elif sum_columns:
                merged[key] = map(operator.add, merged[key], values)
    
    for (key, values) in merged.iteritems():
        record = {}
        record.update(zip(keys, key))
        record.update(zip(sum_columns, values))

        yield record


def merge_table_files(filenames, sum_columns = (), sep = None):
    tables = []
    for filename in filenames:
        with open(filename) as table_file:
            lines = table_file.readlines()
            tables.append(list(parse_padded_table(lines, sep = sep)))
    
    header = lines[0].split(sep)
    header[-1] = header[-1].strip()

    rows = []
    for row in merge_table_records(tables, sum_columns = ("Hits", "NTs")):
        rows.append([try_cast(row[key], int) for key in header])
    rows.sort()
    rows.insert(0, header)
    
    for row in padded_table(rows):
        yield row
