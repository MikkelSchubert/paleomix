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
import types
import itertools
import binascii


def safe_coerce_to_tuple(value):
    """Takes a value which be a single object, or an an iterable and returns the content wrapped in
    a tuple. In the case of strings, the original string object is returned in a tuple, and not as 
    a tuple of chars."""
    if isinstance(value, types.StringTypes):
        return (value,)

    try:
        return tuple(value)
    except TypeError:
        return (value,)


def try_cast(value, cast_to):
    try:
        return cast_to(value)
    except ValueError:
        return value
    

def crc32(data):
    return binascii.crc32(data) & 0xffffffff


def set_in(dictionary, keys, value):
    """Traverses a set of nested dictionaries using the given keys, 
       and assigns the specified value to the inner-most 
       dictionary (obtained from the second-to-last key), using
       the last key in keys. Thus calling set_in is(d, [X, Y, Z], v)
       is equivalent to calling 
         d.setdefault(X, {}).setdefault(Y, {})[Z] = v

       Behavior on non-dictionaries is undefined."""
    keys = list(keys)
    if not keys:
        raise ValueError("No keys passed to 'set_in'!")

    for key in keys[:-1]:
        try:
            dictionary = dictionary[key]
        except KeyError:
            new_dict = {}
            dictionary[key] = new_dict
            dictionary = new_dict

    dictionary[keys[-1]] = value


def get_in(dictionary, keys, default = None):
    """Traverses a set of nested dictionaries using the keys in
       kws, and returns the value assigned to the final keyword
       in the innermost dictionary. Calling get_in(d, [X, Y])
       is equivalent to calling d.get(X).get(Y), with the
       difference that any missing keys causes the default value
       to be returned.

       Behavior on non-dictgionaries is undefined."""
    keys = list(keys)
    for key in keys[:-1]:
        try:
            dictionary = dictionary[key]
        except KeyError:
            return default

    return dictionary.get(keys[-1], default)


def split_before(iterable, pred):
    """Takes a sequence and splits it before every value where pred(v) is true.
    Thus split_before(range(10), key = lambda x: x % 2 == 0) would return the
    sequence [[1], [2,3], [4,5], [6,7], [7,8], [9]]"""
    items = []
    for value in iterable:
        if pred(value) and items:
            yield items
            items = []
        items.append(value)

    if items:
        yield items


def is_strictly_increasing(lst):
    """Returns true if the contents of the list is strictly increasing."""
    return all(x < y for (x, y) in itertools.izip(lst, itertools.islice(lst, 1, None)))


# Copied from the Python 'itertools' module documentation
def grouper(size, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * size
    return itertools.izip_longest(fillvalue=fillvalue, *args)


def group_by_pred(pred, iterable):
    """Splits items in a sequence into two lists, one containing
    items matching the predicate, and another containing those that
    do not."""
    is_true, is_false = [], []
    for item in iterable:
        if pred(item):
            is_true.append(item)
        else:
            is_false.append(item)

    return is_true, is_false
