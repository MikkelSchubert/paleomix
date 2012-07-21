#!/usr/bin/python

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


def crc32(s):
    return binascii.crc32(s) & 0xffffffff


def set_in(d, kws, value):
    """Traverses a set of nested dictionaries using the keys in
       kws, and assigns the specified value to the inner-most 
       dictionary (obtained from the second-to-last key), using
       the last key in kws. Thus calling set_in is(d, [X, Y, Z], v)
       is equivalent to calling 
         d.setdefault(X, {}).setdefault(Y, {})[Z] = v

       Behavior on non-dictionaries is undefined."""
    kws = list(kws)
    if not kws:
        raise ValueError("No keys passed to 'set_in'!")

    for kw in kws[:-1]:
        try:
            d = d[kw]
        except KeyError:
            nd = {}
            d[kw] = nd
            d = nd

    d[kws[-1]] = value


def get_in(d, kws, default = None):
    """Traverses a set of nested dictionaries using the keys in
       kws, and returns the value assigned to the final keyword
       in the innermost dictionary. Calling get_in(d, [X, Y])
       is equivalent to calling d.get(X).get(Y), with the
       difference that any missing keys causes the default value
       to be returned.

       Behavior on non-dictgionaries is undefined."""
    kws = list(kws)
    for kw in kws[:-1]:
        try:
            d = d[kw]
        except KeyError:
            return default

    return d.get(kws[-1], default)


def split_before(it, key):
    """Takes a sequence and splits it before every value where key(v) is true.
    Thus split_before(range(10), key = lambda x: x % 2 == 0) would return the
    sequence [[1], [2,3], [4,5], [6,7], [7,8], [9]]"""
    items = []
    for v in it:
        if key(v) and items:
            yield items
            items = []
        items.append(v)

    if items:
        yield items


def is_strictly_increasing(lst):
    """Returns true if the contents of the list is strictly increasing."""
    return all(x < y for (x, y) in itertools.izip(lst, itertools.islice(lst, 1, None)))


# Copied from the Python 'itertools' module documentation
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)


def group_by_pred(pred, it):
    """Splits items in a sequence into two lists, one containing
    items matching the predicate, and another containing those that
    do not."""
    is_true, is_false = [], []
    for item in it:
        if pred(item):
            is_true.append(item)
        else:
            is_false.append(item)

    return is_true, is_false
