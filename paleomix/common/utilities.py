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
import binascii
import copy
import cPickle
import heapq
import itertools
import pickle
import types


def _safe_coerce(cls):
    def _do_safe_coerce(value):
        if isinstance(value, (types.StringTypes, types.DictType)):
            return cls((value,))

        try:
            return cls(value)
        except TypeError:
            return cls((value,))

    _do_safe_coerce.__doc__ = \
        """Takes a value which be a single object, or an an iterable
        and returns the content wrapped in a {0}. In the case of strings,
        and dictionaries the original string object is returned in a {0},
        and not as a {0} of chars. A TypeError is raised if this is not
        possible (e.g. dict in frozenset).""".format(cls.__name__)
    _do_safe_coerce.__name__ = \
        "safe_coerce_to_{0}".format(cls.__name__)

    return _do_safe_coerce

safe_coerce_to_tuple = _safe_coerce(tuple)
safe_coerce_to_frozenset = _safe_coerce(frozenset)


def try_cast(value, cast_to):
    try:
        return cast_to(value)
    except (ValueError, TypeError):
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


def get_in(dictionary, keys, default=None):
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
    pairs = itertools.izip(lst, itertools.islice(lst, 1, None))

    return all(x < y for (x, y) in pairs)


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


def fragment(size, lstlike):
    """Faster alternative to grouper for lists/strings."""
    return (lstlike[i:i + size] for i in range(0, len(lstlike), size))


def cumsum(lst, initial=0):
    """Yields the cummulative sums of the values in a
    iterable, starting with the specified initial value."""
    for item in lst:
        initial += item
        yield initial


def fast_pickle_test(obj):
    """Attempts to pickle an object, raising a PicklingError
    if the object is unpicklable. This function uses cPickle
    to determine if the object is pickable, but 'pickle' to
    generate the exception, since the python module produces
    more informative error messages."""
    try:
        cPickle.dumps(obj)
    except (TypeError, cPickle.PicklingError):
        pickle.dumps(obj)
        assert False  # pragma: no coverage


def fill_dict(destination, source):
    """Returns a copy of 'destination' after setting missing key-
    pairs with copies of those of 'source' recursively."""
    if not isinstance(destination, dict) or not isinstance(source, dict):
        raise TypeError("Non-dictionary parameters in 'fill_dict'")

    def _fill_dict(cur_dest, cur_src):
        for key in cur_src:
            if isinstance(cur_src[key], dict) \
                    and isinstance(cur_dest.get(key), dict):
                _fill_dict(cur_dest[key], cur_src[key])
            elif key not in cur_dest:
                cur_dest[key] = cur_src[key]
        return cur_dest

    return _fill_dict(copy.deepcopy(destination), copy.deepcopy(source))


def chain_sorted(*sequences, **kwargs):
    """Chains together sorted sequences, and yields the contents
    in the same order, such that the result is also a sorted sequence.
    The function accepts a 'key'-function keyword, following sort().

    chain_sorted is intended for a few long sequences, and not many short
    sequences. Behavior is undefined if the sequences are not sorted.

    Example:
      >>> tuple(chain_sorted((1, 3, 5), (0, 2, 4)))
      (0, 1, 2, 3, 4, 5)
    """
    key = kwargs.pop('key', None)
    if kwargs:
        raise TypeError("chain_sorted expected keyword 'key', got %r"
                        % (', '.join(kwargs)))

    iterators = []
    for index, sequence_iter in enumerate(map(iter, sequences)):
        try:
            current = sequence_iter.next()
            key_value = current if key is None else key(current)

            iterators.append((key_value, index, current, sequence_iter))
        except StopIteration:
            pass

    heapq.heapify(iterators)

    _len, _heappop, _heapreplace = len, heapq.heappop, heapq.heapreplace

    while _len(iterators) > 1:
        last_key_value, index, current, sequence_iter = iterators[0]
        yield current

        for current in sequence_iter:
            key_value = current if key is None else key(current)

            # Optimization for runs of repeated values
            if key_value != last_key_value:
                _heapreplace(iterators,
                             (key_value, index, current, sequence_iter))
                break
            else:
                yield current
        else:
            # No items remaining in top iterator
            _heappop(iterators)

    if _len(iterators) == 1:
        _, _, current, sequence_iter = iterators[0]

        yield current
        for current in sequence_iter:
            yield current


class Immutable(object):
    """Mixin implementing a immutable class; member variables are specified in
    the init function, cannot be changed afterwards; note that this does not
    prevent changes to the member variables themselves (if not immutable)."""

    def __init__(self, **kwargs):
        object.__init__(self)
        for (key, value) in kwargs.iteritems():
            object.__setattr__(self, key, value)

    def __setattr__(self, _name, _value):
        raise NotImplementedError("Object is immutable")

    def __delattr__(self, _name):
        raise NotImplementedError("Object is immutable")


class TotallyOrdered(object):
    """Mixin implementing a rich-comparison interface, provided
    that the subclass implements the less-than operator (__lt__).
    The __lt__ function should return NotImplemented if the other
    object is not the same type.

    The implementation assumes total order:
    http://en.wikipedia.org/wiki/Total_order
    """

    def __lt__(self, other):
        raise NotImplementedError("__lt__ must be implemented!")

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return not ((self < other) or (other < self))

    def __ne__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (self == other)

    def __le__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (other < self)

    def __ge__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (self < other)

    def __gt__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (other < self)

    # Shut up warning; if hashable, then the subclass will have
    # to implement the __hash__ member function.
    __hash__ = None
