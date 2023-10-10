#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import copy
import itertools
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
    Hashable,
    Iterable,
    Iterator,
    Sequence,
    Tuple,
    TypeVar,
    overload,
)

T = TypeVar("T")


def safe_coerce_to_tuple(value: Any) -> Tuple[Any, ...]:
    """Convert value to a tuple, unless it is a string or a non-sequence, in which case
    it is return as a single-element tuple."""
    if isinstance(value, str):
        return (value,)

    try:
        return tuple(value)
    except TypeError:
        return (value,)


def safe_coerce_to_frozenset(value: Any) -> FrozenSet[Any]:
    """Convert value to a tuple, unless it is a string or a non-sequence, in which case
    it is return as a single-element tuple."""
    if isinstance(value, str):
        return frozenset((value,))

    try:
        return frozenset(value)
    except TypeError:
        return frozenset((value,))


def try_cast(value: Any, cast_to: type) -> Any:
    try:
        return cast_to(value)
    except (ValueError, TypeError):
        return value


def set_in(
    dictionary: Dict[object, object],
    keys: Iterable[Hashable],
    value: object,
) -> None:
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
            child = dictionary[key]
            if not isinstance(child, dict):
                raise TypeError(child)

            dictionary = child
        except KeyError:
            new_dict = {}  # type: Dict[object, object]
            dictionary[key] = new_dict
            dictionary = new_dict

    dictionary[keys[-1]] = value


def get_in(
    dictionary: Dict[Any, Any],
    keys: Iterable[Hashable],
    default: Any = None,
) -> Any:
    """Traverses a set of nested dictionaries using the keys in
    kws, and returns the value assigned to the final keyword
    in the innermost dictionary. Calling get_in(d, [X, Y])
    is equivalent to calling d.get(X).get(Y), with the
    difference that any missing keys causes the default value
    to be returned.

    Behavior on non-dictionaries is undefined."""
    keys = list(keys)
    for key in keys[:-1]:
        try:
            dictionary = dictionary[key]
        except KeyError:
            return default

    return dictionary.get(keys[-1], default)


def split_before(iterable: Iterable[T], pred: Callable[[T], bool]) -> Iterator[list[T]]:
    """Takes a sequence and splits it before every value where pred(v) is true.
    Thus split_before(range(10), key = lambda x: x % 2 == 0) would return the
    sequence [[1], [2,3], [4,5], [6,7], [7,8], [9]]"""
    items: list[T] = []
    for value in iterable:
        if pred(value) and items:
            yield items
            items = []
        items.append(value)

    if items:
        yield items


# Copied from the Python 'itertools' module documentation
def grouper(
    size: int,
    iterable: Iterable[T],
    fillvalue: object = None,
) -> itertools.zip_longest[tuple[T]]:
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * size
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def group_by_pred(
    pred: Callable[[T], bool],
    iterable: Iterable[T],
) -> tuple[list[T], list[T]]:
    """Splits items in a sequence into two lists, one containing
    items matching the predicate, and another containing those that
    do not."""
    is_true: list[T] = []
    is_false: list[T] = []
    for item in iterable:
        if pred(item):
            is_true.append(item)
        else:
            is_false.append(item)

    return is_true, is_false


@overload
def fragment(size: int, items: str) -> Iterable[str]:
    ...


@overload
def fragment(size: int, items: bytes) -> Iterable[bytes]:
    ...


@overload
def fragment(size: int, items: Sequence[T]) -> Iterable[Sequence[T]]:
    ...


def fragment(size: int, items: Sequence[T]) -> Iterable[Sequence[T]]:
    """Faster alternative to grouper for lists/strings."""
    return (items[i : i + size] for i in range(0, len(items), size))


def fill_dict(destination: Dict[Any, Any], source: Dict[Any, Any]) -> Dict[Any, Any]:
    """Returns a copy of 'destination' after setting missing key-
    pairs with copies of those of 'source' recursively."""
    if not isinstance(destination, dict) or not isinstance(source, dict):
        raise TypeError("Non-dictionary parameters in 'fill_dict'")

    def _fill_dict(cur_dest: Dict[Any, Any], cur_src: Dict[Any, Any]) -> Dict[Any, Any]:
        for key in cur_src:
            if isinstance(cur_src[key], dict) and isinstance(cur_dest.get(key), dict):
                _fill_dict(cur_dest[key], cur_src[key])
            elif key not in cur_dest:
                cur_dest[key] = cur_src[key]
        return cur_dest

    return _fill_dict(copy.deepcopy(destination), copy.deepcopy(source))


class Immutable:
    """Mixin implementing a immutable class; member variables are specified in
    the init function, cannot be changed afterwards; note that this does not
    prevent changes to the member variables themselves (if not immutable)."""

    def __init__(self, **kwargs: Any):
        object.__init__(self)
        for key, value in kwargs.items():
            object.__setattr__(self, key, value)

    def __setattr__(self, _name: str, _value: Any) -> None:
        raise NotImplementedError("Object is immutable")

    def __delattr__(self, _name: str) -> None:
        raise NotImplementedError("Object is immutable")


class TotallyOrdered:
    """Mixin implementing a rich-comparison interface, provided
    that the subclass implements the less-than operator (__lt__).
    The __lt__ function should return NotImplemented if the other
    object is not the same type.

    The implementation assumes total order:
    http://en.wikipedia.org/wiki/Total_order
    """

    def __lt__(self, other: Any) -> bool:
        raise NotImplementedError("__lt__ must be implemented!")

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not ((self < other) or (other < self))

    def __ne__(self, other: Any) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (self == other)

    def __le__(self, other: Any) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (other < self)

    def __ge__(self, other: Any) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (self < other)

    def __gt__(self, other: Any) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return other < self
