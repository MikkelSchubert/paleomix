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
    Iterable,
    Iterator,
    Sequence,
    TypeVar,
    overload,
)

T = TypeVar("T")


def safe_coerce_to_tuple(value: object) -> tuple[Any, ...]:
    """Convert value to a tuple, unless it is a string or a non-sequence, in which case
    it is return as a single-element tuple."""
    if isinstance(value, (str, bytes)):
        return (value,)

    try:
        return tuple(value)  # pyright: ignore[reportArgumentType]
    except TypeError:
        return (value,)


def safe_coerce_to_frozenset(value: object) -> frozenset[Any]:
    """Convert value to a tuple, unless it is a string or a non-sequence, in which case
    it is return as a single-element tuple."""
    if isinstance(value, (str, bytes)):
        return frozenset((value,))

    try:
        return frozenset(value)  # pyright: ignore[reportArgumentType]
    except TypeError:
        return frozenset((value,))


def try_cast(value: object, cast_to: type) -> object:
    try:
        return cast_to(value)
    except (ValueError, TypeError):
        return value


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
def fragment(size: int, items: str) -> Iterable[str]: ...


@overload
def fragment(size: int, items: bytes) -> Iterable[bytes]: ...


@overload
def fragment(size: int, items: Sequence[T]) -> Iterable[Sequence[T]]: ...


def fragment(size: int, items: Sequence[T]) -> Iterable[Sequence[T]]:
    """Faster alternative to grouper for lists/strings."""
    return (items[i : i + size] for i in range(0, len(items), size))


def fill_dict(destination: dict[Any, Any], source: dict[Any, Any]) -> dict[Any, Any]:
    """Returns a copy of 'destination' after setting missing key-
    pairs with copies of those of 'source' recursively."""
    if not isinstance(destination, dict) or not isinstance(source, dict):
        raise TypeError("Non-dictionary parameters in 'fill_dict'")

    def _fill_dict(cur_dest: dict[Any, Any], cur_src: dict[Any, Any]) -> dict[Any, Any]:
        for key in cur_src:
            if isinstance(cur_src[key], dict) and isinstance(cur_dest.get(key), dict):
                _fill_dict(cur_dest[key], cur_src[key])
            elif key not in cur_dest:
                cur_dest[key] = cur_src[key]
        return cur_dest

    return _fill_dict(copy.deepcopy(destination), copy.deepcopy(source))


# FIXME: Properly handle slots, type-safety
class Immutable:
    """Mixin implementing a immutable class; member variables are specified in
    the init function, cannot be changed afterwards; note that this does not
    prevent changes to the member variables themselves (if not immutable)."""

    def __init__(self, **kwargs: object) -> None:
        if hasattr(self, "__slots__"):
            raise AssertionError("Immutable does not support slots")

        object.__init__(self)
        for key, value in kwargs.items():
            object.__setattr__(self, key, value)

    def __setattr__(self, _name: str, _value: object) -> None:
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

    def __lt__(self, other: object) -> bool:
        raise NotImplementedError("__lt__ must be implemented!")

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not ((self < other) or (other < self))

    def __ne__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (self == other)

    def __le__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (other < self)

    def __ge__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return not (self < other)

    def __gt__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return other < self
