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

import operator
from typing import Callable, Iterable, Sequence, TypeVar

import pytest

import paleomix.common.utilities as utils

T = TypeVar("T")

################################################################################
################################################################################
# Tests for 'safe_coerce_to_tuple'


def test_safe_coerce_to_tuple__str() -> None:
    assert utils.safe_coerce_to_tuple("foo") == ("foo",)


def test_safe_coerce_to_tuple__int() -> None:
    assert utils.safe_coerce_to_tuple(17) == (17,)


def test_safe_coerce_to_tuple__list() -> None:
    assert utils.safe_coerce_to_tuple([1, 3, 2]) == (1, 3, 2)


def test_safe_coerce_to_tuple__tuple() -> None:
    assert utils.safe_coerce_to_tuple((1, 3, 2)) == (1, 3, 2)


def test_safe_coerce_to_tuple__iterable() -> None:
    assert utils.safe_coerce_to_tuple(iter(range(3))) == (0, 1, 2)


def test_safe_coerce_to_tuple__dict() -> None:
    assert utils.safe_coerce_to_tuple({1: 2, 3: 4}) == (1, 3)


###############################################################################
###############################################################################
# Tests for 'safe_coerce_to_frozenset'


def test_safe_coerce_to_frozenset__str() -> None:
    assert utils.safe_coerce_to_frozenset("foo") == frozenset(("foo",))


def test_safe_coerce_to_frozenset__unicode() -> None:
    assert utils.safe_coerce_to_frozenset("foo") == frozenset(("foo",))


def test_safe_coerce_to_frozenset__int() -> None:
    assert utils.safe_coerce_to_frozenset(17) == frozenset((17,))


def test_safe_coerce_to_frozenset__list() -> None:
    assert utils.safe_coerce_to_frozenset([1, 3, 2]) == frozenset((1, 3, 2))


def test_safe_coerce_to_frozenset__tuple() -> None:
    assert utils.safe_coerce_to_frozenset((1, 3, 2)) == frozenset((1, 3, 2))


def test_safe_coerce_to_frozenset__iterable() -> None:
    assert utils.safe_coerce_to_frozenset(range(3)) == frozenset((0, 1, 2))


def test_safe_coerce_to_frozenset__dict() -> None:
    assert utils.safe_coerce_to_frozenset({1: 2, 3: 4}) == frozenset((1, 3))


###############################################################################
###############################################################################
# Tests for 'try_cast'


def test_try_cast__int_to_int() -> None:
    assert utils.try_cast(17, int) == 17


def test_try_cast__float_to_int() -> None:
    assert utils.try_cast(17.3, int) == 17


def test_try_cast__good_str_to_int() -> None:
    assert utils.try_cast("17", int) == 17


def test_try_cast__bad_str_to_int() -> None:
    assert utils.try_cast("x17", int) == "x17"


def test_try_cast__list_to_int() -> None:
    assert utils.try_cast([1, 2, 3], int) == [1, 2, 3]


###############################################################################
###############################################################################
# Tests for 'set_in'


def test_set_in__single_kw_in_empty_dictionary() -> None:
    value: dict[object, object] = {}
    utils.set_in(value, ["Foo"], 17)
    assert value == {"Foo": 17}


def test_set_in__two_kws_in_empty_dictionary() -> None:
    value: dict[object, object] = {}
    utils.set_in(value, ["Foo", 13], 17)
    assert value == {"Foo": {13: 17}}


def test_set_in__three_kws_in_empty_dictionary() -> None:
    value: dict[object, object] = {}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert value == {"Foo": {13: {(1, 2): 17}}}


def test_set_in__three_kws_in_partial_dictionary() -> None:
    value: dict[object, object] = {"Foo": {12: 0}}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert value == {"Foo": {12: 0, 13: {(1, 2): 17}}}

    value: dict[object, object] = {"Foo": {13: {"Bar": None}}}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert value == {"Foo": {13: {(1, 2): 17, "Bar": None}}}


def test_set_in__update_value_one_kw() -> None:
    value: dict[object, object] = {1: None}
    utils.set_in(value, [1], 3.14)
    assert value == {1: 3.14}


def test_set_in__update_value_two_kw() -> None:
    value: dict[object, object] = {1: {2: 3}}
    utils.set_in(value, [1, 2], 365)
    assert value == {1: {2: 365}}


def test_set_in__fail_on_no_kws() -> None:
    with pytest.raises(ValueError, match="No keys passed to 'set_in'"):
        utils.set_in({}, [], 17)


def test_set_in__fail_on_invalid_sub_dictionary_first_level() -> None:
    with pytest.raises(TypeError):
        utils.set_in(None, [1], 17)  # pyright: ignore[reportGeneralTypeIssues]


def test_set_in__fail_on_invalid_sub_dictionary_second_level() -> None:
    with pytest.raises(TypeError):
        utils.set_in({1: None}, [1, 2], 17)


def test_set_in__fail_on_invalid_sub_dictionary_third_level() -> None:
    with pytest.raises(TypeError):
        utils.set_in({1: {2: None}}, [1, 2, 3], 17)


def test_set_in__iteratable_keywords() -> None:
    value: dict[object, object] = {}
    utils.set_in(value, iter(["Foo", 13, (1, 2)]), 17)
    assert value == {"Foo": {13: {(1, 2): 17}}}


###############################################################################
###############################################################################
# Tests for 'get_in'


def test_get_in__get_value_one_keyword() -> None:
    assert utils.get_in({1: 2}, [1]) == 2


def test_get_in__get_value_two_keywords() -> None:
    assert utils.get_in({1: {2: 3}}, [1, 2]) == 3


def test_get_in__get_value_three_keywords() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, [1, 2, 3]) == 4


def test_get_in__get_default_one_keyword() -> None:
    assert utils.get_in({1: 2}, [2]) is None


def test_get_in__get_default_one_keyword_with_default() -> None:
    assert utils.get_in({1: 2}, [2], "other") == "other"


def test_get_in__get_default_three_keywords_fail_at_first() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, [2, 2, 4]) is None


def test_get_in__get_default_three_keywords_fail_at_first_with_default() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, [2, 2, 4], "other") == "other"


def test_get_in__get_default_three_keywords_fail_at_second() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, [1, 3, 4]) is None


def test_get_in__get_default_three_keywords_fail_at_second_with_default() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, [1, 3, 4], "other") == "other"


def test_get_in__get_default_three_keywords_fail_at_third() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, [1, 2, 4]) is None


def test_get_in__get_default_three_keywords_fail_at_third_with_default() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, [1, 2, 4], "other") == "other"


def test_get_in__iterator_keywords() -> None:
    assert utils.get_in({1: {2: {3: 4}}}, iter([1, 2, 3])) == 4


###############################################################################
###############################################################################
# Tests for 'split_before'


def split_before(iterable: Iterable[T], pred: Callable[[T], bool]) -> list[list[T]]:
    # Convertion to list allows the implementation to be
    # lazy, while making comparisons for asserts easier
    return list(utils.split_before(iterable, pred))


def test_split_before__split_empty_list() -> None:
    empty_list: list[int] = []
    assert split_before(empty_list, bool) == []


def test_split_before__split_list_with_no_true_pred() -> None:
    assert split_before(list(range(10)), lambda x: False) == [list(range(10))]


def test_split_before__split_list_true_pred_at_first_position() -> None:
    assert split_before(list(range(4)), lambda x: x % 2 == 0) == [[0, 1], [2, 3]]


def test_split_before__split_list_true_pred_at_second_position() -> None:
    assert split_before(list(range(4)), lambda x: x % 2 == 1) == [[0], [1, 2], [3]]


def test_split_before__split_consequtive_true_pred() -> None:
    assert split_before(list(range(0, 5, 2)), lambda x: x % 2 == 0) == [[0], [2], [4]]


def test_split_before__no_hits() -> None:
    assert split_before(list(range(1, 5)), lambda x: x % 5 == 0) == [list(range(1, 5))]


###############################################################################
###############################################################################
# Tests for 'grouper'


def test_grouper__empty_list() -> None:
    empty_list: list[int] = []
    result = utils.grouper(3, empty_list)
    assert list(result) == []


def test_grouper__non_empty_list() -> None:
    result = utils.grouper(3, list(range(6)))
    expected = [(0, 1, 2), (3, 4, 5)]
    assert list(result) == expected


def test_grouper__non_empty_list_with_trailing() -> None:
    result = utils.grouper(3, list(range(7)))
    expected = [(0, 1, 2), (3, 4, 5), (6, None, None)]
    assert list(result) == expected


def test_grouper__non_empty_list_with_trailing_fill_value() -> None:
    result = utils.grouper(3, list(range(7)), fillvalue=r"\0")
    expected = [(0, 1, 2), (3, 4, 5), (6, r"\0", r"\0")]
    assert list(result) == expected


###############################################################################
###############################################################################
# Tests for 'group_by_pred'


def test_group_by_pred__empty_list() -> None:
    value: list[int] = []
    assert utils.group_by_pred(bool, value) == ([], [])


def test_group_by_pred__always_false() -> None:
    assert utils.group_by_pred(lambda x: False, [1, 2, 3]) == ([], [1, 2, 3])


def test_group_by_pred__always_true() -> None:
    assert utils.group_by_pred(lambda x: True, [1, 2, 3]) == ([1, 2, 3], [])


def test_group_by_pred__is_even() -> None:
    assert utils.group_by_pred(lambda x: x % 2 == 0, [1, 2, 3]) == ([2], [1, 3])


def test_group_by_pred__iterable() -> None:
    assert utils.group_by_pred(lambda x: x % 2 == 0, range(1, 4)) == ([2], [1, 3])


###############################################################################
###############################################################################
# Tests for 'fragment'


def fragment(
    size: int,
    items: Sequence[T],
) -> list[Sequence[T]]:
    """Faster alternative to grouper for lists/strings."""
    return list(utils.fragment(size, items))


def test_fragment__empty() -> None:
    empty_list: list[int] = []
    assert fragment(5, "") == []
    assert fragment(5, empty_list) == []


def test_fragment__partial_fragment() -> None:
    assert fragment(3, "ab") == ["ab"]
    assert fragment(3, ["a", "b"]) == [["a", "b"]]


def test_fragment__single_fragment() -> None:
    assert fragment(3, "abc") == ["abc"]
    assert fragment(3, ["a", "b", "c"]) == [["a", "b", "c"]]


def test_fragment__multiple_fragments() -> None:
    assert fragment(3, "abcdef") == ["abc", "def"]
    assert fragment(3, list("abcdef")) == [list("abc"), list("def")]


def test_fragment__multiple_fragments_partial() -> None:
    assert fragment(3, "abcdefgh") == ["abc", "def", "gh"]
    assert fragment(3, list("abcdefgh")) == [
        list("abc"),
        list("def"),
        list("gh"),
    ]


def test_fragment__range() -> None:
    assert list(utils.fragment(3, range(6))) == [range(3), range(3, 6)]


def test_fragment__iterable() -> None:
    with pytest.raises(TypeError):
        fragment(3, iter(range(6)))  # pyright: ignore[reportGeneralTypeIssues]


def test_fragment__set() -> None:
    with pytest.raises(TypeError):
        fragment(3, set(range(6)))  # pyright: ignore[reportGeneralTypeIssues]


###############################################################################
###############################################################################
# fill_dict


def test_fill_dict__empty_dicts() -> None:
    result = utils.fill_dict({}, {})
    assert result == {}


def test_fill_dict__filling_empty_dict() -> None:
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    expected = {"a": 1, "b": {"c": 2, "d": 3}}
    result = utils.fill_dict({}, source)
    assert result == expected


def test_fill_dict__filling_full_dict() -> None:
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    destination = {"a": 2, "b": {"c": 3, "d": 4}}
    expected = {"a": 2, "b": {"c": 3, "d": 4}}
    result = utils.fill_dict(destination, source)
    assert result == expected


def test_fill_dict__destination_not_modified() -> None:
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    destination = {"b": {"d": 0}}
    utils.fill_dict(destination, source)
    assert destination == {"b": {"d": 0}}


def test_fill_dict__source_not_modified() -> None:
    expected = {"a": 1, "b": {"c": 2, "d": 3}}
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    destination = {"b": {"d": 0}}
    utils.fill_dict(destination, source)
    assert source == expected


def test_fill_dict__destination_must_be_dict() -> None:
    with pytest.raises(TypeError):
        utils.fill_dict([], {})  # pyright: ignore[reportGeneralTypeIssues]


def test_fill_dict__source_must_be_dict() -> None:
    with pytest.raises(TypeError):
        utils.fill_dict({}, [])  # pyright: ignore[reportGeneralTypeIssues]


###############################################################################
###############################################################################
# Immutable


def test_immutable__properties_set() -> None:
    class ImmutableCls(utils.Immutable):
        value: int

        def __init__(self, value: int) -> None:
            utils.Immutable.__init__(self, value=value)

    obj = ImmutableCls(17)
    assert obj.value == 17


@pytest.mark.parametrize(("key", "value"), [("value", 13), ("new_value", "foo")])
def test_immutable__properties_immutable(key: str, value: object) -> None:
    class ImmutableCls(utils.Immutable):
        def __init__(self, value: int) -> None:
            utils.Immutable.__init__(self, value=value)

    obj = ImmutableCls(17)
    with pytest.raises(NotImplementedError):
        setattr(obj, key, value)


def test_immutable__properties_del() -> None:
    class ImmutableCls(utils.Immutable):
        def __init__(self, value: int) -> None:
            utils.Immutable.__init__(self, value=value)

    obj = ImmutableCls(17)
    with pytest.raises(NotImplementedError):
        del obj.value


###############################################################################
###############################################################################
# TotallyOrdered


class SomethingOrdered(utils.TotallyOrdered):
    def __init__(self, value: int) -> None:
        self.value = value

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, SomethingOrdered):
            return NotImplemented

        return self.value < other.value


# Less than
def test_totally_ordered__lt_vs_lt() -> None:
    assert SomethingOrdered(1) < SomethingOrdered(2)


def test_totally_ordered__lt_vs_gt() -> None:
    assert not (SomethingOrdered(1) < SomethingOrdered(0))


def test_totally_ordered__lt_vs_wrong_type() -> None:
    assert SomethingOrdered(1).__lt__("Foo") == NotImplemented


# Less than or equal
def test_totally_ordered__le_vs_le() -> None:
    assert SomethingOrdered(1) <= SomethingOrdered(1)


def test_totally_ordered__le_vs_gt() -> None:
    assert not (SomethingOrdered(1) <= SomethingOrdered(0))


def test_totally_ordered__le_vs_wrong_type() -> None:
    assert SomethingOrdered(1).__le__("Foo") == NotImplemented


# Greater than or equal
def test_totally_ordered__ge_vs_ge() -> None:
    assert SomethingOrdered(1) >= SomethingOrdered(1)


def test_totally_ordered__ge_vs_lt() -> None:
    assert not (SomethingOrdered(0) >= SomethingOrdered(1))


def test_totally_ordered__ge_vs_wrong_type() -> None:
    assert SomethingOrdered(1).__ge__("Foo") == NotImplemented


# Greater than
def test_totally_ordered__gt_vs_gt() -> None:
    assert SomethingOrdered(1) > SomethingOrdered(0)


def test_totally_ordered__gt_vs_eq() -> None:
    assert not (SomethingOrdered(0) > SomethingOrdered(0))


def test_totally_ordered__gt_vs_wrong_type() -> None:
    assert SomethingOrdered(1).__gt__("Foo") == NotImplemented


# Equal to
def test_totally_ordered__eq_vs_eq() -> None:
    assert SomethingOrdered(1) == SomethingOrdered(1)


def test_totally_ordered__eq_vs_ne() -> None:
    assert not (SomethingOrdered(1) == SomethingOrdered(2))  # noqa: SIM201


def test_totally_ordered__eq_vs_wrong_type() -> None:
    assert SomethingOrdered(1).__eq__("Foo") == NotImplemented


# Not equal to
def test_totally_ordered__ne_vs_ne() -> None:
    assert SomethingOrdered(1) != SomethingOrdered(2)


def test_totally_ordered__ne_vs_eq() -> None:
    assert not (SomethingOrdered(1) != SomethingOrdered(1))  # noqa: SIM202


def test_totally_ordered__ne_vs_wrong_type() -> None:
    assert SomethingOrdered(1).__ne__("Foo") == NotImplemented


class SomethingBadlyOrdered(utils.TotallyOrdered):
    def __init__(self, value: int) -> None:
        self.value = value


def test_totally_ordered__missing_implementation() -> None:
    obj_a = SomethingBadlyOrdered(1)
    obj_b = SomethingBadlyOrdered(2)
    with pytest.raises(NotImplementedError):
        operator.gt(obj_a, obj_b)
