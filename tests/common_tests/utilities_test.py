#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import operator

import pytest

import paleomix.common.utilities as utils

################################################################################
################################################################################
# Tests for 'safe_coerce_to_tuple'


def test_safe_coerce_to_tuple__str():
    assert utils.safe_coerce_to_tuple("foo") == ("foo",)


def test_safe_coerce_to_tuple__int():
    assert utils.safe_coerce_to_tuple(17) == (17,)


def test_safe_coerce_to_tuple__list():
    assert utils.safe_coerce_to_tuple([1, 3, 2]) == (1, 3, 2)


def test_safe_coerce_to_tuple__tuple():
    assert utils.safe_coerce_to_tuple((1, 3, 2)) == (1, 3, 2)


def test_safe_coerce_to_tuple__iterable():
    assert utils.safe_coerce_to_tuple(iter(range(3))) == (0, 1, 2)


def test_safe_coerce_to_tuple__dict():
    assert utils.safe_coerce_to_tuple({1: 2, 3: 4}) == ({1: 2, 3: 4},)


###############################################################################
###############################################################################
# Tests for 'safe_coerce_to_frozenset'


def test_safe_coerce_to_frozenset__str():
    assert utils.safe_coerce_to_frozenset("foo") == frozenset(("foo",))


def test_safe_coerce_to_frozenset__unicode():
    assert utils.safe_coerce_to_frozenset("foo") == frozenset(("foo",))


def test_safe_coerce_to_frozenset__int():
    assert utils.safe_coerce_to_frozenset(17) == frozenset((17,))


def test_safe_coerce_to_frozenset__list():
    assert utils.safe_coerce_to_frozenset([1, 3, 2]) == frozenset((1, 3, 2))


def test_safe_coerce_to_frozenset__tuple():
    assert utils.safe_coerce_to_frozenset((1, 3, 2)) == frozenset(((1, 3, 2)))


def test_safe_coerce_to_frozenset__iterable():
    assert utils.safe_coerce_to_frozenset(range(3)) == frozenset((0, 1, 2))


def test_safe_coerce_to_frozenset__dict():
    with pytest.raises(TypeError):
        utils.safe_coerce_to_frozenset({1: 2, 3: 4})


###############################################################################
###############################################################################
# Tests for 'try_cast'


def test_try_cast__int_to_int():
    assert utils.try_cast(17, int) == 17


def test_try_cast__float_to_int():
    assert utils.try_cast(17.3, int) == 17


def test_try_cast__good_str_to_int():
    assert utils.try_cast("17", int) == 17


def test_try_cast__bad_str_to_int():
    assert utils.try_cast("x17", int) == "x17"


def test_try_cast__list_to_int():
    assert utils.try_cast([1, 2, 3], int) == [1, 2, 3]


###############################################################################
###############################################################################
# Tests for 'set_in'


def test_set_in__single_kw_in_empty_dictionary():
    value = {}
    utils.set_in(value, ["Foo"], 17)
    assert value == {"Foo": 17}


def test_set_in__two_kws_in_empty_dictionary():
    value = {}
    utils.set_in(value, ["Foo", 13], 17)
    assert value == {"Foo": {13: 17}}


def test_set_in__three_kws_in_empty_dictionary():
    value = {}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert value == {"Foo": {13: {(1, 2): 17}}}


def test_set_in__three_kws_in_partial_dictionary():
    value = {"Foo": {12: 0}}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert value == {"Foo": {12: 0, 13: {(1, 2): 17}}}

    value = {"Foo": {13: {"Bar": None}}}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert value == {"Foo": {13: {(1, 2): 17, "Bar": None}}}


def test_set_in__update_value_one_kw():
    value = {1: None}
    utils.set_in(value, [1], 3.14)
    assert value == {1: 3.14}


def test_set_in__update_value_two_kw():
    value = {1: {2: 3}}
    utils.set_in(value, [1, 2], 365)
    assert value == {1: {2: 365}}


def test_set_in__fail_on_no_kws():
    with pytest.raises(ValueError):
        utils.set_in({}, [], 17)


def test_set_in__fail_on_invalid_sub_dictionary_first_level():
    with pytest.raises(TypeError):
        utils.set_in(None, [1], 17)


def test_set_in__fail_on_invalid_sub_dictionary_second_level():
    with pytest.raises(TypeError):
        utils.set_in({1: None}, [1, 2], 17)


def test_set_in__fail_on_invalid_sub_dictionary_third_level():
    with pytest.raises(TypeError):
        utils.set_in({1: {2: None}}, [1, 2, 3], 17)


def test_set_in__iteratable_keywords():
    value = {}
    utils.set_in(value, iter(["Foo", 13, (1, 2)]), 17)
    assert value == {"Foo": {13: {(1, 2): 17}}}


###############################################################################
###############################################################################
# Tests for 'get_in'


def test_get_in__get_value_one_keyword():
    assert utils.get_in({1: 2}, [1]) == 2


def test_get_in__get_value_two_keywords():
    assert utils.get_in({1: {2: 3}}, [1, 2]) == 3


def test_get_in__get_value_three_keywords():
    assert utils.get_in({1: {2: {3: 4}}}, [1, 2, 3]) == 4


def test_get_in__get_default_one_keyword():
    assert utils.get_in({1: 2}, [2]) is None


def test_get_in__get_default_one_keyword_with_default():
    assert utils.get_in({1: 2}, [2], "other") == "other"


def test_get_in__get_default_three_keywords_fail_at_first():
    assert utils.get_in({1: {2: {3: 4}}}, [2, 2, 4]) is None


def test_get_in__get_default_three_keywords_fail_at_first_with_default():
    assert utils.get_in({1: {2: {3: 4}}}, [2, 2, 4], "other") == "other"


def test_get_in__get_default_three_keywords_fail_at_second():
    assert utils.get_in({1: {2: {3: 4}}}, [1, 3, 4]) is None


def test_get_in__get_default_three_keywords_fail_at_second_with_default():
    assert utils.get_in({1: {2: {3: 4}}}, [1, 3, 4], "other") == "other"


def test_get_in__get_default_three_keywords_fail_at_third():
    assert utils.get_in({1: {2: {3: 4}}}, [1, 2, 4]) is None


def test_get_in__get_default_three_keywords_fail_at_third_with_default():
    assert utils.get_in({1: {2: {3: 4}}}, [1, 2, 4], "other") == "other"


def test_get_in__iterator_keywords():
    assert utils.get_in({1: {2: {3: 4}}}, iter([1, 2, 3])) == 4


###############################################################################
###############################################################################
# Tests for 'split_before'


def _do_split(lst, key):
    # Conversion to list allows the implementation to be
    # lazy, while making comparisons for asserts easier
    return list(utils.split_before(lst, key))


def test_split_before__split_empty_list():
    assert _do_split([], None) == []


def test_split_before__split_list_with_no_true_pred():
    assert _do_split(list(range(10)), lambda x: False) == [list(range(10))]


def test_split_before__split_list_true_pred_at_first_position():
    assert _do_split(list(range(4)), lambda x: x % 2 == 0) == [[0, 1], [2, 3]]


def test_split_before__split_list_true_pred_at_second_position():
    assert _do_split(list(range(4)), lambda x: x % 2 == 1) == [[0], [1, 2], [3]]


def test_split_before__split_consequtive_true_pred():
    assert _do_split(list(range(0, 5, 2)), lambda x: x % 2 == 0) == [[0], [2], [4]]


def test_split_before__no_hits():
    assert _do_split(list(range(1, 5)), lambda x: x % 5 == 0) == [list(range(1, 5))]


###############################################################################
###############################################################################
# Tests for 'grouper'


def test_grouper__empty_list():
    result = utils.grouper(3, [])
    assert list(result) == []


def test_grouper__non_empty_list():
    result = utils.grouper(3, list(range(6)))
    expected = [(0, 1, 2), (3, 4, 5)]
    assert list(result) == expected


def test_grouper__non_empty_list_with_trailing():
    result = utils.grouper(3, list(range(7)))
    expected = [(0, 1, 2), (3, 4, 5), (6, None, None)]
    assert list(result) == expected


def test_grouper__non_empty_list_with_trailing_fill_value():
    result = utils.grouper(3, list(range(7)), fillvalue=r"\0")
    expected = [(0, 1, 2), (3, 4, 5), (6, r"\0", r"\0")]
    assert list(result) == expected


###############################################################################
###############################################################################
# Tests for 'group_by_pred'


def test_group_by_pred__empty_list():
    assert utils.group_by_pred(id, []) == ([], [])


def test_group_by_pred__always_false():
    assert utils.group_by_pred(lambda x: False, [1, 2, 3]) == ([], [1, 2, 3])


def test_group_by_pred__always_true():
    assert utils.group_by_pred(lambda x: True, [1, 2, 3]) == ([1, 2, 3], [])


def test_group_by_pred__is_even():
    assert utils.group_by_pred(lambda x: x % 2 == 0, [1, 2, 3]) == ([2], [1, 3])


def test_group_by_pred__iterable():
    assert utils.group_by_pred(lambda x: x % 2 == 0, range(1, 4)) == ([2], [1, 3])


###############################################################################
###############################################################################
# Tests for 'fragment'


def test_fragment__empty():
    assert list(utils.fragment(5, "")) == []
    assert list(utils.fragment(5, [])) == []


def test_fragment__partial_fragment():
    assert list(utils.fragment(3, "ab")) == ["ab"]
    assert list(utils.fragment(3, ["a", "b"])) == [["a", "b"]]


def test_fragment__single_fragment():
    assert list(utils.fragment(3, "abc")) == ["abc"]
    assert list(utils.fragment(3, ["a", "b", "c"])) == [["a", "b", "c"]]


def test_fragment__multiple_fragments():
    assert list(utils.fragment(3, "abcdef")) == ["abc", "def"]
    assert list(utils.fragment(3, list("abcdef"))) == [list("abc"), list("def")]


def test_fragment__multiple_fragments_partial():
    assert list(utils.fragment(3, "abcdefgh")) == ["abc", "def", "gh"]
    assert list(utils.fragment(3, list("abcdefgh"))) == [
        list("abc"),
        list("def"),
        list("gh"),
    ]


def test_fragment__range():
    assert list(utils.fragment(3, range(6))) == [range(3), range(3, 6)]


def test_fragment__iterable():
    with pytest.raises(TypeError):
        list(utils.fragment(3, iter(range(6))))


def test_fragment__set():
    with pytest.raises(TypeError):
        list(utils.fragment(3, set(range(6))))


###############################################################################
###############################################################################
# Tests for 'cumsum'


def test_cumsum__empty():
    assert list(utils.cumsum([])) == []


def test_cumsum__integers():
    assert list(utils.cumsum(list(range(-4, 5)))) == [
        -4,
        -7,
        -9,
        -10,
        -10,
        -9,
        -7,
        -4,
        0,
    ]


def test_cumsum__float():
    assert list(utils.cumsum((1.0, 2.0, 3.0))) == [1.0, 3.0, 6.0]


def test_cumsum__initial():
    assert list(utils.cumsum(list(range(5)), -10)) == [-10, -9, -7, -4, 0]


###############################################################################
###############################################################################
# fill_dict


def test_fill_dict__empty_dicts():
    result = utils.fill_dict({}, {})
    assert result == {}


def test_fill_dict__filling_empty_dict():
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    expected = {"a": 1, "b": {"c": 2, "d": 3}}
    result = utils.fill_dict({}, source)
    assert result == expected


def test_fill_dict__filling_full_dict():
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    destination = {"a": 2, "b": {"c": 3, "d": 4}}
    expected = {"a": 2, "b": {"c": 3, "d": 4}}
    result = utils.fill_dict(destination, source)
    assert result == expected


def test_fill_dict__destination_not_modified():
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    destination = {"b": {"d": 0}}
    utils.fill_dict(destination, source)
    assert destination == {"b": {"d": 0}}


def test_fill_dict__source_not_modified():
    expected = {"a": 1, "b": {"c": 2, "d": 3}}
    source = {"a": 1, "b": {"c": 2, "d": 3}}
    destination = {"b": {"d": 0}}
    utils.fill_dict(destination, source)
    assert source == expected


def test_fill_dict__destination_must_be_dict():
    with pytest.raises(TypeError):
        utils.fill_dict([], {})


def test_fill_dict__source_must_be_dict():
    with pytest.raises(TypeError):
        utils.fill_dict({}, [])


###############################################################################
###############################################################################
# chain_sorted


def test_chain_sorted__no_sequences():
    expected = ()
    result = tuple(utils.chain_sorted())
    assert expected == result


def test_chain_sorted__single_sequence():
    sequence = (1, 2, 3)
    result = tuple(utils.chain_sorted(sequence))
    assert sequence == result


_SEQUENTIAL_CONTENT_1 = (1, 2, 3)
_SEQUENTIAL_CONTENT_2 = (4, 5, 6)
_SEQUENTIAL_CONTENT_PERMUTATIONS = (
    (_SEQUENTIAL_CONTENT_1, _SEQUENTIAL_CONTENT_2),
    (_SEQUENTIAL_CONTENT_2, _SEQUENTIAL_CONTENT_1),
)


@pytest.mark.parametrize("seq_a, seq_b", _SEQUENTIAL_CONTENT_PERMUTATIONS)
def test_chain_sorted__sequential_contents(seq_a, seq_b):
    expected = (1, 2, 3, 4, 5, 6)
    result = tuple(utils.chain_sorted(seq_a, seq_b))
    assert expected == result


def test_chain_sorted__mixed_contents():
    sequence_a = (3, 4, 8)
    sequence_c = (0, 1, 6)
    sequence_b = (2, 5, 7)
    expected = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    result = tuple(utils.chain_sorted(sequence_a, sequence_b, sequence_c))
    assert expected == result


def test_chain_sorted__mixed_length_contents():
    sequence_a = (1,)
    sequence_c = (0, 2)
    sequence_b = ()
    expected = (0, 1, 2)
    result = tuple(utils.chain_sorted(sequence_a, sequence_b, sequence_c))
    assert expected == result


def test_chain_sorted__mixed_contents__key():
    sequence_a = (-2, -3, -5)
    sequence_b = (0, -1, -4)
    expected = (0, -1, -2, -3, -4, -5)
    result = tuple(utils.chain_sorted(sequence_a, sequence_b, key=abs))
    assert expected == result


def test_chain_sorted__identical_objects_are_preserved():
    object_a = [1]
    object_b = [1]
    assert object_a is not object_b
    expected = (object_a, object_b)
    result = tuple(utils.chain_sorted([object_a], [object_b]))
    assert expected == result
    assert object_a is result[0] or object_a is result[1]
    assert object_b is result[0] or object_b is result[1]


def test_chain_sorted__stable_sort():
    object_a = [1]
    object_b = [1]
    object_c = [2]
    object_d = [2]
    seq_a = [object_a, object_c]
    seq_b = [object_b, object_d]

    expected = (object_a, object_b, object_c, object_d)
    result = tuple(utils.chain_sorted(seq_a, seq_b))
    assert expected == result
    assert all(a is b for (a, b) in zip(expected, result))

    expected = (object_b, object_a, object_d, object_c)
    result = tuple(utils.chain_sorted(seq_b, seq_a))
    assert expected == result
    assert all(a is b for (a, b) in zip(expected, result))


def test_chain_sorted__runs_of_values():
    object_a = [1]
    object_b = [1]
    object_c = [2]
    object_d = [2]
    seq_a = [object_a, object_b]
    seq_b = [object_c, object_d]

    expected = (object_a, object_b, object_c, object_d)
    result = tuple(utils.chain_sorted(seq_a, seq_b))
    assert expected == result
    assert all(a is b for (a, b) in zip(expected, result))


def test_chain_sorted__invalid_keywords():
    with pytest.raises(TypeError):
        tuple(utils.chain_sorted((1, 2, 3), foobar=None))


###############################################################################
###############################################################################
# Immutable


def test_immutable__properties_set():
    class ImmutableCls(utils.Immutable):
        def __init__(self, value):
            utils.Immutable.__init__(self, value=value)

    obj = ImmutableCls(17)
    assert obj.value == 17


@pytest.mark.parametrize("key, value", (("value", 13), ("new_value", "foo")))
def test_immutable__properties_immutable(key, value):
    class ImmutableCls(utils.Immutable):
        def __init__(self, value):
            utils.Immutable.__init__(self, value=value)

    obj = ImmutableCls(17)
    with pytest.raises(NotImplementedError):
        setattr(obj, key, value)


def test_immutable__properties_del():
    class ImmutableCls(utils.Immutable):
        def __init__(self, value):
            utils.Immutable.__init__(self, value=value)

    def _del_property(obj):
        del obj.value

    obj = ImmutableCls(17)
    with pytest.raises(NotImplementedError):
        _del_property(obj)


###############################################################################
###############################################################################
# TotallyOrdered


class SomethingOrdered(utils.TotallyOrdered):
    def __init__(self, value):
        self.value = value

    def __lt__(self, other):
        if not isinstance(other, SomethingOrdered):
            return NotImplemented

        return self.value < other.value


# Less than
def test_totally_ordered__lt_vs_lt():
    assert SomethingOrdered(1) < SomethingOrdered(2)


def test_totally_ordered__lt_vs_gt():
    assert not (SomethingOrdered(1) < SomethingOrdered(0))


def test_totally_ordered__lt_vs_wrong_type():
    assert SomethingOrdered(1).__lt__("Foo") == NotImplemented


# Less than or equal
def test_totally_ordered__le_vs_le():
    assert SomethingOrdered(1) <= SomethingOrdered(1)


def test_totally_ordered__le_vs_gt():
    assert not (SomethingOrdered(1) <= SomethingOrdered(0))


def test_totally_ordered__le_vs_wrong_type():
    assert SomethingOrdered(1).__le__("Foo") == NotImplemented


# Greater than or equal
def test_totally_ordered__ge_vs_ge():
    assert SomethingOrdered(1) >= SomethingOrdered(1)


def test_totally_ordered__ge_vs_lt():
    assert not (SomethingOrdered(0) >= SomethingOrdered(1))


def test_totally_ordered__ge_vs_wrong_type():
    assert SomethingOrdered(1).__ge__("Foo") == NotImplemented


# Greater than
def test_totally_ordered__gt_vs_gt():
    assert SomethingOrdered(1) > SomethingOrdered(0)


def test_totally_ordered__gt_vs_eq():
    assert not (SomethingOrdered(0) > SomethingOrdered(0))


def test_totally_ordered__gt_vs_wrong_type():
    assert SomethingOrdered(1).__gt__("Foo") == NotImplemented


# Equal to
def test_totally_ordered__eq_vs_eq():
    assert SomethingOrdered(1) == SomethingOrdered(1)


def test_totally_ordered__eq_vs_ne():
    assert not (SomethingOrdered(1) == SomethingOrdered(2))


def test_totally_ordered__eq_vs_wrong_type():
    assert SomethingOrdered(1).__eq__("Foo") == NotImplemented


# Not equal to
def test_totally_ordered__ne_vs_ne():
    assert SomethingOrdered(1) != SomethingOrdered(2)


def test_totally_ordered__ne_vs_eq():
    assert not (SomethingOrdered(1) != SomethingOrdered(1))


def test_totally_ordered__ne_vs_wrong_type():
    assert SomethingOrdered(1).__ne__("Foo") == NotImplemented


class SomethingBadlyOrdered(utils.TotallyOrdered):
    def __init__(self, value):
        self.value = value


def test_totally_ordered__missing_implementation():
    obj_a = SomethingBadlyOrdered(1)
    obj_b = SomethingBadlyOrdered(2)
    with pytest.raises(NotImplementedError):
        operator.gt(obj_a, obj_b)
