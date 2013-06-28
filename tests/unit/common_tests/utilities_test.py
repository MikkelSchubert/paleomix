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
import random
import pickle
import nose.tools
from nose.tools import assert_equal

import pypeline.common.utilities as utils


################################################################################
################################################################################
## Tests for 'safe_coerce_to_tuple'

def test_safe_coerce_to_tuple__str():
    assert_equal(utils.safe_coerce_to_tuple("foo"), ("foo",))

def test_safe_coerce_to_tuple__unicode():
    assert_equal(utils.safe_coerce_to_tuple(u"foo"), (u"foo",))

def test_safe_coerce_to_tuple__int():
    assert_equal(utils.safe_coerce_to_tuple(17), (17,))

def test_safe_coerce_to_tuple__list():
    assert_equal(utils.safe_coerce_to_tuple([1, 3, 2]), (1, 3, 2))

def test_safe_coerce_to_tuple__tuple():
    assert_equal(utils.safe_coerce_to_tuple((1, 3, 2)), (1, 3, 2))

def test_safe_coerce_to_tuple__iterable():
    assert_equal(utils.safe_coerce_to_tuple(xrange(3)), (0, 1, 2))

def test_safe_coerce_to_tuple__dict():
    assert_equal(utils.safe_coerce_to_tuple({1 : 2, 3 : 4}), ({1 : 2, 3 : 4},))



################################################################################
################################################################################
## Tests for 'safe_coerce_to_frozenset'

def test_safe_coerce_to_frozenset__str():
    assert_equal(utils.safe_coerce_to_frozenset("foo"), frozenset(("foo",)))

def test_safe_coerce_to_frozenset__unicode():
    assert_equal(utils.safe_coerce_to_frozenset(u"foo"), frozenset((u"foo",)))

def test_safe_coerce_to_frozenset__int():
    assert_equal(utils.safe_coerce_to_frozenset(17), frozenset((17,)))

def test_safe_coerce_to_frozenset__list():
    assert_equal(utils.safe_coerce_to_frozenset([1, 3, 2]), frozenset((1, 3, 2)))

def test_safe_coerce_to_frozenset__tuple():
    assert_equal(utils.safe_coerce_to_frozenset((1, 3, 2)), frozenset(((1, 3, 2))))

def test_safe_coerce_to_frozenset__iterable():
    assert_equal(utils.safe_coerce_to_frozenset(xrange(3)), frozenset((0, 1, 2)))

@nose.tools.raises(TypeError)
def test_safe_coerce_to_frozenset__dict():
    utils.safe_coerce_to_frozenset({1 : 2, 3 : 4})



################################################################################
################################################################################
## Tests for 'try_cast'

def test_try_cast__int_to_int():
    assert_equal(utils.try_cast(17, int), 17)

def test_try_cast__float_to_int():
    assert_equal(utils.try_cast(17.3, int), 17)

def test_try_cast__good_str_to_int():
    assert_equal(utils.try_cast("17", int), 17)

def test_try_cast__bad_str_to_int():
    assert_equal(utils.try_cast("x17", int), "x17")

def test_try_cast__list_to_int():
    assert_equal(utils.try_cast([1, 2, 3], int), [1, 2, 3])




################################################################################
################################################################################
## Tests for 'crc32'

def test_crc32_is_unsigned():
    # The following is known to produce an negative value for python 2.7.2
    data = "Nobody inspects the spammish repetition"
    assert utils.crc32(data) >= 0


################################################################################
################################################################################
## Tests for 'set_in'

def test_set_in__single_kw_in_empty_dictionary():
    value = {}
    utils.set_in(value, ["Foo"], 17)
    assert_equal(value, {"Foo" : 17})

def test_set_in__two_kws_in_empty_dictionary():
    value = {}
    utils.set_in(value, ["Foo", 13], 17)
    assert_equal(value, {"Foo" : {13: 17}})

def test_set_in__three_kws_in_empty_dictionary():
    value = {}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert_equal(value, {"Foo" : {13: {(1, 2) : 17}}})

def test_set_in__three_kws_in_partial_dictionary():
    value = {"Foo" : {12 : 0 }}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert_equal(value, {"Foo" : {12: 0, 13: {(1, 2) : 17}}})

    value = {"Foo" : {13 : {"Bar" : None }}}
    utils.set_in(value, ["Foo", 13, (1, 2)], 17)
    assert_equal(value, {"Foo" : {13: {(1, 2) : 17, "Bar" : None}}})

def test_set_in__update_value_one_kw():
    value = {1 : None}
    utils.set_in(value, [1], 3.14)
    assert_equal(value, {1 : 3.14})

def test_set_in__update_value_two_kw():
    value = {1 : {2 : 3}}
    utils.set_in(value, [1, 2], 365)
    assert_equal(value, {1 : {2 : 365}})

@nose.tools.raises(ValueError)
def test_set_in__fail_on_no_kws():
    utils.set_in({}, [], 17)

@nose.tools.raises(TypeError)
def test_set_in__fail_on_invalid_sub_dictionary_first_level():
    utils.set_in(None, [1], 17)

@nose.tools.raises(TypeError)
def test_set_in__fail_on_invalid_sub_dictionary_second_level():
    utils.set_in({1 : None}, [1, 2], 17)

@nose.tools.raises(TypeError)
def test_set_in__fail_on_invalid_sub_dictionary_third_level():
    utils.set_in({1 : {2 : None}}, [1, 2, 3], 17)

def test_set_in__iteratable_keywords():
    value = {}
    utils.set_in(value, iter(["Foo", 13, (1, 2)]), 17)
    assert_equal(value, {"Foo" : {13: {(1, 2) : 17}}})




################################################################################
################################################################################
## Tests for 'get_in'

def test_get_in__get_value_one_keyword():
    assert_equal(utils.get_in({1 : 2}, [1]), 2)

def test_get_in__get_value_two_keywords():
    assert_equal(utils.get_in({1 : {2 : 3}}, [1, 2]), 3)

def test_get_in__get_value_three_keywords():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1, 2, 3]), 4)

def test_get_in__get_default_one_keyword():
    assert_equal(utils.get_in({1 : 2}, [2]), None)

def test_get_in__get_default_one_keyword_with_default():
    assert_equal(utils.get_in({1 : 2}, [2], "other"), "other")

def test_get_in__get_default_three_keywords_fail_at_first():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [2, 2, 4]), None)

def test_get_in__get_default_three_keywords_fail_at_first_with_default():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [2, 2, 4], "other"), "other")

def test_get_in__get_default_three_keywords_fail_at_second():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1, 3, 4]), None)

def test_get_in__get_default_three_keywords_fail_at_second_with_default():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1, 3, 4], "other"), "other")

def test_get_in__get_default_three_keywords_fail_at_third():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1, 2, 4]), None)

def test_get_in__get_default_three_keywords_fail_at_third_with_default():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1, 2, 4], "other"), "other")

def test_get_in__iterator_keywords():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, iter([1, 2, 3])), 4)




################################################################################
################################################################################
## Tests for 'split_before'

def _do_split(lst, key):
    # Convertion to list allows the implementation to be
    # lazy, while making comparisons for asserts easier
    return list(utils.split_before(lst, key))

def test_split_before__split_empty_list():
    assert_equal(_do_split([], None), [])

def test_split_before__split_list_with_no_true_pred():
    assert_equal(_do_split(range(10), lambda x: False), [range(10)])

def test_split_before__split_list_true_pred_at_first_position():
    assert_equal(_do_split(range(4), lambda x: x % 2 == 0), [[0, 1], [2, 3]])

def test_split_before__split_list_true_pred_at_second_position():
    assert_equal(_do_split(range(4), lambda x: x % 2 == 1), [[0], [1, 2], [3]])

def test_split_before__split_consequtive_true_pred():
    assert_equal(_do_split(range(0, 5, 2), lambda x: x % 2 == 0), [[0], [2], [4]])

def test_split_before__no_hits():
    assert_equal(_do_split(range(1, 5), lambda x: x % 5 == 0), [range(1, 5)])




################################################################################
################################################################################
## Tests for 'is_strictly_increasing'

def test_is_strictly_increasing__increasing_sequence():
    assert utils.is_strictly_increasing(range(100))

def test_is_strictly_increasing__non_increasing_sequence():
    lst = range(100)
    first, second = random.sample(lst, 2)
    lst[first], lst[second] = lst[second], lst[first]

    assert not utils.is_strictly_increasing(lst)




################################################################################
################################################################################
## Tests for 'grouper'

def test_grouper__empty_list():
    result = utils.grouper(3, [])
    assert_equal(list(result), [])

def test_grouper__non_empty_list():
    result = utils.grouper(3, range(6))
    expected = [(0, 1, 2), (3, 4, 5)]
    assert_equal(list(result), expected)

def test_grouper__non_empty_list_with_trailing():
    result = utils.grouper(3, range(7))
    expected = [(0, 1, 2), (3, 4, 5), (6, None, None)]
    assert_equal(list(result), expected)

def test_grouper__non_empty_list_with_trailing_fill_value():
    result = utils.grouper(3, range(7), fillvalue = r'\0')
    expected = [(0, 1, 2), (3, 4, 5), (6, r'\0', r'\0')]
    assert_equal(list(result), expected)




################################################################################
################################################################################
## Tests for 'group_by_pred'

def test_group_by_pred__empty_list():
    assert_equal(utils.group_by_pred(id, []), ([], []))

def test_group_by_pred__always_false():
    assert_equal(utils.group_by_pred(lambda x: False, [1, 2, 3]), ([], [1, 2, 3]))

def test_group_by_pred__always_true():
    assert_equal(utils.group_by_pred(lambda x: True, [1, 2, 3]), ([1, 2, 3], []))

def test_group_by_pred__is_even():
    assert_equal(utils.group_by_pred(lambda x: x % 2 == 0, [1, 2, 3]), ([2], [1, 3]))

def test_group_by_pred__iterable():
    assert_equal(utils.group_by_pred(lambda x: x % 2 == 0, xrange(1, 4)), ([2], [1, 3]))




################################################################################
################################################################################
## Tests for 'fragment'

def test_fragment__empty():
    assert_equal(list(utils.fragment(5, "")), [])
    assert_equal(list(utils.fragment(5, [])), [])

def test_fragment__partial_fragment():
    assert_equal(list(utils.fragment(3, "ab")), ["ab"])
    assert_equal(list(utils.fragment(3, ["a", "b"])), [["a", "b"]])

def test_fragment__single_fragment():
    assert_equal(list(utils.fragment(3, "abc")), ["abc"])
    assert_equal(list(utils.fragment(3, ["a", "b", "c"])), [["a", "b", "c"]])

def test_fragment__multiple_fragments():
    assert_equal(list(utils.fragment(3, "abcdef")), ["abc", "def"])
    assert_equal(list(utils.fragment(3, list("abcdef"))), [list("abc"), list("def")])

def test_fragment__multiple_fragments_partial():
    assert_equal(list(utils.fragment(3, "abcdefgh")), ["abc", "def", "gh"])
    assert_equal(list(utils.fragment(3, list("abcdefgh"))), [list("abc"), list("def"), list("gh")])

@nose.tools.raises(TypeError)
def test_fragment__iterable():
    list(utils.fragment(3, xrange(6)))

@nose.tools.raises(TypeError)
def test_fragment__set():
    list(utils.fragment(3, set(range(6))))




################################################################################
################################################################################
## Tests for 'cumsum'

def test_cumsum__empty():
    assert_equal(list(utils.cumsum([])), [])

def test_cumsum__integers():
    assert_equal(list(utils.cumsum(range(-4, 5))), [-4, -7, -9, -10, -10, -9, -7, -4, 0])

def test_cumsum__float():
    assert_equal(list(utils.cumsum((1.0, 2.0, 3.0))), [1.0, 3.0, 6.0])

def test_cumsum__initial():
    assert_equal(list(utils.cumsum(range(5), -10)), [-10, -9, -7, -4, 0])




################################################################################
################################################################################
## Tests for 'fast_pickle_test'

def test_fast_pickle_test__picklable():
    utils.fast_pickle_test(1)
    utils.fast_pickle_test({})
    utils.fast_pickle_test(test_cumsum__empty)

@nose.tools.raises(pickle.PicklingError)
def test_fast_pickle_test__unpicklable_1():
    _func = lambda: None # pragma: no coverage
    utils.fast_pickle_test(_func)

@nose.tools.raises(pickle.PicklingError)
def test_fast_pickle_test__unpicklable_2():
    def _func():
        return None # pragma: no coverage
    utils.fast_pickle_test(_func)

