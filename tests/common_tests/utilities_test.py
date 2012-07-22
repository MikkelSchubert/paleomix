#!/usr/bin/python

import random
import nose.tools
from nose.tools import assert_equal

import pypeline.common.utilities as utils


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
    dd = {}
    utils.set_in(dd, ["Foo"], 17)
    assert_equal(dd, {"Foo" : 17})

def test_set_in__two_kws_in_empty_dictionary():
    dd = {}
    utils.set_in(dd, ["Foo", 13], 17)
    assert_equal(dd, {"Foo" : {13: 17}})

def test_set_in__three_kws_in_empty_dictionary():
    dd = {}
    utils.set_in(dd, ["Foo", 13, (1, 2)], 17)
    assert_equal(dd, {"Foo" : {13: {(1,2) : 17}}})

def test_set_in__three_kws_in_partial_dictionary():
    dd = {"Foo" : {12 : 0 }}
    utils.set_in(dd, ["Foo", 13, (1, 2)], 17)
    assert_equal(dd, {"Foo" : {12: 0, 13: {(1,2) : 17}}})

    dd = {"Foo" : {13 : {"Bar" : None }}}
    utils.set_in(dd, ["Foo", 13, (1, 2)], 17)
    assert_equal(dd, {"Foo" : {13: {(1,2) : 17, "Bar" : None}}})

def test_set_in__update_value_one_kw():
    dd = {1 : None}
    utils.set_in(dd, [1], 3.14)
    assert_equal(dd, {1 : 3.14})

def test_set_in__update_value_two_kw():
    dd = {1 : {2 : 3}}
    utils.set_in(dd, [1, 2], 365)
    assert_equal(dd, {1 : {2 : 365}})       

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
    dd = {}
    utils.set_in(dd, iter(["Foo", 13, (1, 2)]), 17)
    assert_equal(dd, {"Foo" : {13: {(1,2) : 17}}})

        


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
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [2,2,4]), None)

def test_get_in__get_default_three_keywords_fail_at_first_with_default():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [2,2,4], "other"), "other")

def test_get_in__get_default_three_keywords_fail_at_second():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1,3,4]), None)

def test_get_in__get_default_three_keywords_fail_at_second_with_default():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1,3,4], "other"), "other")

def test_get_in__get_default_three_keywords_fail_at_third():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1,2,4]), None)

def test_get_in__get_default_three_keywords_fail_at_third_with_default():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, [1,2,4], "other"), "other")

def test_get_in__iterator_keywords():
    assert_equal(utils.get_in({1 : {2 : {3 : 4}}}, iter([1, 2, 3])), 4)




################################################################################
################################################################################
## Tests for 'split_before'

def _do_split(lst, key):
    # Convertion to list allows the implementation to be
    # lazy, and makes comparisons easier
    return map(list, utils.split_before(lst, key))

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




################################################################################
################################################################################
## Tests for 'is_strictly_increasing'

def test_is_strictly_increasing__increasing_sequence():
    assert utils.is_strictly_increasing(range(100))

def test_is_strictly_increasing__non_increasing_sequence():
    lst = range(100)
    a, b = random.sample(lst, 2)
    lst[a], lst[b] = lst[b], lst[a]

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
    result = utils.grouper(3, range(7), fillvalue = '\0')
    expected = [(0, 1, 2), (3, 4, 5), (6, '\0', '\0')]
    assert_equal(list(result), expected)


       
################################################################################
################################################################################
## Tests for 'group_by_pred'

def test_group_by_pred__empty_list():
    assert_equal(utils.group_by_pred(id, []), ([], []))

def test_group_by_pred__always_false():
    assert_equal(utils.group_by_pred(lambda x: False, [1,2,3]), ([], [1,2,3]))

def test_group_by_pred__always_true():
    assert_equal(utils.group_by_pred(lambda x: True, [1,2,3]), ([1,2,3], []))

def test_group_by_pred__is_even():
    assert_equal(utils.group_by_pred(lambda x: x % 2 == 0, [1,2,3]), ([2], [1,3]))

def test_group_by_pred__iterable():
    assert_equal(utils.group_by_pred(lambda x: x % 2 == 0, xrange(1,4)), ([2], [1,3]))
