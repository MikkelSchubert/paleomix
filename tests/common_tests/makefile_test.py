#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import os
import types

from nose.tools import \
    assert_is, \
    assert_equal, \
    assert_raises, \
    assert_raises_regexp

from paleomix.common.makefile import \
    DEFAULT_NOT_SET, \
    REQUIRED_VALUE, \
    MakefileError, \
    MakefileSpec, \
    read_makefile, \
    process_makefile, \
    WithoutDefaults, \
    IsInt, \
    IsUnsignedInt, \
    IsFloat, \
    IsBoolean, \
    IsStr, \
    IsNone, \
    ValueLT, \
    ValueLE, \
    ValueGE, \
    ValueGT, \
    ValueIn, \
    ValuesIntersect, \
    ValuesSubsetOf, \
    ValueMissing, \
    And, \
    Or, \
    Xor, \
    Not, \
    StringIn, \
    StringsIntersect, \
    StringsSubsetOf, \
    StringIsUppercase, \
    StringStartsWith, \
    StringEndsWith, \
    IsListOf, \
    IsDictOf, \
    PreProcessMakefile


# Dummy value for the path parameters
_DUMMY_PATH = ("a", "random", "path")
_DUMMY_PATH_STR = ":".join(_DUMMY_PATH)


###############################################################################
###############################################################################
# Setup timestamps for test files

def test_dir():
    return os.path.dirname(os.path.dirname(__file__))


def test_file(*args):
    return os.path.join(test_dir(), "data", *args)


def setup_module():
    timestamps = {test_file("simple.yaml"): 1120719000}

    for filename, timestamp in timestamps.iteritems():
        # Set atime and mtime
        os.utime(filename, (timestamp, timestamp))


###############################################################################
###############################################################################
# MakefileSpec

def test_makefilespec__description_is_set():
    desc = "a random description"
    spec = MakefileSpec(description=desc)
    assert_equal(spec.description, desc)


def test_makefilespec__meets_spec_must_be_implemented():
    spec = MakefileSpec(description="some description")
    assert_raises(NotImplementedError, spec, _DUMMY_PATH, 1)


###############################################################################
###############################################################################
# IsInt

def test_is_int__accepts_integers():
    spec = IsInt()
    spec(_DUMMY_PATH, 1234)
    spec(_DUMMY_PATH, 0)
    spec(_DUMMY_PATH, -1234)


def test_is_int__accepts_longs():
    spec = IsInt()
    spec(_DUMMY_PATH, 1234L)
    spec(_DUMMY_PATH, 0L)
    spec(_DUMMY_PATH, -1234L)


def test_is_int__rejects_not_int():
    def _reject_not_str(value):
        spec = IsInt()
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_str, None
    yield _reject_not_str, False
    yield _reject_not_str, ()
    yield _reject_not_str, {}


def test_is_int__default_description():
    spec = IsInt()
    assert_equal(spec.description, "an integer")


def test_is_int__custom_description():
    custom_desc = "any old integer"
    spec = IsInt(description=custom_desc)
    assert_equal(spec.description, custom_desc)


def test_is_int__default_not_set():
    spec = IsInt()
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_int__default_set__valid_value():
    spec = IsInt(default=7913)
    assert_equal(spec.default, 7913)


def test_is_int__default_set__must_meet_spec():
    assert_raises(ValueError, IsInt, default="abc")


###############################################################################
###############################################################################
# IsUnsignedInt

def test_is_unsigned_int__accepts_non_negative_integers():
    spec = IsUnsignedInt()
    spec(_DUMMY_PATH, 1234)
    spec(_DUMMY_PATH, 0)


def test_is_unsigned_int__accepts_longs():
    spec = IsUnsignedInt()
    spec(_DUMMY_PATH, 1234L)
    spec(_DUMMY_PATH, 0L)


def test_is_unsigned_int__rejects_not_unsigned_int():
    def _reject_not_str(value):
        spec = IsUnsignedInt()
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_str, -1
    yield _reject_not_str, -1L
    yield _reject_not_str, None
    yield _reject_not_str, False
    yield _reject_not_str, ()
    yield _reject_not_str, {}


def test_is_unsigned_int__default_description():
    spec = IsUnsignedInt()
    assert_equal(spec.description, "an unsigned integer")


def test_is_unsigned_int__custom_description():
    custom_desc = "any old unsigned integer"
    spec = IsUnsignedInt(description=custom_desc)
    assert_equal(spec.description, custom_desc)


def test_is_unsigned_int__default_not_set():
    spec = IsUnsignedInt()
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_unsigned_int__default_set__valid_value():
    spec = IsUnsignedInt(default=7913)
    assert_equal(spec.default, 7913)


def test_is_unsigned_int__default_set__must_meet_spec():
    assert_raises(ValueError, IsUnsignedInt, default=-3)


###############################################################################
###############################################################################
# IsFloat

def test_is_float__accepts_float():
    spec = IsFloat()
    spec(_DUMMY_PATH, 1.0)


def test_is_float__rejects_not_float():
    def _reject_not_str(value):
        spec = IsFloat()
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_str, 0
    yield _reject_not_str, None
    yield _reject_not_str, False
    yield _reject_not_str, ()
    yield _reject_not_str, {}


def test_is_float__default_description():
    spec = IsFloat()
    assert_equal(spec.description, "a float")


def test_is_float__custom_description():
    custom_desc = "a floaty, float"
    spec = IsFloat(description=custom_desc)
    assert_equal(spec.description, custom_desc)


def test_is_float__default_not_set():
    spec = IsFloat()
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_float__default_set__valid_value():
    spec = IsFloat(default=3.14)
    assert_equal(spec.default, 3.14)


def test_is_float__default_set__must_meet_spec():
    assert_raises(ValueError, IsFloat, default="abc")


###############################################################################
###############################################################################
# IsBoolean

def test_is_boolean__accepts_boolean():
    spec = IsBoolean()
    spec(_DUMMY_PATH, False)


def test_is_boolean__rejects_not_boolean():
    def _reject_not_str(value):
        spec = IsBoolean()
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_str, None
    yield _reject_not_str, 0
    yield _reject_not_str, ()
    yield _reject_not_str, {}


def test_is_boolean__default_description():
    spec = IsBoolean()
    assert_equal(spec.description, "a boolean")


def test_is_boolean__custom_description():
    custom_desc = "True or False"
    spec = IsBoolean(description=custom_desc)
    assert_equal(spec.description, custom_desc)


def test_is_boolean__default_not_set():
    spec = IsBoolean()
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_boolean__default_set__valid_value():
    spec = IsBoolean(default=True)
    assert_equal(spec.default, True)


def test_is_boolean__default_set__must_meet_spec():
    assert_raises(ValueError, IsBoolean, default="abc")


###############################################################################
###############################################################################
# IsStr

def test_is_str__accepts_standard_str():
    spec = IsStr()
    spec(_DUMMY_PATH, "abc")


def test_is_str__accepts_unicode_str():
    spec = IsStr()
    spec(_DUMMY_PATH, u"def")


def test_is_str__rejects_empty_str():
    spec = IsStr()
    assert_raises(MakefileError, spec, _DUMMY_PATH, "")


def test_is_str__rejects_not_str():
    def _reject_not_str(value):
        spec = IsStr()
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_str, None
    yield _reject_not_str, 1
    yield _reject_not_str, ()
    yield _reject_not_str, {}


def test_is_str__default_description():
    spec = IsStr()
    assert_equal(spec.description, "a non-empty string")


def test_is_str__custom_description():
    custom_desc = "a ball of string"
    spec = IsStr(description=custom_desc)
    assert_equal(spec.description, custom_desc)


def test_is_str__default_not_set():
    spec = IsStr()
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_str__default_set__valid_value():
    spec = IsStr(default="abc")
    assert_equal(spec.default, "abc")


def test_is_str__default_set__must_meet_spec():
    assert_raises(ValueError, IsStr, default=17)


###############################################################################
###############################################################################
# IsNone

def test_is_none__accepts_none():
    spec = IsNone()
    spec(_DUMMY_PATH, None)


def test_is_none__rejects_not_none():
    def _reject_not_none(value):
        spec = IsNone()
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_none, ""
    yield _reject_not_none, 0
    yield _reject_not_none, ()
    yield _reject_not_none, {}


def test_is_none__default_description():
    spec = IsNone()
    assert_equal(spec.description, "None or not set")


def test_is_none__custom_description():
    custom_desc = "NOTHING!"
    spec = IsNone(description=custom_desc)
    assert_equal(spec.description, custom_desc)


def test_is_none__default_not_set():
    spec = IsNone()
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_none__default_not_implemented_for_is_none():
    assert_raises(NotImplementedError, IsNone, default=None)


###############################################################################
###############################################################################
# ValueLT

def test_value_lt__accepts_value_lt():
    spec = ValueLT(7)
    spec(_DUMMY_PATH, 6)


def test_value_lt__rejects_value_eq():
    spec = ValueLT(7)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 7)


def test_value_lt__rejects_value_gt():
    spec = ValueLT(7)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 8)


def test_value_lt__accepts_value_lt__with_key():
    spec = ValueLT(7, key=len)
    spec(_DUMMY_PATH, "abcdef")


def test_value_lt__rejects_value_eq__with_key():
    spec = ValueLT(7, key=len)
    assert_raises(MakefileError, spec, _DUMMY_PATH, "abcdefg")


def test_value_lt__rejects_value_gt__with_key():
    spec = ValueLT(7, key=len)
    assert_raises(MakefileError, spec, _DUMMY_PATH, "abcdefgh")


def test_value_lt__default_description():
    spec = ValueLT('Foo')
    assert_equal(spec.description, "value < 'Foo'")


def test_value_lt__custom_description():
    spec = ValueLT('Bar', description='anything less than {rvalue}')
    assert_equal(spec.description, "anything less than 'Bar'")


def test_is_value_lt__default_not_set():
    spec = ValueLT(10)
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_value_lt__default_set__valid_value():
    spec = ValueLT(10, default=9)
    assert_equal(spec.default, 9)


def test_is_value_lt__default_set__must_meet_spec():
    assert_raises(ValueError, ValueLT, 10, default=17)


###############################################################################
###############################################################################
# ValueLE

def test_value_le__accepts_value_lt():
    spec = ValueLE(7)
    spec(_DUMMY_PATH, 6)


def test_value_le__accepts_value_eq():
    spec = ValueLE(7)
    spec(_DUMMY_PATH, 7)


def test_value_le__rejects_value_gt():
    spec = ValueLE(7)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 8)


def test_value_le__accepts_value_lt__with_key():
    spec = ValueLE(7, key=len)
    spec(_DUMMY_PATH, "abcdef")


def test_value_le__accepts_value_eq__with_key():
    spec = ValueLE(7, key=len)
    spec(_DUMMY_PATH, "abcdefg")


def test_value_le__rejects_value_gt__with_key():
    spec = ValueLE(7, key=len)
    assert_raises(MakefileError, spec, _DUMMY_PATH, "abcdefgh")


def test_value_le__default_description():
    spec = ValueLE('Foo')
    assert_equal(spec.description, "value <= 'Foo'")


def test_value_le__custom_description():
    spec = ValueLE('Bar', description='no more than {rvalue}')
    assert_equal(spec.description, "no more than 'Bar'")


def test_is_value_le__default_not_set():
    spec = ValueLE(10)
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_value_le__default_set__valid_value():
    spec = ValueLE(10, default=10)
    assert_equal(spec.default, 10)


def test_is_value_le__default_set__must_meet_spec():
    assert_raises(ValueError, ValueLE, 10, default=17)


###############################################################################
###############################################################################
# ValueGE

def test_value_ge__rejects_value_lt():
    spec = ValueGE(7)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 6)


def test_value_ge__accepts_value_eq():
    spec = ValueGE(7)
    spec(_DUMMY_PATH, 7)


def test_value_ge__accepts_value_gt():
    spec = ValueGE(7)
    spec(_DUMMY_PATH, 8)


def test_value_ge__accepts_value_lt__with_key():
    spec = ValueGE(7, key=len)
    assert_raises(MakefileError, spec, _DUMMY_PATH, "abcdef")


def test_value_ge__accepts_value_eq__with_key():
    spec = ValueGE(7, key=len)
    spec(_DUMMY_PATH, "abcdefg")


def test_value_ge__accepts_value_gt__with_key():
    spec = ValueGE(7, key=len)
    spec(_DUMMY_PATH, "abcdefgh")


def test_value_ge__default_description():
    spec = ValueGE('Foo')
    assert_equal(spec.description, "value >= 'Foo'")


def test_value_ge__custom_description():
    spec = ValueGE('Bar', description='no less than {rvalue}')
    assert_equal(spec.description, "no less than 'Bar'")


def test_is_value_ge__default_not_set():
    spec = ValueGE(10)
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_value_ge__default_set__valid_value():
    spec = ValueGE(10, default=10)
    assert_equal(spec.default, 10)


def test_is_value_ge__default_set__must_meet_spec():
    assert_raises(ValueError, ValueGE, 10, default=7)


###############################################################################
###############################################################################
# ValueGT

def test_value_gt__rejects_value_lt():
    spec = ValueGT(7)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 6)


def test_value_gt__rejects_value_eq():
    spec = ValueGT(7)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 7)


def test_value_gt__accepts_value_gt():
    spec = ValueGT(7)
    spec(_DUMMY_PATH, 8)


def test_value_gt__accepts_value_lt__with_key():
    spec = ValueGT(7, key=len)
    assert_raises(MakefileError, spec, _DUMMY_PATH, "abcdef")


def test_value_gt__accepts_value_eq__with_key():
    spec = ValueGT(7, key=len)
    assert_raises(MakefileError, spec, _DUMMY_PATH, "abcdefg")


def test_value_gt__accepts_value_gt__with_key():
    spec = ValueGT(7, key=len)
    spec(_DUMMY_PATH, "abcdefgh")


def test_value_gt__default_description():
    spec = ValueGT('Foo')
    assert_equal(spec.description, "value > 'Foo'")


def test_value_gt__custom_description():
    spec = ValueGT('Bar', description='more than {rvalue}')
    assert_equal(spec.description, "more than 'Bar'")


def test_is_value_gt__default_not_set():
    spec = ValueGT(10)
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_value_gt__default_set__valid_value():
    spec = ValueGT(10, default=11)
    assert_equal(spec.default, 11)


def test_is_value_gt__default_set__must_meet_spec():
    assert_raises(ValueError, ValueGT, 10, default=10)


###############################################################################
###############################################################################
# ValueIn

def test_value_in__single_value_in_set():
    spec = ValueIn(range(5))
    spec(_DUMMY_PATH, 1)


def test_value_in__single_value_not_in_set():
    spec = ValueIn(range(5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, 5)


def test_value_in__single_value_in_set__with_key():
    spec = ValueIn(range(5), key=len)
    spec(_DUMMY_PATH, "a")


def test_value_in__single_value_not_in_set__with_key():
    spec = ValueIn(range(5), key=len)
    assert_raises(MakefileError, spec, _DUMMY_PATH, "abcde")


def test_value_in__case_sensitive__value_in_set():
    spec = ValueIn(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, "bCe")


def test_value_in__case_sensitive__value_in_not_set():
    spec = ValueIn(("Abc", "bCe", "cdE"))
    assert_raises(MakefileError, spec, _DUMMY_PATH, "Bce")


def test_value_in__default_description():
    spec = ValueIn(("Abc", "bCe", "cdE"))
    assert_equal(spec.description, "value in 'Abc', 'bCe', or 'cdE'")


def test_value_in__custom_description():
    spec = ValueIn(("Abc", "bCe", "cdE"), description="One of {rvalue}")
    assert_equal(spec.description, "One of 'Abc', 'bCe', or 'cdE'")


def test_is_value_in__default_not_set():
    spec = ValueIn(range(5))
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_value_in__default_set__valid_value():
    spec = ValueIn(range(5), default=4)
    assert_equal(spec.default, 4)


def test_is_value_in__default_set__must_meet_spec():
    assert_raises(ValueError, ValueGT, range(5), default=5)


###############################################################################
###############################################################################
# ValuesIntersects

def test_intersects__single_value_in_set():
    spec = ValuesIntersect(range(5))
    spec(_DUMMY_PATH, [1])


def test_intersects__multiple_values_in_set():
    spec = ValuesIntersect(range(5))
    spec(_DUMMY_PATH, [1, 4])


def test_intersects__single_value_not_in_set():
    spec = ValuesIntersect(range(5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, [5])


def test_intersects__some_values_in_set():
    spec = ValuesIntersect(range(5))
    spec(_DUMMY_PATH, [4, 5])


def test_intersects__empty_set():
    spec = ValuesIntersect(range(5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, [])


def test_intersects__case_sensitive__value_in_set():
    spec = ValuesIntersect(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["bCe"])


def test_intersects__case_sensitive__value_in_not_set():
    spec = ValuesIntersect(("Abc", "bCe", "cdE"))
    assert_raises(MakefileError, spec, _DUMMY_PATH, ["Bce"])


def test_intersects__chars__case_sensitive():
    spec = ValuesIntersect("abcdefghijkl")
    spec(_DUMMY_PATH, "a big deal")


def test_intersects__chars__case_sensitive__rejects_differences_in_case():
    spec = ValuesIntersect("abcdefghijkl")
    assert_raises(MakefileError, spec, _DUMMY_PATH, "A BIG DEAL")


def test_intersects__rejects_dictionary():
    spec = ValuesIntersect("abc")
    assert_raises(MakefileError, spec, _DUMMY_PATH, {"a": 1, "d": 2})


def test_intersects__default_not_set():
    spec = ValuesIntersect(range(5))
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_intersects__default_set__valid_value():
    spec = ValuesIntersect(range(5), default=[3, 4])
    assert_equal(spec.default, [3, 4])


def test_intersects__default_set__must_meet_spec():
    assert_raises(ValueError, ValuesIntersect, range(5), default=[5])


###############################################################################
###############################################################################
# ValueSubsetOf

def test_subset_of__single_value_in_set():
    spec = ValuesSubsetOf(range(5))
    spec(_DUMMY_PATH, [1])


def test_subset_of__multiple_values_in_set():
    spec = ValuesSubsetOf(range(5))
    spec(_DUMMY_PATH, [1, 4])


def test_subset_of__single_value_not_in_set():
    spec = ValuesSubsetOf(range(5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, [5])


def test_subset_of__multiple_values_not_in_set():
    spec = ValuesSubsetOf(range(5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, [4, 5])


def test_subset_of__empty_set():
    spec = ValuesSubsetOf(range(5))
    spec(_DUMMY_PATH, [])


def test_subset_of__case_sensitive__value_in_set():
    spec = ValuesSubsetOf(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["bCe"])


def test_subset_of__case_sensitive__value_in_not_set():
    spec = ValuesSubsetOf(("Abc", "bCe", "cdE"))
    assert_raises(MakefileError, spec, _DUMMY_PATH, ["Bce"])


def test_subset_of__chars__case_sensitive():
    spec = ValuesSubsetOf("abcdefghijkl ")
    spec(_DUMMY_PATH, "a big deal")


def test_subset_of__chars__case_sensitive__rejects_differences_in_case():
    spec = ValuesSubsetOf("abcdefghijkl ")
    assert_raises(MakefileError, spec, _DUMMY_PATH, "A big DEAL")


def test_subset_of__rejects_dictionary():
    spec = ValuesIntersect("abc")
    assert_raises(MakefileError, spec, _DUMMY_PATH, {"a": 1, "b": 2})


def test_subset_of__default_not_set():
    spec = ValuesSubsetOf(range(5))
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_subset_of__default_set__valid_value():
    spec = ValuesSubsetOf(range(5), default=[3, 4])
    assert_equal(spec.default, [3, 4])


def test_subset_of__default_set__must_meet_spec():
    assert_raises(ValueError, ValuesSubsetOf, range(5), default=[4, 5])


###############################################################################
###############################################################################
# And

def test_and__accepts_when_all_true():
    spec = And(IsFloat, ValueLT(1.5))
    spec(_DUMMY_PATH, 0.0)


def test_and__rejects_when_first_is_false():
    spec = And(IsFloat, ValueLT(1.5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, 0)


def test_and__rejects_when_second_is_false():
    spec = And(IsFloat, ValueLT(1.5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, 2.0)


def test_and__rejects_when_both_is_false():
    spec = And(IsFloat, ValueLT(1.5))
    assert_raises(MakefileError, spec, _DUMMY_PATH, 2)


def test_and__rejects_no_tests():
    assert_raises(ValueError, And)


def test_and__rejects_non_spec_tests():
    assert_raises(TypeError, And, id)


def test_and__default_not_set():
    spec = And(IsInt, ValueGT(10))
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_and__default_set__valid_value():
    spec = And(IsInt, ValueGT(10), default=20)
    assert_equal(spec.default, 20)


def test_and__default_set__must_meet_spec():
    assert_raises(ValueError, And, IsInt, ValueGT(10), default=5)


def test_and__defaults_not_set_in_specs():
    assert_raises(ValueError, And, IsInt(default=10), ValueGT(10))


###############################################################################
###############################################################################
# Or

def test_or__accepts_first_test():
    spec = Or(IsStr, IsBoolean)
    spec(_DUMMY_PATH, "Foo")


def test_or__accepts_second_test():
    spec = Or(IsStr, IsBoolean)
    spec(_DUMMY_PATH, False)


def test_or__rejects_if_both_specs_fail():
    spec = Or(IsStr, IsBoolean)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 1)


def test_or__rejects_no_tests():
    assert_raises(ValueError, Or)


def test_or__rejects_non_spec_tests():
    assert_raises(TypeError, Or, id)


def test_or__default_not_set():
    spec = Or(IsInt, ValueGT(10))
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_or__default_set__valid_value():
    spec = Or(IsInt, ValueGT(10), default=17)
    assert_equal(spec.default, 17)


def test_or__default_set__must_meet_spec():
    assert_raises(ValueError, Or, IsInt, ValueGT(10), default=5.5)


def test_or__defaults_not_set_in_specs():
    assert_raises(ValueError, Or, IsInt(default=10), ValueGT(10))


###############################################################################
###############################################################################
# Xor

def test_xor__rejects_when_all_true():
    spec = Xor(IsFloat, ValueLT(1))
    assert_raises(MakefileError, spec, _DUMMY_PATH, 0.0)


def test_xor__accepts_when_first_is_false():
    spec = Xor(IsFloat, ValueLT(1))
    spec(_DUMMY_PATH, 0)


def test_xor__accepts_when_second_is_false():
    spec = Xor(IsFloat, ValueLT(1.0))
    spec(_DUMMY_PATH, 2.0)


def test_xor__rejects_when_both_is_false():
    spec = Xor(IsFloat, ValueLT(1.0))
    assert_raises(MakefileError, spec, _DUMMY_PATH, 2)


def test_xor__rejects_no_tests():
    assert_raises(ValueError, Xor)


def test_xor__rejects_one_test():
    assert_raises(ValueError, Xor, IsInt)


def test_xor__rejects_three_tests():
    assert_raises(ValueError, Xor, IsInt, IsFloat, IsStr)


def test_xor__rejects_non_spec_tests():
    assert_raises(TypeError, Xor, id, id)


def test_xor__default_not_set():
    spec = Xor(IsInt, ValueGT(10))
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_xor__default_set__valid_value():
    spec = Xor(IsInt, ValueGT(10), default=5)
    assert_equal(spec.default, 5)


def test_xor__default_set__must_meet_spec():
    assert_raises(ValueError, Xor, IsInt, ValueGT(10), default=17)


def test_xor__defaults_not_set_in_specs():
    assert_raises(ValueError, Xor, IsInt(default=10), ValueGT(10))


###############################################################################
###############################################################################
# Not

def test_not__accepts_when_test_is_false():
    spec = Not(IsInt)
    spec(_DUMMY_PATH, True)


def test_not__rejects_when_test_is_true():
    spec = Not(IsInt)
    assert_raises(MakefileError, spec, _DUMMY_PATH, 1)


def test_not__defaults_not_set_in_specs():
    assert_raises(ValueError, Not, IsInt(default=10))


###############################################################################
###############################################################################
# StringIn

def test_string_in__case_sensitive__value_in_set():
    spec = StringIn(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, "bCe")


def test_string_in__case_insensitive__value_in_set():
    spec = StringIn(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, "Bce")


def test_string_in__case_insensitive__value_not_set():
    spec = StringIn(("Abc", "bCe", "cdE"))
    assert_raises(MakefileError, spec, _DUMMY_PATH, "ABce")


def test_string_in__case_insensitive__mixed_string__non_string_found():
    spec = StringIn(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, 1)


def test_string_in__case_insensitive__mixed_string__string_found():
    spec = StringIn(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, "a")


def test_string_in__default_not_set():
    spec = StringIn("ABCDEFGH")
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_string_in__default_set__valid_value():
    spec = StringIn("ABCDEFGH", default="e")
    assert_equal(spec.default, "e")


def test_string_in__default_set__must_meet_spec():
    assert_raises(ValueError, StringIn, "ABCDEFGH", default="i")


###############################################################################
###############################################################################
# StringsIntersect

def test_strings_intersect__case_insensitive__value_in_set():
    spec = StringsIntersect(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["Bce"])


def test_strings_intersect__chars__case_insensitive__accepts_differt_in_case():
    spec = StringsIntersect("abcdefghijkl")
    spec(_DUMMY_PATH, "A BIG DEAL")


def test_strings_intersect__case_insensitive__mixed_string__non_string_found():
    spec = StringsIntersect(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, [1])


def test_strings_intersect__case_insensitive__mixed_string_string_found():
    spec = StringsIntersect(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, "a")


def test_strings_intersect__rejects_dictionary():
    spec = StringsIntersect("abc")
    assert_raises(MakefileError, spec, _DUMMY_PATH, {"a": 1, "b": 2})


def test_strings_intersect__default_not_set():
    spec = StringsIntersect("ABCDEFGH")
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_strings_intersect__default_set__valid_value():
    spec = StringsIntersect("ABCDEFGH", default="eabi")
    assert_equal(spec.default, "eabi")


def test_strings_intersect__default_set__must_meet_spec():
    assert_raises(ValueError, StringsIntersect, "ABCDEFGH", default=[1, 2, 3])


###############################################################################
###############################################################################
# StringsSubsetOf

def test_subset_of__case_insensitive__value_in_set():
    spec = StringsSubsetOf(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["Bce"])


def test_subset_of__chars__case_insensitive__accepts_differences_in_case():
    spec = StringsSubsetOf("abcdefghijkl ")
    spec(_DUMMY_PATH, "A big DEAL")


def test_subset_of__case_insensitive__mixed_string__non_string_found():
    spec = StringsSubsetOf(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, [1])


def test_subset_of__case_insensitive__mixed_string_string_found():
    spec = StringsSubsetOf(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, "a")


def test_strings_subset_of__rejects_dictionary():
    spec = StringsSubsetOf("abc")
    assert_raises(MakefileError, spec, _DUMMY_PATH, {"a": 1, "b": 2})


def test_strings_subset_of__default_not_set():
    spec = StringsSubsetOf("ABCDEFGH")
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_strings_subset_of__default_set__valid_value():
    spec = StringsSubsetOf("ABCDEFGH", default="adFg")
    assert_equal(spec.default, "adFg")


def test_string_subset_of__default_set__must_meet_spec():
    assert_raises(ValueError, StringsSubsetOf, "ABCDEFGH", default=[1, 2, 3])


###############################################################################
###############################################################################
# StringIsUppercase

def test_string_is_uppercase__accepts_standard_str():
    spec = StringIsUppercase()
    spec(_DUMMY_PATH, "ABC")


def test_string_is_uppercase__accepts_unicode_str():
    spec = StringIsUppercase()
    spec(_DUMMY_PATH, u"DEF")


def test_string_is_uppercase__rejects_empty_string():
    spec = StringIsUppercase()
    assert_raises(MakefileError, spec, _DUMMY_PATH, "")


def test_string_is_uppercase__rejects_not_uppercase_str():
    def _reject_not_uppercase_str(value):
        spec = StringIsUppercase()
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_uppercase_str, "AcEf"
    yield _reject_not_uppercase_str, None
    yield _reject_not_uppercase_str, 1
    yield _reject_not_uppercase_str, ()
    yield _reject_not_uppercase_str, {}


def test_string_is_uppercase__default_not_set():
    spec = StringIsUppercase()
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_string_is_uppercase__default_set__valid_value():
    spec = StringIsUppercase(default="FOO")
    assert_equal(spec.default, "FOO")


def test_string_is_uppercase__default_set__must_meet_spec():
    assert_raises(ValueError, StringIsUppercase, default="foo")


###############################################################################
###############################################################################
# StringStartsWith

def test_string_starts_with__accepts_standard_str():
    spec = StringStartsWith("A_")
    spec(_DUMMY_PATH, "A_BC")


def test_string_starts_with__accepts_unicode_str():
    spec = StringStartsWith("A_")
    spec(_DUMMY_PATH, u"A_DEF")


def test_string_starts_with__rejects_string_without_prefix():
    spec = StringStartsWith("A_")
    assert_raises(MakefileError, spec, _DUMMY_PATH, "B_GHI")


def test_string_starts_with__rejects_not_uppercase_str():
    spec = StringStartsWith("Foo")

    def _reject_not_str_with_prefix(value):
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_str_with_prefix, None
    yield _reject_not_str_with_prefix, 1
    yield _reject_not_str_with_prefix, ()
    yield _reject_not_str_with_prefix, {}


def test_string_starts_with__default_not_set():
    spec = StringStartsWith("Foo")
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_string_starts_with__default_set__valid_value():
    spec = StringStartsWith("Foo", default="FooBar")
    assert_equal(spec.default, "FooBar")


def test_string_starts_with__default_set__must_meet_spec():
    assert_raises(ValueError, StringStartsWith, "FooBar", default="BarFoo")


###############################################################################
###############################################################################
# StringEndsWith

def test_string_ends_with__accepts_standard_str():
    spec = StringEndsWith("_A")
    spec(_DUMMY_PATH, "BC_A")


def test_string_ends_with__accepts_unicode_str():
    spec = StringEndsWith("_A")
    spec(_DUMMY_PATH, u"DEF_A")


def test_string_ends_with__rejects_string_without_prefix():
    spec = StringEndsWith("_A")
    assert_raises(MakefileError, spec, _DUMMY_PATH, "GHI_B")


def test_string_ends_with__rejects_not_uppercase_str():
    spec = StringEndsWith("Foo")

    def _reject_not_str_with_postfix(value):
        assert_raises(MakefileError, spec, _DUMMY_PATH, value)

    yield _reject_not_str_with_postfix, None
    yield _reject_not_str_with_postfix, 1
    yield _reject_not_str_with_postfix, ()
    yield _reject_not_str_with_postfix, {}


def test_string_ends_with__default_not_set():
    spec = StringEndsWith("Foo")
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_string_ends_with__default_set__valid_value():
    spec = StringEndsWith("Bar", default="FooBar")
    assert_equal(spec.default, "FooBar")


def test_string_ends_with__default_set__must_meet_spec():
    assert_raises(ValueError, StringEndsWith, "FooBar", default="BarFoo")


###############################################################################
###############################################################################
# IsListOf

def test_is_list_of__empty_list_always_ok():
    spec = IsListOf(IsInt)
    spec(_DUMMY_PATH, [])


def test_is_list_of__list_of_ints_accepted():
    spec = IsListOf(IsInt)
    spec(_DUMMY_PATH, [1, 2, 3])


def test_is_list_of__list_of_non_ints_rejected():
    spec = IsListOf(IsInt)
    assert_raises(MakefileError, spec, _DUMMY_PATH, ['a', 'b', 'c'])


def test_is_list_of__mixed_list_rejected():
    spec = IsListOf(IsInt)
    assert_raises(MakefileError, spec, _DUMMY_PATH, [1, 'b', 3])


def test_is_list_of__default_description():
    spec = IsListOf(IsInt, IsFloat)
    assert_equal(spec.description, "[(an integer) or (a float), ...]")


def test_is_list_of__non_list_rejected():
    spec = IsListOf(IsInt)
    assert_raises(MakefileError, spec, _DUMMY_PATH, {1: 2})


def test_is_list_of__default_not_set():
    spec = IsListOf(IsInt)
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_list_of__default_set__valid_value():
    spec = IsListOf(IsInt, default=range(5))
    assert_equal(spec.default, range(5))


def test_is_list_of__default_set__must_meet_spec():
    assert_raises(ValueError, IsListOf, IsInt, default=17)


def test_is_list_of__defaults_not_set_in_specs():
    assert_raises(ValueError, IsListOf, IsInt(default=10))


###############################################################################
###############################################################################
# IsDictOf

def test_is_dict_of__empty_dict_always_ok():
    spec = IsDictOf(IsInt, IsStr)
    spec(_DUMMY_PATH, {})


def test_is_dict_of__correct_key_and_value_accepted():
    spec = IsDictOf(IsInt, IsStr)
    spec(_DUMMY_PATH, {1: "foo"})


def test_is_dict_of__wrong_key_and_correct_value_rejected():
    spec = IsDictOf(IsInt, IsStr)
    assert_raises(MakefileError, spec, _DUMMY_PATH, {1.5: "foo"})


def test_is_dict_of__correct_key_and_wrong_value_rejected():
    spec = IsDictOf(IsInt, IsStr)
    assert_raises(MakefileError, spec, _DUMMY_PATH, {1: 1})


def test_is_dict_of__mixed_rejected():
    spec = IsDictOf(IsInt, IsStr)
    assert_raises(MakefileError, spec, _DUMMY_PATH, {1: 1, 2: "foo"})


def test_is_dict_of__default_description():
    spec = IsDictOf(IsInt, IsStr)
    assert_equal(spec.description, "{(an integer) : (a non-empty string)}")


def test_is_dict_of__rejects_non_dict():
    spec = IsDictOf(IsInt, IsStr)
    assert_raises(MakefileError, spec, _DUMMY_PATH, [])


def test_is_dict_of__default_not_set():
    spec = IsDictOf(IsInt, IsInt)
    assert_is(spec.default, DEFAULT_NOT_SET)


def test_is_dict_of__default_set__valid_value():
    spec = IsDictOf(IsInt, IsInt, default={1: 2})
    assert_equal(spec.default, {1: 2})


def test_is_dict_of__default_set__must_meet_spec():
    assert_raises(ValueError, IsDictOf, IsInt, IsInt, default={1: "b"})


def test_is_dict_of__defaults_not_set_in_key_specs():
    assert_raises(ValueError, IsDictOf, IsInt(default=10), IsInt)


def test_is_dict_of__defaults_not_set_in_value_specs():
    assert_raises(ValueError, IsDictOf, IsInt, IsInt(default=10))


###############################################################################
###############################################################################
# Path is displayed in exception

def test_specs__path_is_displayed_in_exception():
    def _path_is_displayed_in_exception(spec, value):
        assert_raises_regexp(MakefileError, _DUMMY_PATH_STR,
                             spec, _DUMMY_PATH, value)

    yield _path_is_displayed_in_exception, IsInt(), "foo"
    yield _path_is_displayed_in_exception, IsUnsignedInt(), -1
    yield _path_is_displayed_in_exception, IsFloat(), "abc"
    yield _path_is_displayed_in_exception, IsBoolean(), 1
    yield _path_is_displayed_in_exception, IsStr(), 1
    yield _path_is_displayed_in_exception, IsNone(), 1
    yield _path_is_displayed_in_exception, ValueLT(0), 1
    yield _path_is_displayed_in_exception, ValueLE(0), 1
    yield _path_is_displayed_in_exception, ValueGE(0), -1
    yield _path_is_displayed_in_exception, ValueGT(0), -1
    yield _path_is_displayed_in_exception, ValueIn([1]), 2
    yield _path_is_displayed_in_exception, ValuesIntersect([1]), [2]
    yield _path_is_displayed_in_exception, ValuesSubsetOf([1]), [2]
    yield _path_is_displayed_in_exception, ValueMissing(), True
    yield _path_is_displayed_in_exception, And(IsStr), 1
    yield _path_is_displayed_in_exception, Or(IsStr), 1
    yield _path_is_displayed_in_exception, Xor(IsStr, IsInt), True
    yield _path_is_displayed_in_exception, Not(IsInt), 1
    yield _path_is_displayed_in_exception, StringIn("abc"), 1
    yield _path_is_displayed_in_exception, StringsIntersect("abc"), [1]
    yield _path_is_displayed_in_exception, StringsSubsetOf("abc"), [1]
    yield _path_is_displayed_in_exception, StringIsUppercase(), 1
    yield _path_is_displayed_in_exception, StringStartsWith("FOO"), 1
    yield _path_is_displayed_in_exception, StringEndsWith("FOO"), 1
    yield _path_is_displayed_in_exception, IsListOf(IsInt), "foo"
    yield _path_is_displayed_in_exception, IsDictOf(IsInt, IsInt), 1


###############################################################################
###############################################################################
# process_makefile

def test_process_makefile__dict_keys_found():
    def _dict_keys_found(current, specs):
        process_makefile(current, specs)

    # String keys
    yield _dict_keys_found, {"B": 7}, {"A": IsInt, "B": IsInt}
    # Spec keys
    yield _dict_keys_found, {1: "Abc"}, {IsStr: IsInt, IsInt: IsStr}
    # Spec keys, instantiated
    yield _dict_keys_found, {1: "Abc"}, {IsStr(): IsInt, IsInt: IsStr()}
    # Mixed keys, spec matches
    yield _dict_keys_found, {3: 14}, {IsInt: IsInt, "A": IsInt}
    # Mixed keys, key matches
    yield _dict_keys_found, {"A": 23}, {IsInt: IsInt, "A": IsInt}


def test_process_makefile__dict_keys_not_found():
    def _dict_keys_missing(current, specs):
        assert_raises(MakefileError, process_makefile, current, specs)

    # String keys
    yield _dict_keys_missing, {"C": 7}, {"A": IsInt, "B": IsInt}
    # Spec keys
    yield _dict_keys_missing, {1.3: "Abc"}, {IsStr: IsInt, IsInt: IsStr}
    # Spec keys, instantiated
    yield _dict_keys_missing, {1.3: "Abc"}, {IsStr(): IsInt, IsInt: IsStr()}
    # Mixed keys, spec matches
    yield _dict_keys_missing, {"C": 14}, {IsInt: IsInt, "A": IsInt}
    yield _dict_keys_missing, {"A": 23}, {}


def test_validate_makefile__unexpected_type_in_reference():
    current = {1: 2}
    specs = {IsInt: 2}
    assert_raises(TypeError, process_makefile, current, specs)


def test_validate_makefile__unexpected_type_in_current():
    current = {1: []}
    specs = {IsInt: {IsInt: IsInt}}
    assert_raises(MakefileError, process_makefile, current, specs)


def test_process_makefile__sets_missing_keys():
    current = {"A": 1}
    specs = {"A": IsInt(default=0),
             "B": IsInt(default=-1),
             "C": IsInt(default=-2)}
    expected = {"A": 1, "B": -1, "C": -2}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__mixed_keys():
    current = {"A": 1}
    specs = {IsStr: IsInt,
             "B": IsInt(default=-1),
             "C": IsInt(default=-2)}
    expected = {"A": 1, "B": -1, "C": -2}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__sets_missing_recursive():
    current = {"A": 1, "B": {"C":  2}}
    specs = {"A": IsInt(default=0),
             "B": {"C": IsInt(default=-1),
                   "D": IsInt(default=-2)}}
    expected = {"A": 1, "B": {"C":  2, "D": -2}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__sets_missing_recursive__with_missing_substructure():
    current = {"A": 1}
    specs = {"A": IsInt(default=0),
             "B": {"C": IsInt(default=-1),
                   "D": IsInt(default=-2)}}
    expected = {"A": 1, "B": {"C": -1, "D": -2}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__shared_subtrees_with_defaults():
    subtree = {"A": IsInt(default=1234),
               "B": IsInt(default=5678)}
    specs = {"A": subtree,
             "B": subtree}
    current = {"A": {"B": 17},
               "B": {"A": 71}}
    expected = {"A": {"A": 1234, "B": 17},
                "B": {"A": 71, "B": 5678}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__shared_subtrees_with_defaults__defaults_disabled():
    subtree = {"A": IsInt(default=1234),
               "B": IsInt(default=5678)}
    specs = {"A": subtree,
             "B": WithoutDefaults(subtree)}
    current = {"A": {"B": 17},
               "B": {"A": 71}}
    expected = {"A": {"A": 1234, "B": 17},
                "B": {"A": 71}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__accept_when_required_value_is_set():
    current = {"A": 1, "B": {"C": 3}}
    expected = {"A": 1, "B": {"C": 3}}
    specs = {"A": IsInt, "B": {"C": IsInt(default=REQUIRED_VALUE)}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__fails_when_required_value_not_set():
    current = {"A": 1}
    specs = {"A": IsInt, "B": {"C": IsInt(default=REQUIRED_VALUE)}}
    assert_raises(MakefileError, process_makefile, current, specs)


def test_process_makefile__fails_required_value_not_set_in_dynamic_subtree():
    current = {"A": 1, "B": {}}
    specs = {"A": IsInt, IsStr: {"C": IsInt(default=REQUIRED_VALUE)}}
    assert_raises(MakefileError, process_makefile, current, specs)


def test_process_makefile__accept_missing_value_if_in_implicit_subtree():
    current = {"A": 1}
    expected = {"A": 1}
    specs = {"A": IsInt, IsStr: {"C": IsInt(default=REQUIRED_VALUE)}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__path_shown_in_exception_for_list():
    assert_raises_regexp(MakefileError, _DUMMY_PATH_STR,
                         process_makefile, {}, [], _DUMMY_PATH)


def test_process_makefile__path_shown_in_exception_for_dict():
    assert_raises_regexp(MakefileError, _DUMMY_PATH_STR,
                         process_makefile, [], {}, _DUMMY_PATH)


def test_process_makefile__implicit_subdict_is_allowed():
    current = {"A": 1, "B": None}
    expected = {"A": 1, "B": {"C": 3}}
    specs = {"A": IsInt, "B": {"C": IsInt(default=3)}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


###############################################################################
###############################################################################
# process_makefile -- lists

def test_process_makefile__list_types_accepted():
    current = {"A": 1, "B": [17, "Foo"]}
    expected = {"A": 1, "B": [17, "Foo"]}
    specs = {"A": IsInt, "B": [IsInt, IsStr]}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__wrong_list_types():
    current = {"A": 1, "B": [17, "foo"]}
    specs = {"A": IsInt, "B": [IsInt]}
    assert_raises(MakefileError, process_makefile, current, specs)


def test_process_makefile__missing_list_defaults_to_empty():
    current = {"A": 1}
    expected = {"A": 1, "B": {"C": []}}
    specs = {"A": IsInt, "B": {"C": [IsInt]}}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__missing_list_default_value():
    current = {"A": 1}
    expected = {"A": 1, "B": [1, 2, 3]}
    specs = {"A": IsInt, "B": IsListOf(IsInt, default=[1, 2, 3])}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__key_specified_but_no_entries():
    current = {"A": 1, "B": None}
    expected = {"A": 1, "B": []}
    specs = {"A": IsInt, "B": [IsInt]}
    result = process_makefile(current, specs)
    assert_equal(result, expected)


def test_process_makefile__list_spec_must_contain_specs():
    specs = {"A": IsInt, "B": [1, 2, 3]}
    assert_raises(TypeError, process_makefile, {}, specs)


def test_process_makefile__list_spec_must_contain_only_specs():
    specs = {"A": IsInt, "B": [1, 2, IsStr]}
    assert_raises(TypeError, process_makefile, {}, specs)


###############################################################################
###############################################################################
# read_makefile

def test_read_makefile__missing_file():
    assert_raises(IOError, read_makefile, "does_not_exist.yaml", {})


def test_read_makefile__not_a_yaml_file():
    fpath = test_file("fasta_file.fasta")
    assert_raises(MakefileError, read_makefile, fpath, {})


def test_read_makefile__missing_simple_file():
    specs = {"Defaults": {"First": IsFloat, "Second": IsStr}}
    expected = {
        "Makefile": {"Defaults": {"First": 1e-4,
                                  "Second": "a string"}},
        "Statistics": {
            "Filename": test_file("simple.yaml"),
            "Hash": "563a2052b67dcde9f193fbe8d51fa2b6f0806505",
            "MTime": "2005-07-07 08:50:00",
        }
    }
    result = read_makefile(test_file("simple.yaml"), specs)
    assert_equal(expected, result)


###############################################################################
###############################################################################
# PreProcessMakefile

class _PreProcess(PreProcessMakefile):
    def __call__(self, path, value):
        if isinstance(value, types.StringTypes):
            return int(value), IsInt

        return value, IsInt


def test__preprocess_makefile__missing_value():
    spec = {"Key": _PreProcess()}
    assert_equal({}, process_makefile({}, spec))


def test__preprocess_makefile__expected_value():
    spec = {"Key": _PreProcess()}
    assert_equal({"Key": 13}, process_makefile({"Key": 13}, spec))


def test__preprocess_makefile__processed_value():
    spec = {"Key": _PreProcess()}
    assert_equal({"Key": 14}, process_makefile({"Key": "14"}, spec))


def test__preprocess_makefile__invalid_value():
    spec = {"Key": _PreProcess()}
    assert_raises(MakefileError, process_makefile, {"Key": False}, spec)


def test__preprocess_makefile__invalid_string():
    spec = {"Key": _PreProcess()}
    # Failures in processing should propagate out
    assert_raises(ValueError, process_makefile, {"Key": "x14"}, spec)


class _PreProcessWithDefault(PreProcessMakefile):
    def __init__(self, default):
        self._default = default

    def __call__(self, path, value):
        if isinstance(value, types.StringTypes):
            return int(value), IsInt

        return value, IsInt(default=self._default)


def test__preprocess_makefile__with_default__missing_value():
    spec = {"Key": _PreProcessWithDefault(314)}
    assert_equal({"Key": 314}, process_makefile({}, spec))


def test__preprocess_makefile__with_default__expected_value():
    spec = {"Key": _PreProcessWithDefault(314)}
    assert_equal({"Key": 14}, process_makefile({"Key": 14}, spec))
