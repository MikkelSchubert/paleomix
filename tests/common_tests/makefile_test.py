#!/usr/bin/python3
#
# Copyright (c) 2013 Mikkel Schubert <MikkelSch@gmail.com>
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
import pytest

from paleomix.common.makefile import (
    DEFAULT_NOT_SET,
    REQUIRED_VALUE,
    MakefileError,
    MakefileSpec,
    read_makefile,
    process_makefile,
    WithoutDefaults,
    IsInt,
    IsUnsignedInt,
    IsFloat,
    IsBoolean,
    IsStr,
    IsNone,
    ValueIn,
    ValuesIntersect,
    ValuesSubsetOf,
    ValueMissing,
    DeprecatedOption,
    RemovedOption,
    And,
    Or,
    Not,
    StringIn,
    StringStartsWith,
    StringEndsWith,
    IsListOf,
    IsDictOf,
    PreProcessMakefile,
)


# Dummy value for the path parameters
_DUMMY_PATH = ("a", "random", "path")
_DUMMY_PATH_STR = " :: ".join(_DUMMY_PATH)


class Unhashable:
    __hash__ = None


_COMMON_INVALID_VALUES = {
    None: None,
    False: False,
    (): (),
    "list_1": [],
    "dict_1": {},
    "no_hash_1": [Unhashable()],
    "no_hash_2": (Unhashable(),),
    "dict_2": {None: Unhashable()},
    "obj_1": object,
    "obj_2": object(),
}


_COMMON_VALUES = [
    True,
    False,
    None,
    {},
    {"foo": 1},
    [],
    ["label"],
    1.7,
    10,
    "test value",
]


def _common_invalid_values(exclude=(), extra=()):
    selection = list(extra)
    for key, value in _COMMON_INVALID_VALUES.items():
        if key not in exclude:
            selection.append(value)

    return selection


###############################################################################
###############################################################################
# MakefileSpec


def test_makefilespec__description_is_set():
    desc = "a random description"
    spec = MakefileSpec(description=desc)
    assert spec.description == desc


def test_makefilespec__meets_spec_must_be_implemented():
    spec = MakefileSpec(description="some description")
    with pytest.raises(NotImplementedError):
        spec(_DUMMY_PATH, 1)


###############################################################################
###############################################################################
# IsInt


def test_is_int__accepts_integers():
    spec = IsInt()
    spec(_DUMMY_PATH, 1234)
    spec(_DUMMY_PATH, 0)
    spec(_DUMMY_PATH, -1234)


@pytest.mark.parametrize("value", _common_invalid_values())
def test_is_int__rejects_not_int(value):
    spec = IsInt()
    with pytest.raises(MakefileError, match="Expected value: an integer"):
        spec(_DUMMY_PATH, value)


def test_is_int__default_description():
    spec = IsInt()
    assert spec.description == "an integer"


def test_is_int__custom_description():
    custom_desc = "any old integer"
    spec = IsInt(description=custom_desc)
    assert spec.description == custom_desc


def test_is_int__default_not_set():
    spec = IsInt()
    assert spec.default is DEFAULT_NOT_SET


def test_is_int__default_set__valid_value():
    spec = IsInt(default=7913)
    assert spec.default == 7913


def test_is_int__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        IsInt(default="abc")


###############################################################################
###############################################################################
# IsUnsignedInt


def test_is_unsigned_int__accepts_non_negative_integers():
    spec = IsUnsignedInt()
    spec(_DUMMY_PATH, 1234)
    spec(_DUMMY_PATH, 0)


@pytest.mark.parametrize("value", _common_invalid_values(extra=(-1,)))
def test_is_unsigned_int__rejects_not_unsigned_int(value):
    spec = IsUnsignedInt()
    with pytest.raises(MakefileError, match="Expected value: an unsigned integer"):
        spec(_DUMMY_PATH, value)


def test_is_unsigned_int__default_description():
    spec = IsUnsignedInt()
    assert spec.description == "an unsigned integer"


def test_is_unsigned_int__custom_description():
    custom_desc = "any old unsigned integer"
    spec = IsUnsignedInt(description=custom_desc)
    assert spec.description == custom_desc


def test_is_unsigned_int__default_not_set():
    spec = IsUnsignedInt()
    assert spec.default is DEFAULT_NOT_SET


def test_is_unsigned_int__default_set__valid_value():
    spec = IsUnsignedInt(default=7913)
    assert spec.default == 7913


def test_is_unsigned_int__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        IsUnsignedInt(default=-3)


###############################################################################
###############################################################################
# IsFloat


def test_is_float__accepts_float():
    spec = IsFloat()
    spec(_DUMMY_PATH, 1.0)


@pytest.mark.parametrize("value", _common_invalid_values(extra=(0,)))
def test_is_float__rejects_not_float(value):
    spec = IsFloat()
    with pytest.raises(MakefileError, match="Expected value: a float"):
        spec(_DUMMY_PATH, value)


def test_is_float__default_description():
    spec = IsFloat()
    assert spec.description == "a float"


def test_is_float__custom_description():
    custom_desc = "a floaty, float"
    spec = IsFloat(description=custom_desc)
    assert spec.description == custom_desc


def test_is_float__default_not_set():
    spec = IsFloat()
    assert spec.default is DEFAULT_NOT_SET


def test_is_float__default_set__valid_value():
    spec = IsFloat(default=3.14)
    assert spec.default == 3.14


def test_is_float__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        IsFloat(default="abc")


###############################################################################
###############################################################################
# IsBoolean


def test_is_boolean__accepts_boolean():
    spec = IsBoolean()
    spec(_DUMMY_PATH, False)


@pytest.mark.parametrize("value", _common_invalid_values(exclude=(False,), extra=(0,)))
def test_is_boolean__rejects_not_boolean(value):
    spec = IsBoolean()
    with pytest.raises(MakefileError, match="Expected value: a boolean"):
        spec(_DUMMY_PATH, value)


def test_is_boolean__default_description():
    spec = IsBoolean()
    assert spec.description == "a boolean"


def test_is_boolean__custom_description():
    custom_desc = "True or False"
    spec = IsBoolean(description=custom_desc)
    assert spec.description == custom_desc


def test_is_boolean__default_not_set():
    spec = IsBoolean()
    assert spec.default is DEFAULT_NOT_SET


def test_is_boolean__default_set__valid_value():
    spec = IsBoolean(default=True)
    assert spec.default


def test_is_boolean__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        IsBoolean(default="abc")


###############################################################################
###############################################################################
# IsStr


def test_is_str__accepts_standard_str():
    spec = IsStr()
    spec(_DUMMY_PATH, "abc")


def test_is_str__rejects_empty_str():
    with pytest.raises(MakefileError):
        IsStr()(_DUMMY_PATH, "")
    with pytest.raises(MakefileError):
        IsStr(min_len=1)(_DUMMY_PATH, "")
    with pytest.raises(MakefileError):
        IsStr(min_len=2)(_DUMMY_PATH, "")
    with pytest.raises(MakefileError):
        IsStr(min_len=3)(_DUMMY_PATH, "")


def test_is_str__accepts_empty_str():
    spec = IsStr(min_len=0)
    spec(_DUMMY_PATH, "")


def test_is_str__rejects_negative_min_len():
    with pytest.raises(ValueError):
        IsStr(min_len=-1)


@pytest.mark.parametrize("value", _common_invalid_values(extra=(1,)))
def test_is_str__rejects_not_str(value):
    spec = IsStr()
    with pytest.raises(MakefileError, match="Expected value: a non-empty string"):
        spec(_DUMMY_PATH, value)


def test_is_str__default_description():
    spec = IsStr()
    assert spec.description == "a non-empty string"


def test_is_str__custom_description():
    custom_desc = "a ball of string"
    spec = IsStr(description=custom_desc)
    assert spec.description == custom_desc


def test_is_str__default_not_set():
    spec = IsStr()
    assert spec.default is DEFAULT_NOT_SET


def test_is_str__default_set__valid_value():
    spec = IsStr(default="abc")
    assert spec.default == "abc"


def test_is_str__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        IsStr(default=17)


def test_is_str__min_len_0():
    spec = IsStr(min_len=0)
    spec(_DUMMY_PATH, "")
    spec(_DUMMY_PATH, "a")
    spec(_DUMMY_PATH, "ab")


def test_is_str__min_len_1():
    spec = IsStr(min_len=1)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "")
    spec(_DUMMY_PATH, "a")
    spec(_DUMMY_PATH, "ab")


def test_is_str__min_len_2():
    spec = IsStr(min_len=2)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "a")
    spec(_DUMMY_PATH, "ab")
    spec(_DUMMY_PATH, "abc")


###############################################################################
###############################################################################
# IsNone


def test_is_none__accepts_none():
    spec = IsNone()
    spec(_DUMMY_PATH, None)


@pytest.mark.parametrize(
    "value", _common_invalid_values(exclude=(None,), extra=(0, ""))
)
def test_is_none__rejects_not_none(value):
    spec = IsNone()
    with pytest.raises(MakefileError, match="Expected value: null or not set"):
        spec(_DUMMY_PATH, value)


def test_is_none__default_description():
    spec = IsNone()
    assert spec.description == "null or not set"


def test_is_none__custom_description():
    custom_desc = "NOTHING!"
    spec = IsNone(description=custom_desc)
    assert spec.description == custom_desc


def test_is_none__default_not_set():
    spec = IsNone()
    assert spec.default is DEFAULT_NOT_SET


def test_is_none__default_not_implemented_for_is_none():
    with pytest.raises(NotImplementedError):
        IsNone(default=None)


###############################################################################
###############################################################################
# ValueMissing


@pytest.mark.parametrize("value", _COMMON_VALUES)
def test_value_missing(value):
    spec = ValueMissing()

    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# Deprecated option

_DEPRECATED_TMPL = "option has been deprecated and will be removed in the future: %s"


def test_deprecated_option__properties():
    spec = DeprecatedOption(IsInt(default=17))

    assert spec.default == 17
    assert spec.description == IsInt(default=17).description


def test_deprecated_option__meets_spec():
    spec = DeprecatedOption(IsInt(default=17))

    assert spec.meets_spec(10)
    assert not spec.meets_spec("10")


def test_deprecated_option__logs_on_success(caplog):
    expected = _DEPRECATED_TMPL % (_DUMMY_PATH_STR,)
    spec = DeprecatedOption(IsInt(default=17))
    spec(_DUMMY_PATH, 10)

    assert expected in caplog.text


def test_deprecated_option__no_log_on_failure(caplog):
    expected = _DEPRECATED_TMPL % (_DUMMY_PATH_STR,)
    spec = DeprecatedOption(IsInt(default=17))

    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "10")

    assert expected not in caplog.text


###############################################################################
###############################################################################
# RemovedOption


@pytest.mark.parametrize("value", _COMMON_VALUES)
def test_deprecated_option(value):
    spec = RemovedOption()

    spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# ValueIn


def test_value_in__single_value_in_set():
    spec = ValueIn(list(range(5)))
    spec(_DUMMY_PATH, 1)


def test_value_in__single_value_not_in_set():
    spec = ValueIn(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 5)


def test_value_in__single_value_in_set__with_key():
    spec = ValueIn(list(range(5)), key=len)
    spec(_DUMMY_PATH, "a")


def test_value_in__single_value_not_in_set__with_key():
    spec = ValueIn(list(range(5)), key=len)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "abcde")


def test_value_in__case_sensitive__value_in_set():
    spec = ValueIn(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, "bCe")


def test_value_in__case_sensitive__value_in_not_set():
    spec = ValueIn(("Abc", "bCe", "cdE"))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "Bce")


def test_value_in__default_description():
    spec = ValueIn(("Abc", "bCe", "cdE"))
    assert spec.description == "value in 'Abc', 'bCe', or 'cdE'"


def test_value_in__custom_description():
    spec = ValueIn(("Abc", "bCe", "cdE"), description="One of {rvalue}")
    assert spec.description == "One of 'Abc', 'bCe', or 'cdE'"


def test_is_value_in__default_not_set():
    spec = ValueIn(list(range(5)))
    assert spec.default is DEFAULT_NOT_SET


def test_is_value_in__default_set__valid_value():
    spec = ValueIn(list(range(5)), default=4)
    assert spec.default == 4


def test_is_value_in__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        ValueIn(list(range(5)), default=5)


@pytest.mark.parametrize("value", _common_invalid_values(extra=("foo",)))
def test_is_value_in__handles_types(value):
    spec = ValueIn((1, 2, 3, 4))
    with pytest.raises(MakefileError, match="Expected value: value in 1, 2, 3, or 4"):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# ValuesIntersects


def test_intersects__single_value_in_set():
    spec = ValuesIntersect(list(range(5)))
    spec(_DUMMY_PATH, [1])


def test_intersects__multiple_values_in_set():
    spec = ValuesIntersect(list(range(5)))
    spec(_DUMMY_PATH, [1, 4])


def test_intersects__single_value_not_in_set():
    spec = ValuesIntersect(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [5])


def test_intersects__some_values_in_set():
    spec = ValuesIntersect(list(range(5)))
    spec(_DUMMY_PATH, [4, 5])


def test_intersects__empty_set():
    spec = ValuesIntersect(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [])


def test_intersects__case_sensitive__value_in_set():
    spec = ValuesIntersect(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["bCe"])


def test_intersects__case_sensitive__value_in_not_set():
    spec = ValuesIntersect(("Abc", "bCe", "cdE"))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, ["Bce"])


def test_intersects__chars__case_sensitive():
    spec = ValuesIntersect("abcdefghijkl")
    spec(_DUMMY_PATH, "a big deal")


def test_intersects__chars__case_sensitive__rejects_differences_in_case():
    spec = ValuesIntersect("abcdefghijkl")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "A BIG DEAL")


def test_intersects__rejects_dictionary():
    spec = ValuesIntersect("abc")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {"a": 1, "d": 2})


def test_intersects__default_not_set():
    spec = ValuesIntersect(list(range(5)))
    assert spec.default is DEFAULT_NOT_SET


def test_intersects__default_set__valid_value():
    spec = ValuesIntersect(list(range(5)), default=[3, 4])
    assert spec.default == [3, 4]


def test_intersects__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        ValuesIntersect(list(range(5)), default=[5])


@pytest.mark.parametrize("value", _common_invalid_values(extra=("foo",)))
def test_intersects__handles_types(value):
    spec = ValuesIntersect(list(range(5)))
    with pytest.raises(
        MakefileError, match="Expected value: one or more of 0, 1, 2, 3, and 4"
    ):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# ValueSubsetOf


def test_subset_of__single_value_in_set():
    spec = ValuesSubsetOf(list(range(5)))
    spec(_DUMMY_PATH, [1])


def test_subset_of__multiple_values_in_set():
    spec = ValuesSubsetOf(list(range(5)))
    spec(_DUMMY_PATH, [1, 4])


def test_subset_of__empty_set_is_subset():
    spec = ValuesSubsetOf(list(range(5)))
    spec(_DUMMY_PATH, [])


def test_subset_of__single_value_not_in_set():
    spec = ValuesSubsetOf(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [5])


def test_subset_of__multiple_values_not_in_set():
    spec = ValuesSubsetOf(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [4, 5])


def test_subset_of__empty_set():
    spec = ValuesSubsetOf(list(range(5)))
    spec(_DUMMY_PATH, [])


def test_subset_of__case_sensitive__value_in_set():
    spec = ValuesSubsetOf(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["bCe"])


def test_subset_of__case_sensitive__value_in_not_set():
    spec = ValuesSubsetOf(("Abc", "bCe", "cdE"))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, ["Bce"])


def test_subset_of__chars__case_sensitive():
    spec = ValuesSubsetOf("abcdefghijkl ")
    spec(_DUMMY_PATH, "a big deal")


def test_subset_of__chars__case_sensitive__rejects_differences_in_case():
    spec = ValuesSubsetOf("abcdefghijkl ")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "A big DEAL")


def test_subset_of__rejects_dictionary():
    spec = ValuesIntersect("abc")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {"a": 1, "b": 2})


def test_subset_of__default_not_set():
    spec = ValuesSubsetOf(list(range(5)))
    assert spec.default is DEFAULT_NOT_SET


def test_subset_of__default_set__valid_value():
    spec = ValuesSubsetOf(list(range(5)), default=[3, 4])
    assert spec.default == [3, 4]


def test_subset_of__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        ValuesSubsetOf(list(range(5)), default=[4, 5])


@pytest.mark.parametrize(
    "value", _common_invalid_values(extra=("foo",), exclude=("list_1", ()))
)
def test_subset_of__handles_types(value):
    spec = ValuesSubsetOf(list(range(5)))
    with pytest.raises(
        MakefileError, match="Expected value: subset of 0, 1, 2, 3, and 4"
    ):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# And


def test_and__accepts_when_all_true():
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    spec(_DUMMY_PATH, 0.0)


def test_and__rejects_when_first_is_false():
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 1)


def test_and__rejects_when_second_is_false():
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 3.0)


def test_and__rejects_when_both_is_false():
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 3)


def test_and__rejects_no_tests():
    with pytest.raises(ValueError):
        And()


def test_and__rejects_non_spec_tests():
    with pytest.raises(TypeError):
        And(id)


def test_and__default_not_set():
    spec = And(IsInt, ValueIn(list(range(10))))
    assert spec.default is DEFAULT_NOT_SET


def test_and__default_set__valid_value():
    spec = And(IsInt, ValueIn(list(range(30))), default=20)
    assert spec.default == 20


def test_and__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        And(IsInt, ValueIn((1,)), default=5)


def test_and__defaults_not_set_in_specs():
    with pytest.raises(ValueError):
        And(IsInt(default=10), ValueIn((list(range(100)))))


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
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 1)


def test_or__rejects_no_tests():
    with pytest.raises(ValueError):
        Or()


def test_or__rejects_non_spec_tests():
    with pytest.raises(TypeError):
        Or(id)


def test_or__default_not_set():
    spec = Or(IsInt, ValueIn((10,)))
    assert spec.default is DEFAULT_NOT_SET


def test_or__default_set__valid_value():
    spec = Or(IsInt, ValueIn((10,)), default=17)
    assert spec.default == 17


def test_or__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        Or(IsInt, ValueIn((10,)), default=5.5)


def test_or__defaults_not_set_in_specs():
    with pytest.raises(ValueError):
        Or(IsInt(default=10), ValueIn((10,)))


###############################################################################
###############################################################################
# Not


def test_not__accepts_when_test_is_false():
    spec = Not(IsInt)
    spec(_DUMMY_PATH, True)


def test_not__rejects_when_test_is_true():
    spec = Not(IsInt)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 1)


def test_not__defaults_not_set_in_specs():
    with pytest.raises(ValueError):
        Not(IsInt(default=10))


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
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "ABce")


def test_string_in__case_insensitive__mixed_string__non_string_found():
    spec = StringIn(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, 1)


def test_string_in__case_insensitive__mixed_string__string_found():
    spec = StringIn(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, "a")


def test_string_in__default_not_set():
    spec = StringIn("ABCDEFGH")
    assert spec.default is DEFAULT_NOT_SET


def test_string_in__default_set__valid_value():
    spec = StringIn("ABCDEFGH", default="e")
    assert spec.default == "e"


def test_string_in__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        StringIn("ABCDEFGH", default="i")


@pytest.mark.parametrize("value", _common_invalid_values(extra=("foo",)))
def test_string_in__handles_types(value):
    spec = StringIn("ABC")
    with pytest.raises(MakefileError, match="Expected value: one of 'A', 'B', or 'C'"):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# StringStartsWith


def test_string_starts_with__accepts_standard_str():
    spec = StringStartsWith("A_")
    spec(_DUMMY_PATH, "A_BC")


def test_string_starts_with__rejects_string_without_prefix():
    spec = StringStartsWith("A_")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "B_GHI")


@pytest.mark.parametrize("value", _common_invalid_values(extra=(1,)))
def test_string_starts_with__rejects_not_uppercase_str(value):
    spec = StringStartsWith("Foo")
    with pytest.raises(
        MakefileError, match="Expected value: a string with prefix 'Foo'"
    ):
        spec(_DUMMY_PATH, value)


def test_string_starts_with__default_not_set():
    spec = StringStartsWith("Foo")
    assert spec.default is DEFAULT_NOT_SET


def test_string_starts_with__default_set__valid_value():
    spec = StringStartsWith("Foo", default="FooBar")
    assert spec.default == "FooBar"


def test_string_starts_with__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        StringStartsWith("FooBar", default="BarFoo")


###############################################################################
###############################################################################
# StringEndsWith


def test_string_ends_with__accepts_standard_str():
    spec = StringEndsWith("_A")
    spec(_DUMMY_PATH, "BC_A")


def test_string_ends_with__rejects_string_without_prefix():
    spec = StringEndsWith("_A")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "GHI_B")


@pytest.mark.parametrize("value", _common_invalid_values(extra=(1,)))
def test_string_ends_with__rejects_not_uppercase_str(value):
    spec = StringEndsWith("Foo")
    with pytest.raises(
        MakefileError, match="Expected value: a string with postfix 'Foo'"
    ):
        spec(_DUMMY_PATH, value)


def test_string_ends_with__default_not_set():
    spec = StringEndsWith("Foo")
    assert spec.default is DEFAULT_NOT_SET


def test_string_ends_with__default_set__valid_value():
    spec = StringEndsWith("Bar", default="FooBar")
    assert spec.default == "FooBar"


def test_string_ends_with__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        StringEndsWith("FooBar", default="BarFoo")


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
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, ["a", "b", "c"])


def test_is_list_of__mixed_list_rejected():
    spec = IsListOf(IsInt)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [1, "b", 3])


def test_is_list_of__default_description():
    spec = IsListOf(IsInt, IsFloat)
    assert spec.description == "[(an integer) or (a float), ...]"


def test_is_list_of__non_list_rejected():
    spec = IsListOf(IsInt)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1: 2})


def test_is_list_of__default_not_set():
    spec = IsListOf(IsInt)
    assert spec.default is DEFAULT_NOT_SET


def test_is_list_of__default_set__valid_value():
    spec = IsListOf(IsInt, default=list(range(5)))
    assert spec.default == list(range(5))


def test_is_list_of__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        IsListOf(IsInt, default=17)


def test_is_list_of__defaults_not_set_in_specs():
    with pytest.raises(ValueError):
        IsListOf(IsInt(default=10))


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
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1.5: "foo"})


def test_is_dict_of__correct_key_and_wrong_value_rejected():
    spec = IsDictOf(IsInt, IsStr)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1: 1})


def test_is_dict_of__mixed_rejected():
    spec = IsDictOf(IsInt, IsStr)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1: 1, 2: "foo"})


def test_is_dict_of__default_description():
    spec = IsDictOf(IsInt, IsStr)
    assert spec.description == "{(an integer) : (a non-empty string)}"


def test_is_dict_of__rejects_non_dict():
    spec = IsDictOf(IsInt, IsStr)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [])


def test_is_dict_of__default_not_set():
    spec = IsDictOf(IsInt, IsInt)
    assert spec.default is DEFAULT_NOT_SET


def test_is_dict_of__default_set__valid_value():
    spec = IsDictOf(IsInt, IsInt, default={1: 2})
    assert spec.default == {1: 2}


def test_is_dict_of__default_set__must_meet_spec():
    with pytest.raises(ValueError):
        IsDictOf(IsInt, IsInt, default={1: "b"})


def test_is_dict_of__defaults_not_set_in_key_specs():
    with pytest.raises(ValueError):
        IsDictOf(IsInt(default=10), IsInt)


def test_is_dict_of__defaults_not_set_in_value_specs():
    with pytest.raises(ValueError):
        IsDictOf(IsInt, IsInt(default=10))


###############################################################################
###############################################################################
# Path is displayed in exception

_PATH_IN_EXCEPTION_VALUES = (
    (IsInt(), "foo"),
    (IsUnsignedInt(), -1),
    (IsFloat(), "abc"),
    (IsBoolean(), 1),
    (IsStr(), 1),
    (IsNone(), 1),
    (ValueIn([1]), 2),
    (ValuesIntersect([1]), [2]),
    (ValuesSubsetOf([1]), [2]),
    (ValueMissing(), True),
    (And(IsStr), 1),
    (Or(IsStr), 1),
    (Not(IsInt), 1),
    (StringIn("abc"), 1),
    (StringStartsWith("FOO"), 1),
    (StringEndsWith("FOO"), 1),
    (IsListOf(IsInt), "foo"),
    (IsDictOf(IsInt, IsInt), 1),
)


@pytest.mark.parametrize("spec, value", _PATH_IN_EXCEPTION_VALUES)
def test_specs__path_is_displayed_in_exception(spec, value):
    with pytest.raises(MakefileError, match=_DUMMY_PATH_STR):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# process_makefile

_MAKEFILE_SPEC_MET = (
    ({"B": 7}, {"A": IsInt, "B": IsInt}),  # String keys
    ({1: "Abc"}, {IsStr: IsInt, IsInt: IsStr}),  # Spec keys
    ({1: "Abc"}, {IsStr(): IsInt, IsInt: IsStr()}),  # Spec keys, instantiated
    ({3: 14}, {IsInt: IsInt, "A": IsInt}),  # Mixed keys, spec matches
    ({"A": 23}, {IsInt: IsInt, "A": IsInt}),  # Mixed keys, key matches
)


@pytest.mark.parametrize("makefile, spec", _MAKEFILE_SPEC_MET)
def test_process_makefile__dict_keys_found(makefile, spec):
    process_makefile(makefile, spec)


_MAKEFILE_SPEC_NOT_MET = (
    ({"C": 7}, {"A": IsInt, "B": IsInt}),  # String keys
    ({1.3: "Abc"}, {IsStr: IsInt, IsInt: IsStr}),  # Spec keys
    ({1.3: "Abc"}, {IsStr(): IsInt, IsInt: IsStr()}),  # Spec keys, instantiated
    ({"C": 14}, {IsInt: IsInt, "A": IsInt}),  # Mixed keys, spec matches
    ({"A": 23}, {}),  # Extra keys
)


@pytest.mark.parametrize("makefile, spec", _MAKEFILE_SPEC_NOT_MET)
def test_process_makefile__dict_keys_not_found(makefile, spec):
    with pytest.raises(MakefileError):
        process_makefile(makefile, spec)


def test_validate_makefile__unexpected_type_in_reference():
    current = {1: 2}
    specs = {IsInt: 2}
    with pytest.raises(TypeError):
        process_makefile(current, specs)


def test_validate_makefile__unexpected_type_in_current():
    current = {1: []}
    specs = {IsInt: {IsInt: IsInt}}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__sets_missing_keys():
    current = {"A": 1}
    specs = {"A": IsInt(default=0), "B": IsInt(default=-1), "C": IsInt(default=-2)}
    expected = {"A": 1, "B": -1, "C": -2}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__mixed_keys():
    current = {"A": 1}
    specs = {IsStr: IsInt, "B": IsInt(default=-1), "C": IsInt(default=-2)}
    expected = {"A": 1, "B": -1, "C": -2}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__sets_missing_recursive():
    current = {"A": 1, "B": {"C": 2}}
    specs = {
        "A": IsInt(default=0),
        "B": {"C": IsInt(default=-1), "D": IsInt(default=-2)},
    }
    expected = {"A": 1, "B": {"C": 2, "D": -2}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__sets_missing_recursive__with_missing_substructure():
    current = {"A": 1}
    specs = {
        "A": IsInt(default=0),
        "B": {"C": IsInt(default=-1), "D": IsInt(default=-2)},
    }
    expected = {"A": 1, "B": {"C": -1, "D": -2}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__shared_subtrees_with_defaults():
    subtree = {"A": IsInt(default=1234), "B": IsInt(default=5678)}
    specs = {"A": subtree, "B": subtree}
    current = {"A": {"B": 17}, "B": {"A": 71}}
    expected = {"A": {"A": 1234, "B": 17}, "B": {"A": 71, "B": 5678}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__shared_subtrees_with_defaults__defaults_disabled():
    subtree = {"A": IsInt(default=1234), "B": IsInt(default=5678)}
    specs = {"A": subtree, "B": WithoutDefaults(subtree)}
    current = {"A": {"B": 17}, "B": {"A": 71}}
    expected = {"A": {"A": 1234, "B": 17}, "B": {"A": 71}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__accept_when_required_value_is_set():
    current = {"A": 1, "B": {"C": 3}}
    expected = {"A": 1, "B": {"C": 3}}
    specs = {"A": IsInt, "B": {"C": IsInt(default=REQUIRED_VALUE)}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__fails_when_required_value_not_set():
    current = {"A": 1}
    specs = {"A": IsInt, "B": {"C": IsInt(default=REQUIRED_VALUE)}}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__fails_required_value_not_set_in_dynamic_subtree():
    current = {"A": 1, "B": {}}
    specs = {"A": IsInt, IsStr: {"C": IsInt(default=REQUIRED_VALUE)}}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__accept_missing_value_if_in_implicit_subtree():
    current = {"A": 1}
    expected = {"A": 1}
    specs = {"A": IsInt, IsStr: {"C": IsInt(default=REQUIRED_VALUE)}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__path_shown_in_exception_for_list():
    with pytest.raises(MakefileError, match=_DUMMY_PATH_STR):
        process_makefile({}, [], _DUMMY_PATH)


def test_process_makefile__path_shown_in_exception_for_dict():
    with pytest.raises(MakefileError, match=_DUMMY_PATH_STR):
        process_makefile([], {}, _DUMMY_PATH)


def test_process_makefile__implicit_subdict_is_allowed():
    current = {"A": 1, "B": None}
    expected = {"A": 1, "B": {"C": 3}}
    specs = {"A": IsInt, "B": {"C": IsInt(default=3)}}
    result = process_makefile(current, specs)
    assert result == expected


###############################################################################
###############################################################################
# process_makefile -- lists


def test_process_makefile__list_types_accepted():
    current = {"A": 1, "B": [17, "Foo"]}
    expected = {"A": 1, "B": [17, "Foo"]}
    specs = {"A": IsInt, "B": [IsInt, IsStr]}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__wrong_list_types():
    current = {"A": 1, "B": [17, "foo"]}
    specs = {"A": IsInt, "B": [IsInt]}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__missing_list_defaults_to_empty():
    current = {"A": 1}
    expected = {"A": 1, "B": {"C": []}}
    specs = {"A": IsInt, "B": {"C": [IsInt]}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__missing_list_default_value():
    current = {"A": 1}
    expected = {"A": 1, "B": [1, 2, 3]}
    specs = {"A": IsInt, "B": IsListOf(IsInt, default=[1, 2, 3])}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__key_specified_but_no_entries():
    current = {"A": 1, "B": None}
    expected = {"A": 1, "B": []}
    specs = {"A": IsInt, "B": [IsInt]}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__list_spec_must_contain_specs():
    specs = {"A": IsInt, "B": [1, 2, 3]}
    with pytest.raises(TypeError):
        process_makefile({}, specs)


def test_process_makefile__list_spec_must_contain_only_specs():
    specs = {"A": IsInt, "B": [1, 2, IsStr]}
    with pytest.raises(TypeError):
        process_makefile({}, specs)


###############################################################################
###############################################################################
# read_makefile


def test_read_makefile__missing_file():
    with pytest.raises(IOError):
        read_makefile("does_not_exist.yaml", {})


def test_read_makefile__not_a_yaml_file(tmp_path):
    fpath = tmp_path / "file.fasta"
    fpath.write_text(">Sequence\nACGTTAGATAC\n")

    with pytest.raises(MakefileError):
        read_makefile(fpath, {})


def test_read_makefile__simple_file(tmp_path):
    fpath = tmp_path / "test.yaml"
    fpath.write_text('Defaults:\n  "First": 1e-4\n  "Second": "a string"\n')

    specs = {"Defaults": {"First": IsFloat, "Second": IsStr}}

    assert read_makefile(fpath, specs) == {
        "Defaults": {"First": 1e-4, "Second": "a string"}
    }


###############################################################################
###############################################################################
# PreProcessMakefile


class _PreProcess(PreProcessMakefile):
    def __call__(self, path, value):
        if isinstance(value, str):
            return int(value), IsInt

        return value, IsInt


def test__preprocess_makefile__missing_value():
    spec = {"Key": _PreProcess()}
    assert {} == process_makefile({}, spec)


def test__preprocess_makefile__expected_value():
    spec = {"Key": _PreProcess()}
    assert {"Key": 13} == process_makefile({"Key": 13}, spec)


def test__preprocess_makefile__processed_value():
    spec = {"Key": _PreProcess()}
    assert {"Key": 14} == process_makefile({"Key": "14"}, spec)


def test__preprocess_makefile__invalid_value():
    spec = {"Key": _PreProcess()}
    with pytest.raises(MakefileError):
        process_makefile({"Key": False}, spec)


def test__preprocess_makefile__invalid_string():
    spec = {"Key": _PreProcess()}
    # Failures in processing should propagate out
    with pytest.raises(ValueError):
        process_makefile({"Key": "x14"}, spec)


class _PreProcessWithDefault(PreProcessMakefile):
    def __init__(self, default):
        self._default = default

    def __call__(self, path, value):
        if isinstance(value, str):
            return int(value), IsInt

        return value, IsInt(default=self._default)


def test__preprocess_makefile__with_default__missing_value():
    spec = {"Key": _PreProcessWithDefault(314)}
    assert {"Key": 314} == process_makefile({}, spec)


def test__preprocess_makefile__with_default__expected_value():
    spec = {"Key": _PreProcessWithDefault(314)}
    assert {"Key": 14} == process_makefile({"Key": 14}, spec)
