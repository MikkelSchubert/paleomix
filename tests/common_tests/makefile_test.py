# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import pytest

from paleomix.common.makefile import (
    DEFAULT_NOT_SET,
    REQUIRED_VALUE,
    And,
    DeprecatedOption,
    FASTQPath,
    IsBoolean,
    IsDictOf,
    IsFloat,
    IsInt,
    IsListOf,
    IsNone,
    IsStr,
    IsUnsignedInt,
    MakefileError,
    MakefileSpec,
    Not,
    Or,
    PreProcessMakefile,
    RemovedOption,
    SpecPath,
    SpecTree,
    StringEndsWith,
    StringStartsWith,
    ValueIn,
    ValueMissing,
    ValuesIntersect,
    ValuesSubsetOf,
    WithoutDefaults,
    process_makefile,
    read_makefile,
)

# Dummy value for the path parameters
_DUMMY_PATH = ("a", "random", "path")
_DUMMY_PATH_STR = " :: ".join(_DUMMY_PATH)


class Unhashable:
    __hash__ = None  # pyright: ignore[reportAssignmentType]


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


def _common_invalid_values(
    exclude: Sequence[Any] = (),
    extra: Iterable[Any] = (),
) -> list[Any]:
    selection = list(extra)
    for key, value in _COMMON_INVALID_VALUES.items():
        if key not in exclude:
            selection.append(value)

    return selection


###############################################################################
###############################################################################
# MakefileSpec


def test_makefilespec__description_is_set() -> None:
    desc = "a random description"
    spec = MakefileSpec(description=desc)
    assert spec.description == desc


def test_makefilespec__meets_spec_must_be_implemented() -> None:
    spec = MakefileSpec(description="some description")
    with pytest.raises(NotImplementedError):
        spec(_DUMMY_PATH, 1)


###############################################################################
###############################################################################
# IsInt


def test_is_int__accepts_integers() -> None:
    spec = IsInt()
    assert spec(_DUMMY_PATH, 1234) == 1234
    assert spec(_DUMMY_PATH, 0) == 0
    assert spec(_DUMMY_PATH, -1234) == -1234


@pytest.mark.parametrize("value", _common_invalid_values())
def test_is_int__rejects_not_int(value: object) -> None:
    spec = IsInt()
    with pytest.raises(MakefileError, match="Expected value: an integer"):
        spec(_DUMMY_PATH, value)


def test_is_int__default_description() -> None:
    spec = IsInt()
    assert spec.description == "an integer"


def test_is_int__custom_description() -> None:
    custom_desc = "any old integer"
    spec = IsInt(description=custom_desc)
    assert spec.description == custom_desc


def test_is_int__default_not_set() -> None:
    spec = IsInt()
    assert spec.default is DEFAULT_NOT_SET


def test_is_int__default_set__valid_value() -> None:
    spec = IsInt(default=7913)
    assert spec.default == 7913


def test_is_int__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        IsInt(default="abc")


###############################################################################
###############################################################################
# IsUnsignedInt


def test_is_unsigned_int__accepts_non_negative_integers() -> None:
    spec = IsUnsignedInt()
    assert spec(_DUMMY_PATH, 1234) == 1234
    assert spec(_DUMMY_PATH, 0) == 0


@pytest.mark.parametrize("value", _common_invalid_values(extra=(-1,)))
def test_is_unsigned_int__rejects_not_unsigned_int(value: object) -> None:
    spec = IsUnsignedInt()
    with pytest.raises(MakefileError, match="Expected value: an unsigned integer"):
        spec(_DUMMY_PATH, value)


def test_is_unsigned_int__default_description() -> None:
    spec = IsUnsignedInt()
    assert spec.description == "an unsigned integer"


def test_is_unsigned_int__custom_description() -> None:
    custom_desc = "any old unsigned integer"
    spec = IsUnsignedInt(description=custom_desc)
    assert spec.description == custom_desc


def test_is_unsigned_int__default_not_set() -> None:
    spec = IsUnsignedInt()
    assert spec.default is DEFAULT_NOT_SET


def test_is_unsigned_int__default_set__valid_value() -> None:
    spec = IsUnsignedInt(default=7913)
    assert spec.default == 7913


def test_is_unsigned_int__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        IsUnsignedInt(default=-3)


###############################################################################
###############################################################################
# IsFloat


def test_is_float__accepts_float() -> None:
    spec = IsFloat()
    assert spec(_DUMMY_PATH, 1.0) == 1.0


@pytest.mark.parametrize("value", _common_invalid_values(extra=(0,)))
def test_is_float__rejects_not_float(value: object) -> None:
    spec = IsFloat()
    with pytest.raises(MakefileError, match="Expected value: a float"):
        spec(_DUMMY_PATH, value)


def test_is_float__default_description() -> None:
    spec = IsFloat()
    assert spec.description == "a float"


def test_is_float__custom_description() -> None:
    custom_desc = "a floaty, float"
    spec = IsFloat(description=custom_desc)
    assert spec.description == custom_desc


def test_is_float__default_not_set() -> None:
    spec = IsFloat()
    assert spec.default is DEFAULT_NOT_SET


def test_is_float__default_set__valid_value() -> None:
    spec = IsFloat(default=3.14)
    assert spec.default == 3.14


def test_is_float__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        IsFloat(default="abc")


###############################################################################
###############################################################################
# IsBoolean


def test_is_boolean__accepts_boolean() -> None:
    spec = IsBoolean()
    assert spec(_DUMMY_PATH, False) is False
    assert spec(_DUMMY_PATH, True) is True


@pytest.mark.parametrize("value", _common_invalid_values(exclude=(False,), extra=(0,)))
def test_is_boolean__rejects_not_boolean(value: object) -> None:
    spec = IsBoolean()
    with pytest.raises(MakefileError, match="Expected value: a boolean"):
        spec(_DUMMY_PATH, value)


def test_is_boolean__default_description() -> None:
    spec = IsBoolean()
    assert spec.description == "a boolean"


def test_is_boolean__custom_description() -> None:
    custom_desc = "True or False"
    spec = IsBoolean(description=custom_desc)
    assert spec.description == custom_desc


def test_is_boolean__default_not_set() -> None:
    spec = IsBoolean()
    assert spec.default is DEFAULT_NOT_SET


def test_is_boolean__default_set__valid_value() -> None:
    spec = IsBoolean(default=True)
    assert spec.default


def test_is_boolean__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        IsBoolean(default="abc")


###############################################################################
###############################################################################
# IsStr


def test_is_str__accepts_standard_str() -> None:
    spec = IsStr()
    assert spec(_DUMMY_PATH, "abc") == "abc"


def test_is_str__rejects_empty_str() -> None:
    with pytest.raises(MakefileError):
        IsStr()(_DUMMY_PATH, "")
    with pytest.raises(MakefileError):
        IsStr(min_len=1)(_DUMMY_PATH, "")
    with pytest.raises(MakefileError):
        IsStr(min_len=2)(_DUMMY_PATH, "")
    with pytest.raises(MakefileError):
        IsStr(min_len=3)(_DUMMY_PATH, "")


def test_is_str__accepts_empty_str() -> None:
    spec = IsStr(min_len=0)
    spec(_DUMMY_PATH, "")


def test_is_str__rejects_negative_min_len() -> None:
    with pytest.raises(ValueError, match="min_len must be non-negative"):
        IsStr(min_len=-1)


@pytest.mark.parametrize("value", _common_invalid_values(extra=(1,)))
def test_is_str__rejects_not_str(value: object) -> None:
    spec = IsStr()
    with pytest.raises(MakefileError, match="Expected value: a non-empty string"):
        spec(_DUMMY_PATH, value)


def test_is_str__default_description() -> None:
    spec = IsStr()
    assert spec.description == "a non-empty string"


def test_is_str__custom_description() -> None:
    custom_desc = "a ball of string"
    spec = IsStr(description=custom_desc)
    assert spec.description == custom_desc


def test_is_str__default_not_set() -> None:
    spec = IsStr()
    assert spec.default is DEFAULT_NOT_SET


def test_is_str__default_set__valid_value() -> None:
    spec = IsStr(default="abc")
    assert spec.default == "abc"


def test_is_str__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        IsStr(default=17)


def test_is_str__min_len_0() -> None:
    spec = IsStr(min_len=0)
    spec(_DUMMY_PATH, "")
    spec(_DUMMY_PATH, "a")
    spec(_DUMMY_PATH, "ab")


def test_is_str__min_len_1() -> None:
    spec = IsStr(min_len=1)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "")
    spec(_DUMMY_PATH, "a")
    spec(_DUMMY_PATH, "ab")


def test_is_str__min_len_2() -> None:
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


def test_is_none__accepts_none() -> None:
    spec = IsNone()
    assert spec(_DUMMY_PATH, None) is None


@pytest.mark.parametrize(
    "value", _common_invalid_values(exclude=(None,), extra=(0, ""))
)
def test_is_none__rejects_not_none(value: object) -> None:
    spec = IsNone()
    with pytest.raises(MakefileError, match="Expected value: null or not set"):
        spec(_DUMMY_PATH, value)


def test_is_none__default_description() -> None:
    spec = IsNone()
    assert spec.description == "null or not set"


def test_is_none__custom_description() -> None:
    custom_desc = "NOTHING!"
    spec = IsNone(description=custom_desc)
    assert spec.description == custom_desc


def test_is_none__default_not_set() -> None:
    spec = IsNone()
    assert spec.default is DEFAULT_NOT_SET


def test_is_none__default_not_implemented_for_is_none() -> None:
    with pytest.raises(NotImplementedError):
        IsNone(default=None)


###############################################################################
###############################################################################
# ValueMissing


@pytest.mark.parametrize("value", _COMMON_VALUES)
def test_value_missing(value: object) -> None:
    spec = ValueMissing()

    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# Deprecated option

_DEPRECATED_TMPL = "option has been deprecated and will be removed in the future: %s"
_REMOVED_TMPL = "option has been removed and no longer has any effect: %s"


def test_deprecated_option__properties() -> None:
    spec = DeprecatedOption(IsInt(default=17))

    assert spec.default == 17
    assert spec.description == IsInt(default=17).description


def test_deprecated_option__meets_spec() -> None:
    spec = DeprecatedOption(IsInt(default=17))

    assert spec.meets_spec(10)
    assert not spec.meets_spec("10")


def test_deprecated_option__logs_on_success(caplog: pytest.LogCaptureFixture) -> None:
    expected = _DEPRECATED_TMPL % (_DUMMY_PATH_STR,)
    spec = DeprecatedOption(IsInt(default=17))
    assert spec(_DUMMY_PATH, 10) == 10

    assert expected in caplog.text


def test_deprecated_option__no_log_on_failure(caplog: pytest.LogCaptureFixture) -> None:
    expected = _DEPRECATED_TMPL % (_DUMMY_PATH_STR,)
    spec = DeprecatedOption(IsInt(default=17))

    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "10")

    assert expected not in caplog.text


###############################################################################
###############################################################################
# RemovedOption


@pytest.mark.parametrize("value", _COMMON_VALUES)
def test_deprecated_option(caplog: pytest.LogCaptureFixture, value: object) -> None:
    expected = _REMOVED_TMPL % (_DUMMY_PATH_STR,)
    spec = RemovedOption()

    assert spec(_DUMMY_PATH, value) == value
    assert expected in caplog.text


###############################################################################
###############################################################################
# ValueIn


def test_value_in__single_value_in_set() -> None:
    spec = ValueIn(list(range(5)))
    assert spec(_DUMMY_PATH, 1) == 1


def test_value_in__single_value_not_in_set() -> None:
    spec = ValueIn(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 5)


def test_value_in__case_sensitive__value_in_set() -> None:
    spec = ValueIn(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, "bCe")


def test_value_in__case_sensitive__case_is_normalized() -> None:
    spec = ValueIn(("Abc", "bCe", "cdE"))
    assert spec(_DUMMY_PATH, "Bce") == "bCe"


def test_value_in__default_description() -> None:
    spec = ValueIn(("Abc", "bCe", "cdE"))
    assert spec.description == "value in 'Abc', 'bCe', or 'cdE'"


def test_value_in__custom_description() -> None:
    spec = ValueIn(("Abc", "bCe", "cdE"), description="One of {rvalue}")
    assert spec.description == "One of 'Abc', 'bCe', or 'cdE'"


def test_is_value_in__default_not_set() -> None:
    spec = ValueIn(list(range(5)))
    assert spec.default is DEFAULT_NOT_SET


def test_is_value_in__default_set__valid_value() -> None:
    spec = ValueIn(list(range(5)), default=4)
    assert spec.default == 4


def test_is_value_in__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        ValueIn(list(range(5)), default=5)


@pytest.mark.parametrize("value", _common_invalid_values(extra=("foo",)))
def test_is_value_in__handles_types(value: object) -> None:
    spec = ValueIn((1, 2, 3, 4))
    with pytest.raises(MakefileError, match="Expected value: value in 1, 2, 3, or 4"):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# ValuesIntersects


def test_intersects__single_value_in_set() -> None:
    spec = ValuesIntersect(list(range(5)))
    assert spec(_DUMMY_PATH, [1]) == [1]


def test_intersects__multiple_values_in_set() -> None:
    spec = ValuesIntersect(list(range(5)))
    spec(_DUMMY_PATH, [1, 4])


def test_intersects__single_value_not_in_set() -> None:
    spec = ValuesIntersect(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [5])


def test_intersects__some_values_in_set() -> None:
    spec = ValuesIntersect(list(range(5)))
    spec(_DUMMY_PATH, [4, 5])


def test_intersects__empty_set() -> None:
    spec = ValuesIntersect(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [])


def test_intersects__case_sensitive__value_in_set() -> None:
    spec = ValuesIntersect(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["bCe"])


def test_intersects__case_sensitive__value_in_not_set() -> None:
    spec = ValuesIntersect(("Abc", "bCe", "cdE"))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, ["Bce"])


def test_intersects__chars__case_sensitive() -> None:
    spec = ValuesIntersect("abcdefghijkl")
    spec(_DUMMY_PATH, "a big deal")


def test_intersects__chars__case_sensitive__rejects_differences_in_case() -> None:
    spec = ValuesIntersect("abcdefghijkl")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "A BIG DEAL")


def test_intersects__rejects_dictionary() -> None:
    spec = ValuesIntersect("abc")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {"a": 1, "d": 2})


def test_intersects__default_not_set() -> None:
    spec = ValuesIntersect(list(range(5)))
    assert spec.default is DEFAULT_NOT_SET


def test_intersects__default_set__valid_value() -> None:
    spec = ValuesIntersect(list(range(5)), default=[3, 4])
    assert spec.default == [3, 4]


def test_intersects__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        ValuesIntersect(list(range(5)), default=[5])


@pytest.mark.parametrize("value", _common_invalid_values(extra=("foo",)))
def test_intersects__handles_types(value: object) -> None:
    spec = ValuesIntersect(list(range(5)))
    with pytest.raises(
        MakefileError, match="Expected value: one or more of 0, 1, 2, 3, and 4"
    ):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# ValueSubsetOf


def test_subset_of__single_value_in_set() -> None:
    spec = ValuesSubsetOf(list(range(5)))
    assert spec(_DUMMY_PATH, [1]) == [1]


def test_subset_of__multiple_values_in_set() -> None:
    spec = ValuesSubsetOf(list(range(5)))
    spec(_DUMMY_PATH, [1, 4])


def test_subset_of__empty_set_is_subset() -> None:
    spec = ValuesSubsetOf(list(range(5)))
    spec(_DUMMY_PATH, [])


def test_subset_of__single_value_not_in_set() -> None:
    spec = ValuesSubsetOf(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [5])


def test_subset_of__multiple_values_not_in_set() -> None:
    spec = ValuesSubsetOf(list(range(5)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [4, 5])


def test_subset_of__empty_set() -> None:
    spec = ValuesSubsetOf(list(range(5)))
    spec(_DUMMY_PATH, [])


def test_subset_of__case_sensitive__value_in_set() -> None:
    spec = ValuesSubsetOf(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, ["bCe"])


def test_subset_of__case_sensitive__value_in_not_set() -> None:
    spec = ValuesSubsetOf(("Abc", "bCe", "cdE"))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, ["Bce"])


def test_subset_of__chars__case_sensitive() -> None:
    spec = ValuesSubsetOf("abcdefghijkl ")
    spec(_DUMMY_PATH, "a big deal")


def test_subset_of__chars__case_sensitive__rejects_differences_in_case() -> None:
    spec = ValuesSubsetOf("abcdefghijkl ")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "A big DEAL")


def test_subset_of__rejects_dictionary() -> None:
    spec = ValuesIntersect("abc")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {"a": 1, "b": 2})


def test_subset_of__default_not_set() -> None:
    spec = ValuesSubsetOf(list(range(5)))
    assert spec.default is DEFAULT_NOT_SET


def test_subset_of__default_set__valid_value() -> None:
    spec = ValuesSubsetOf(list(range(5)), default=[3, 4])
    assert spec.default == [3, 4]


def test_subset_of__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        ValuesSubsetOf(list(range(5)), default=[4, 5])


@pytest.mark.parametrize(
    "value", _common_invalid_values(extra=("foo",), exclude=("list_1", ()))
)
def test_subset_of__handles_types(value: object) -> None:
    spec = ValuesSubsetOf(list(range(5)))
    with pytest.raises(
        MakefileError, match="Expected value: subset of 0, 1, 2, 3, and 4"
    ):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# And


def test_and__accepts_when_all_true() -> None:
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    assert spec(_DUMMY_PATH, 0.0) == 0.0


def test_and__rejects_when_first_is_false() -> None:
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 1)


def test_and__rejects_when_second_is_false() -> None:
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 3.0)


def test_and__rejects_when_both_is_false() -> None:
    spec = And(IsFloat, ValueIn((0.0, 1, 2)))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 3)


def test_and__rejects_no_tests() -> None:
    with pytest.raises(ValueError, match="No specification given"):
        And()


def test_and__rejects_non_spec_tests() -> None:
    with pytest.raises(TypeError):
        And(id)  # pyright: ignore[reportArgumentType]


def test_and__default_not_set() -> None:
    spec = And(IsInt, ValueIn(list(range(10))))
    assert spec.default is DEFAULT_NOT_SET


def test_and__default_set__valid_value() -> None:
    spec = And(IsInt, ValueIn(list(range(30))), default=20)
    assert spec.default == 20


def test_and__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        And(IsInt, ValueIn((1,)), default=5)


def test_and__defaults_not_set_in_specs() -> None:
    with pytest.raises(
        ValueError,
        match="Default values cannot be set in specs given to logical operators",
    ):
        And(IsInt(default=10), ValueIn(list(range(100))))


###############################################################################
###############################################################################
# Or


def test_or__accepts_first_test() -> None:
    spec = Or(IsStr, IsBoolean)
    assert spec(_DUMMY_PATH, "Foo") == "Foo"


def test_or__accepts_second_test() -> None:
    spec = Or(IsStr, IsBoolean)
    spec(_DUMMY_PATH, False)


def test_or__rejects_if_both_specs_fail() -> None:
    spec = Or(IsStr, IsBoolean)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 1)


def test_or__rejects_no_tests() -> None:
    with pytest.raises(ValueError, match="No specification given"):
        Or()


def test_or__rejects_non_spec_tests() -> None:
    with pytest.raises(TypeError):
        Or(id)  # pyright: ignore[reportArgumentType]


def test_or__default_not_set() -> None:
    spec = Or(IsInt, ValueIn((10,)))
    assert spec.default is DEFAULT_NOT_SET


def test_or__default_set__valid_value() -> None:
    spec = Or(IsInt, ValueIn((10,)), default=17)
    assert spec.default == 17


def test_or__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        Or(IsInt, ValueIn((10,)), default=5.5)


def test_or__defaults_not_set_in_specs() -> None:
    with pytest.raises(
        ValueError,
        match="Default values cannot be set in specs given to logical operators",
    ):
        Or(IsInt(default=10), ValueIn((10,)))


###############################################################################
###############################################################################
# Not


def test_not__accepts_when_test_is_false() -> None:
    spec = Not(IsInt)
    assert spec(_DUMMY_PATH, True) is True


def test_not__rejects_when_test_is_true() -> None:
    spec = Not(IsInt)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, 1)


def test_not__defaults_not_set_in_specs() -> None:
    with pytest.raises(
        ValueError,
        match="Default values cannot be set in specs given to logical operators",
    ):
        Not(IsInt(default=10))


###############################################################################
###############################################################################
# ValueIn for strings (formerly ValueIn)


def test_string_in__case_sensitive__value_in_set() -> None:
    spec = ValueIn(("Abc", "bCe", "cdE"))
    assert spec(_DUMMY_PATH, "bCe") == "bCe"


def test_string_in__case_insensitive__value_in_set() -> None:
    spec = ValueIn(("Abc", "bCe", "cdE"))
    spec(_DUMMY_PATH, "Bce")


def test_string_in__case_insensitive__value_not_set() -> None:
    spec = ValueIn(("Abc", "bCe", "cdE"))
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "ABce")


def test_string_in__case_insensitive__mixed_string__non_string_found() -> None:
    spec = ValueIn(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, 1)


def test_string_in__case_insensitive__mixed_string__string_found() -> None:
    spec = ValueIn(("A", "c", "B", 1, 2, 3))
    spec(_DUMMY_PATH, "a")


def test_string_in__default_not_set() -> None:
    spec = ValueIn("ABCDEFGH")
    assert spec.default is DEFAULT_NOT_SET


def test_string_in__default_set__valid_value() -> None:
    spec = ValueIn("ABCDEFGH", default="e")
    assert spec.default == "E"


def test_string_in__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        ValueIn("ABCDEFGH", default="i")


@pytest.mark.parametrize("value", _common_invalid_values(extra=("foo",)))
def test_string_in__handles_types(value: object) -> None:
    spec = ValueIn("ABC")
    with pytest.raises(
        MakefileError,
        match="Expected value: value in 'A', 'B', or 'C'",
    ):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# StringStartsWith


def test_string_starts_with__accepts_standard_str() -> None:
    spec = StringStartsWith("A_")
    assert spec(_DUMMY_PATH, "A_BC") == "A_BC"


def test_string_starts_with__rejects_string_without_prefix() -> None:
    spec = StringStartsWith("A_")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "B_GHI")


@pytest.mark.parametrize("value", _common_invalid_values(extra=(1,)))
def test_string_starts_with__rejects_not_uppercase_str(value: object) -> None:
    spec = StringStartsWith("Foo")
    with pytest.raises(
        MakefileError, match="Expected value: a string with prefix 'Foo'"
    ):
        spec(_DUMMY_PATH, value)


def test_string_starts_with__default_not_set() -> None:
    spec = StringStartsWith("Foo")
    assert spec.default is DEFAULT_NOT_SET


def test_string_starts_with__default_set__valid_value() -> None:
    spec = StringStartsWith("Foo", default="FooBar")
    assert spec.default == "FooBar"


def test_string_starts_with__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        StringStartsWith("FooBar", default="BarFoo")


###############################################################################
###############################################################################
# StringEndsWith


def test_string_ends_with__accepts_standard_str() -> None:
    spec = StringEndsWith("_A")
    assert spec(_DUMMY_PATH, "BC_A") == "BC_A"


def test_string_ends_with__rejects_string_without_prefix() -> None:
    spec = StringEndsWith("_A")
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, "GHI_B")


@pytest.mark.parametrize("value", _common_invalid_values(extra=(1,)))
def test_string_ends_with__rejects_not_uppercase_str(value: object) -> None:
    spec = StringEndsWith("Foo")
    with pytest.raises(
        MakefileError, match="Expected value: a string with postfix 'Foo'"
    ):
        spec(_DUMMY_PATH, value)


def test_string_ends_with__default_not_set() -> None:
    spec = StringEndsWith("Foo")
    assert spec.default is DEFAULT_NOT_SET


def test_string_ends_with__default_set__valid_value() -> None:
    spec = StringEndsWith("Bar", default="FooBar")
    assert spec.default == "FooBar"


def test_string_ends_with__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        StringEndsWith("FooBar", default="BarFoo")


###############################################################################
###############################################################################
# FASTQPath


def test_fastq_path__path_must_meet_spec() -> None:
    se_reads = "/path/read.fq.gz"
    pe_reads = "/path/read_{pair}.fq.gz"

    accept_pe = FASTQPath(paired_end=True)
    accept_se = FASTQPath(paired_end=False)
    accept_any = FASTQPath(paired_end=None)

    assert accept_se(_DUMMY_PATH, se_reads) == se_reads
    assert accept_pe(_DUMMY_PATH, pe_reads) == pe_reads
    assert accept_any(_DUMMY_PATH, se_reads) == se_reads
    assert accept_any(_DUMMY_PATH, pe_reads) == pe_reads

    with pytest.raises(MakefileError, match="Expected value: a path without a {pair}"):
        accept_se(_DUMMY_PATH, pe_reads)

    with pytest.raises(MakefileError, match="Expected value: a path with a {pair} key"):
        accept_pe(_DUMMY_PATH, se_reads)


@pytest.mark.parametrize("value", _common_invalid_values(extra=(1,)))
def test_fastq_path__rejects_non_strings(value: object) -> None:
    spec = FASTQPath()
    with pytest.raises(MakefileError, match="Expected value: a path without a {pair}"):
        spec(_DUMMY_PATH, value)


def test_fastq_path__default_not_set() -> None:
    spec = FASTQPath()
    assert spec.default is DEFAULT_NOT_SET


def test_fastq_path__default_set__valid_value() -> None:
    spec = FASTQPath(default="FooBar")
    assert spec.default == "FooBar"


def test_fastq_path__default_set__must_meet_spec() -> None:
    se_reads = "/path/read.fq.gz"
    pe_reads = "/path/read_{pair}.fq.gz"

    FASTQPath(paired_end=None, default=pe_reads)
    FASTQPath(paired_end=None, default=pe_reads)

    FASTQPath(paired_end=True, default=pe_reads)
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        FASTQPath(paired_end=True, default=se_reads)

    FASTQPath(paired_end=False, default=se_reads)
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        FASTQPath(paired_end=False, default=pe_reads)


###############################################################################
###############################################################################
# IsListOf


def test_is_list_of__empty_list_always_ok() -> None:
    spec = IsListOf(IsInt)
    spec(_DUMMY_PATH, [])


def test_is_list_of__list_of_ints_accepted() -> None:
    spec = IsListOf(IsInt)
    spec(_DUMMY_PATH, [1, 2, 3])


def test_is_list_of__list_of_non_ints_rejected() -> None:
    spec = IsListOf(IsInt)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, ["a", "b", "c"])


def test_is_list_of__mixed_list_rejected() -> None:
    spec = IsListOf(IsInt)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [1, "b", 3])


def test_is_list_of__default_description() -> None:
    spec = IsListOf(IsInt, IsFloat)
    assert spec.description == "[(an integer) or (a float), ...]"


def test_is_list_of__non_list_rejected() -> None:
    spec = IsListOf(IsInt)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1: 2})


def test_is_list_of__default_not_set() -> None:
    spec = IsListOf(IsInt)
    assert spec.default is DEFAULT_NOT_SET


def test_is_list_of__default_set__valid_value() -> None:
    spec = IsListOf(IsInt, default=list(range(5)))
    assert spec.default == list(range(5))


def test_is_list_of__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        IsListOf(IsInt, default=17)


def test_is_list_of__defaults_not_set_in_specs() -> None:
    with pytest.raises(
        ValueError,
        match="Default values cannot be set in specs given to logical operators",
    ):
        IsListOf(IsInt(default=10))


###############################################################################
###############################################################################
# IsDictOf


def test_is_dict_of__empty_dict_always_ok() -> None:
    spec = IsDictOf(IsInt, IsStr)
    spec(_DUMMY_PATH, {})


def test_is_dict_of__correct_key_and_value_accepted() -> None:
    spec = IsDictOf(IsInt, IsStr)
    spec(_DUMMY_PATH, {1: "foo"})


def test_is_dict_of__wrong_key_and_correct_value_rejected() -> None:
    spec = IsDictOf(IsInt, IsStr)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1.5: "foo"})


def test_is_dict_of__correct_key_and_wrong_value_rejected() -> None:
    spec = IsDictOf(IsInt, IsStr)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1: 1})


def test_is_dict_of__mixed_rejected() -> None:
    spec = IsDictOf(IsInt, IsStr)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, {1: 1, 2: "foo"})


def test_is_dict_of__default_description() -> None:
    spec = IsDictOf(IsInt, IsStr)
    assert spec.description == "{(an integer) : (a non-empty string)}"


def test_is_dict_of__rejects_non_dict() -> None:
    spec = IsDictOf(IsInt, IsStr)
    with pytest.raises(MakefileError):
        spec(_DUMMY_PATH, [])


def test_is_dict_of__default_not_set() -> None:
    spec = IsDictOf(IsInt, IsInt)
    assert spec.default is DEFAULT_NOT_SET


def test_is_dict_of__default_set__valid_value() -> None:
    spec = IsDictOf(IsInt, IsInt, default={1: 2})
    assert spec.default == {1: 2}


def test_is_dict_of__default_set__must_meet_spec() -> None:
    with pytest.raises(ValueError, match="Default value does not meet requirements"):
        IsDictOf(IsInt, IsInt, default={1: "b"})


def test_is_dict_of__defaults_not_set_in_key_specs() -> None:
    with pytest.raises(ValueError, match="Default values cannot be set"):
        IsDictOf(IsInt(default=10), IsInt)


def test_is_dict_of__defaults_not_set_in_value_specs() -> None:
    with pytest.raises(ValueError, match="Default values cannot be set"):
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
    (ValueIn("abc"), 1),
    (StringStartsWith("FOO"), 1),
    (StringEndsWith("FOO"), 1),
    (IsListOf(IsInt), "foo"),
    (IsDictOf(IsInt, IsInt), 1),
)


@pytest.mark.parametrize(("spec", "value"), _PATH_IN_EXCEPTION_VALUES)
def test_specs__path_is_displayed_in_exception(
    spec: MakefileSpec, value: object
) -> None:
    with pytest.raises(MakefileError, match=_DUMMY_PATH_STR):
        spec(_DUMMY_PATH, value)


###############################################################################
###############################################################################
# process_makefile

_MAKEFILE_SPEC_MET: tuple[tuple[object, SpecTree], ...] = (
    ({"B": 7}, {"A": IsInt, "B": IsInt}),  # String keys
    ({1: "Abc"}, {IsStr: IsInt, IsInt: IsStr}),  # Spec keys
    ({1: "Abc"}, {IsStr(): IsInt, IsInt: IsStr()}),  # Spec keys, instantiated
    ({3: 14}, {IsInt: IsInt, "A": IsInt}),  # Mixed keys, spec matches
    ({"A": 23}, {IsInt: IsInt, "A": IsInt}),  # Mixed keys, key matches
)


@pytest.mark.parametrize(("makefile", "spec"), _MAKEFILE_SPEC_MET)
def test_process_makefile__dict_keys_found(makefile: object, spec: SpecTree) -> None:
    process_makefile(makefile, spec)


_MAKEFILE_SPEC_NOT_MET: tuple[tuple[object, SpecTree], ...] = (
    ({"C": 7}, {"A": IsInt, "B": IsInt}),  # String keys
    ({1.3: "Abc"}, {IsStr: IsInt, IsInt: IsStr}),  # Spec keys
    ({1.3: "Abc"}, {IsStr(): IsInt, IsInt: IsStr()}),  # Spec keys, instantiated
    ({"C": 14}, {IsInt: IsInt, "A": IsInt}),  # Mixed keys, spec matches
    ({"A": 23}, {}),  # Extra keys
)


@pytest.mark.parametrize(("makefile", "spec"), _MAKEFILE_SPEC_NOT_MET)
def test_process_makefile__dict_keys_not_found(
    makefile: object,
    spec: SpecTree,
) -> None:
    with pytest.raises(MakefileError):
        process_makefile(makefile, spec)


def test_validate_makefile__unexpected_type_in_reference() -> None:
    current = {1: 2}
    specs: SpecTree = {IsInt: 2}
    with pytest.raises(TypeError):
        process_makefile(current, specs)


def test_validate_makefile__unexpected_type_in_current() -> None:
    current: object = {1: []}
    specs: SpecTree = {IsInt: {IsInt: IsInt}}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__sets_missing_keys() -> None:
    current = {"A": 1}
    specs: SpecTree = {
        "A": IsInt(default=0),
        "B": IsInt(default=-1),
        "C": IsInt(default=-2),
    }
    expected = {"A": 1, "B": -1, "C": -2}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__mixed_keys() -> None:
    current = {"A": 1}
    specs: SpecTree = {IsStr: IsInt, "B": IsInt(default=-1), "C": IsInt(default=-2)}
    expected = {"A": 1, "B": -1, "C": -2}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__sets_missing_recursive() -> None:
    current = {"A": 1, "B": {"C": 2}}
    specs: SpecTree = {
        "A": IsInt(default=0),
        "B": {"C": IsInt(default=-1), "D": IsInt(default=-2)},
    }
    expected = {"A": 1, "B": {"C": 2, "D": -2}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__sets_missing_recursive__with_missing_substructure() -> None:
    current = {"A": 1}
    specs: SpecTree = {
        "A": IsInt(default=0),
        "B": {"C": IsInt(default=-1), "D": IsInt(default=-2)},
    }
    expected = {"A": 1, "B": {"C": -1, "D": -2}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__shared_subtrees_with_defaults() -> None:
    subtree: SpecTree = {"A": IsInt(default=1234), "B": IsInt(default=5678)}
    specs: SpecTree = {"A": subtree, "B": subtree}
    current = {"A": {"B": 17}, "B": {"A": 71}}
    expected = {"A": {"A": 1234, "B": 17}, "B": {"A": 71, "B": 5678}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__shared_subtrees_with_defaults__defaults_disabled() -> None:
    subtree: SpecTree = {"A": IsInt(default=1234), "B": IsInt(default=5678)}
    specs: SpecTree = {"A": subtree, "B": WithoutDefaults(subtree)}
    current = {"A": {"B": 17}, "B": {"A": 71}}
    expected = {"A": {"A": 1234, "B": 17}, "B": {"A": 71}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__accept_when_required_value_is_set() -> None:
    current = {"A": 1, "B": {"C": 3}}
    expected = {"A": 1, "B": {"C": 3}}
    specs: SpecTree = {"A": IsInt, "B": {"C": IsInt(default=REQUIRED_VALUE)}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__fails_when_required_value_not_set() -> None:
    current = {"A": 1}
    specs: SpecTree = {"A": IsInt, "B": {"C": IsInt(default=REQUIRED_VALUE)}}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__fails_required_value_not_set_in_dynamic_subtree() -> None:
    current: object = {"A": 1, "B": {}}
    specs: SpecTree = {"A": IsInt, IsStr: {"C": IsInt(default=REQUIRED_VALUE)}}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__accept_missing_value_if_in_implicit_subtree() -> None:
    current = {"A": 1}
    expected = {"A": 1}
    specs: SpecTree = {"A": IsInt, IsStr: {"C": IsInt(default=REQUIRED_VALUE)}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__path_shown_in_exception_for_list() -> None:
    with pytest.raises(MakefileError, match=_DUMMY_PATH_STR):
        process_makefile({}, [], _DUMMY_PATH)


def test_process_makefile__path_shown_in_exception_for_dict() -> None:
    with pytest.raises(MakefileError, match=_DUMMY_PATH_STR):
        process_makefile([], {}, _DUMMY_PATH)


def test_process_makefile__implicit_subdict_is_allowed() -> None:
    current = {"A": 1, "B": None}
    expected = {"A": 1, "B": {"C": 3}}
    specs: SpecTree = {"A": IsInt, "B": {"C": IsInt(default=3)}}
    result = process_makefile(current, specs)
    assert result == expected


###############################################################################
###############################################################################
# process_makefile -- lists


def test_process_makefile__list_types_accepted() -> None:
    current = {"A": 1, "B": [17, "Foo"]}
    expected = {"A": 1, "B": [17, "Foo"]}
    specs: SpecTree = {"A": IsInt, "B": [IsInt, IsStr]}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__wrong_list_types() -> None:
    current = {"A": 1, "B": [17, "foo"]}
    specs: SpecTree = {"A": IsInt, "B": [IsInt]}
    with pytest.raises(MakefileError):
        process_makefile(current, specs)


def test_process_makefile__missing_list_defaults_to_empty() -> None:
    current = {"A": 1}
    expected: object = {"A": 1, "B": {"C": []}}
    specs: SpecTree = {"A": IsInt, "B": {"C": [IsInt]}}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__missing_list_default_value() -> None:
    current = {"A": 1}
    expected = {"A": 1, "B": [1, 2, 3]}
    specs: SpecTree = {"A": IsInt, "B": IsListOf(IsInt, default=[1, 2, 3])}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__key_specified_but_no_entries() -> None:
    current = {"A": 1, "B": None}
    expected: object = {"A": 1, "B": []}
    specs: SpecTree = {"A": IsInt, "B": [IsInt]}
    result = process_makefile(current, specs)
    assert result == expected


def test_process_makefile__list_spec_must_contain_specs() -> None:
    specs: SpecTree = {"A": IsInt, "B": [1, 2, 3]}
    with pytest.raises(TypeError):
        process_makefile({}, specs)


def test_process_makefile__list_spec_must_contain_only_specs() -> None:
    specs: SpecTree = {"A": IsInt, "B": [1, 2, IsStr]}
    with pytest.raises(TypeError):
        process_makefile({}, specs)


###############################################################################
###############################################################################
# read_makefile


def test_read_makefile__missing_file() -> None:
    with pytest.raises(FileNotFoundError):
        read_makefile("does_not_exist.yaml", {})


def test_read_makefile__not_a_yaml_file(tmp_path: Path) -> None:
    fpath = tmp_path / "file.fasta"
    fpath.write_text(">Sequence\nACGTTAGATAC\n")

    with pytest.raises(MakefileError):
        read_makefile(fpath, {})


def test_read_makefile__simple_file(tmp_path: Path) -> None:
    fpath = tmp_path / "test.yaml"
    fpath.write_text('Defaults:\n  "First": 1e-4\n  "Second": "a string"\n')

    specs: SpecTree = {"Defaults": {"First": IsFloat, "Second": IsStr}}

    assert read_makefile(fpath, specs) == {
        "Defaults": {"First": 1e-4, "Second": "a string"}
    }


###############################################################################
###############################################################################
# PreProcessMakefile


class _PreProcess(PreProcessMakefile):
    def __call__(self, _path: SpecPath, value: object) -> tuple[object, SpecTree]:
        if isinstance(value, str):
            return int(value), IsInt

        return value, IsInt


def test__preprocess_makefile__missing_value() -> None:
    spec: SpecTree = {"Key": _PreProcess()}
    assert process_makefile({}, spec) == {}


def test__preprocess_makefile__expected_value() -> None:
    spec: SpecTree = {"Key": _PreProcess()}
    assert process_makefile({"Key": 13}, spec) == {"Key": 13}


def test__preprocess_makefile__processed_value() -> None:
    spec: SpecTree = {"Key": _PreProcess()}
    assert process_makefile({"Key": "14"}, spec) == {"Key": 14}


def test__preprocess_makefile__invalid_value() -> None:
    spec: SpecTree = {"Key": _PreProcess()}
    with pytest.raises(MakefileError):
        process_makefile({"Key": False}, spec)


def test__preprocess_makefile__invalid_string() -> None:
    spec: SpecTree = {"Key": _PreProcess()}
    # Failures in processing should propagate out
    with pytest.raises(ValueError, match="invalid literal"):
        process_makefile({"Key": "x14"}, spec)


class _PreProcessWithDefault(PreProcessMakefile):
    def __init__(self, default: object) -> None:
        self._default = default

    def __call__(self, _path: SpecPath, value: object) -> tuple[object, SpecTree]:
        if isinstance(value, str):
            return int(value), IsInt

        return value, IsInt(default=self._default)


def test__preprocess_makefile__with_default__missing_value() -> None:
    spec: SpecTree = {"Key": _PreProcessWithDefault(314)}
    assert process_makefile({}, spec) == {"Key": 314}


def test__preprocess_makefile__with_default__expected_value() -> None:
    spec: SpecTree = {"Key": _PreProcessWithDefault(314)}
    assert process_makefile({"Key": "14"}, spec) == {"Key": 14}
