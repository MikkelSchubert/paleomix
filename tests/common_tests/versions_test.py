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
import operator
import pickle
import sys

import pytest

import paleomix.common.versions as versions

###############################################################################
###############################################################################
# Check class


def test_check__func_must_be_callable():
    with pytest.raises(TypeError):
        versions.Check("FooBar", 3, 7, 5)


def test_check_str():
    obj = versions.Check("FooBar", operator.lt, 3, 7, 5)
    assert str(obj) == "FooBar"


###############################################################################
###############################################################################
# Check class -- hash and comparisons


def test_check__eq_same_func_desc_and_version():
    obj_1 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    obj_2 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    assert hash(obj_1) == hash(obj_2)
    assert obj_1 == obj_2


def test_check__not_eq_for_diff_func_same_desc_and_version():
    obj_1 = versions.Check("Desc {}", operator.gt, 1, 2, 3)
    obj_2 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    assert hash(obj_1) != hash(obj_2)
    assert obj_1 != obj_2


def test_check__not_eq_for_diff_desc_same_func_and_version():
    obj_1 = versions.Check("Desc1 {}", operator.lt, 1, 2, 3)
    obj_2 = versions.Check("Desc2 {}", operator.lt, 1, 2, 3)
    assert hash(obj_1) != hash(obj_2)
    assert obj_1 != obj_2


def test_check__not_eq_for_same_func_desc_diff_version():
    obj_1 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    obj_2 = versions.Check("Desc {}", operator.lt, 1, 3, 3)
    assert hash(obj_1) != hash(obj_2)
    assert obj_1 != obj_2


###############################################################################
###############################################################################
# EQ class


def test_eq__str__one_value():
    obj = versions.EQ(1)
    assert str(obj) == "v1"


def test_eq__str__two_values():
    obj = versions.EQ(2, 1)
    assert str(obj) == "v2.1"


def test_eq__check_values__equal():
    obj = versions.EQ(2, 3)
    assert obj((2, 3))


def test_eq__check_values__not_equal():
    obj = versions.EQ(2, 3)
    assert not obj((1, 3))
    assert not obj((2, 2))
    assert not obj((1, 4))


def test_eq__check_values__equal_truncated():
    obj = versions.EQ(2, 3)
    assert obj((2, 3, 1))


def test_eq__check_values__equal_too_few_values():
    obj = versions.EQ(2, 3)
    with pytest.raises(ValueError):
        obj((2,))


def test_eq__check_values__not_equal_too_few_values():
    obj = versions.EQ(2, 3)
    with pytest.raises(ValueError):
        obj((1,))


###############################################################################
###############################################################################
# GE class


def test_ge__str__one_value():
    obj = versions.GE(1)
    assert str(obj) == "at least v1"


def test_ge__str__two_values():
    obj = versions.GE(2, 1)
    assert str(obj) == "at least v2.1"


def test_ge__check_values__greater_than_or_equal():
    obj = versions.GE(2, 3)
    assert obj((2, 3))
    assert obj((2, 4))
    assert obj((3, 0))


def test_ge__check_values__not_greater_than_or_equal():
    obj = versions.GE(2, 3)
    assert not obj((1, 3))
    assert not obj((2, 2))


def test_ge__check_values__greater_than_or_equal_truncated():
    obj = versions.GE(2, 3)
    assert obj((2, 3, 1))
    assert obj((2, 4, 2))


def test_ge__check_values__equal_too_few_values():
    obj = versions.GE(2, 3)
    with pytest.raises(ValueError):
        obj((2,))


def test_ge__check_values__not_equal_too_few_values():
    obj = versions.GE(2, 3)
    with pytest.raises(ValueError):
        obj((1,))


###############################################################################
###############################################################################
# LT class


def test_lt__str__one_value():
    obj = versions.LT(1)
    assert str(obj) == "prior to v1"


def test_lt__str__two_values():
    obj = versions.LT(2, 1)
    assert str(obj) == "prior to v2.1"


def test_lt__check_values__less_than():
    obj = versions.LT(2, 3)
    assert obj((2, 2))
    assert obj((1, 9))


def test_lt__check_values__not_less_than():
    obj = versions.LT(2, 3)
    assert not obj((2, 3))
    assert not obj((2, 4))


def test_lt__check_values__less_than_truncated():
    obj = versions.LT(2, 3)
    assert obj((2, 2, 1))
    assert obj((2, 1, 2))


def test_lt__check_values__less_than_too_few_values():
    obj = versions.LT(2, 3)
    with pytest.raises(ValueError):
        obj((1,))


def test_lt__check_values__not_less_than_too_few_values():
    obj = versions.LT(2, 3)
    with pytest.raises(ValueError):
        obj((3,))


###############################################################################
###############################################################################
# Any class


def test_any__str():
    obj = versions.Any()
    assert str(obj) == "any version"


def test_lt__check_values__always_true():
    obj = versions.Any()
    assert obj((1,))
    assert obj((2, 3))
    assert obj((4, 5, 6))
    assert obj((5, 6, 7, 8))


###############################################################################
###############################################################################
# And class


def test_and__init__non_check_value():
    with pytest.raises(ValueError):
        versions.And(versions.LT(2), None)


###############################################################################
###############################################################################
# And class -- str


def test_and__str__single_item():
    obj = versions.And(versions.GE(1))
    assert str(obj) == "at least v1"


def test_and__str__two_items():
    obj_ge = versions.GE(1, 2)
    obj_lt = versions.LT(3, 4)
    obj = versions.And(obj_ge, obj_lt)

    assert str(obj) == "at least v1.2 and prior to v3.4"


def test_and__str__two_items__first_is_operator():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.LT(3, 4)
    obj = versions.And(obj_1, obj_2)

    assert str(obj) == "(at least v1.2 and prior to v2.0) and prior to v3.4"


def test_and__str__two_items__second_is_operator():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.Or(versions.GE(2, 0), versions.LT(3, 4))
    obj = versions.And(obj_1, obj_2)

    assert str(obj) == "at least v1.2 and (at least v2.0 or prior to v3.4)"


###############################################################################
###############################################################################
# And class -- check_version


def test_and__check_version__both_true():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.LT(2, 0)
    obj = versions.And(obj_1, obj_2)
    assert obj((1, 3))


def test_and__check_version__first_true():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.And(versions.GE(2, 3), versions.LT(3, 0))
    obj = versions.And(obj_1, obj_2)
    assert not obj((1, 3))


def test_and__check_version__second_true():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.And(versions.GE(2, 3), versions.LT(3, 0))
    obj = versions.And(obj_1, obj_2)
    assert not obj((2, 3))


def test_and__check_version__neither_true():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.And(versions.GE(2, 3), versions.LT(3, 0))
    obj = versions.And(obj_1, obj_2)
    assert not obj((2, 2))


_TRUNCATED_VERSIONS = (
    (versions.GE(1, 2), versions.LT(2, 0)),
    (versions.GE(1, 2, 2), versions.LT(2, 0)),
    (versions.GE(1, 2), versions.LT(2, 0, 1)),
)


@pytest.mark.parametrize("obj_1, obj_2", _TRUNCATED_VERSIONS)
def test_and__check_version__truncated(obj_1, obj_2):
    obj = versions.And(obj_1, obj_2)
    assert obj((1, 3, 3))


_INSUFFICIENT_VALUES = (
    (versions.GE(1, 2, 2), versions.LT(2, 0)),
    (versions.GE(1, 2), versions.LT(2, 0, 1)),
    (versions.GE(1, 2, 2), versions.LT(2, 0, 1)),
)


@pytest.mark.parametrize("obj_1, obj_2", _INSUFFICIENT_VALUES)
def test_and__check_version__insufficient_number_of_values(obj_1, obj_2):
    obj = versions.And(obj_1, obj_2)
    with pytest.raises(ValueError):
        obj((1, 3))


###############################################################################
###############################################################################
# Or class


def test_or__init__non_check_value():
    with pytest.raises(ValueError):
        versions.Or(versions.LT(2), None)


###############################################################################
###############################################################################
# Or class -- str


def test_or__str__single_item():
    obj = versions.Or(versions.GE(1))
    assert str(obj) == "at least v1"


def test_or__str__two_items():
    obj_ge = versions.GE(1, 2)
    obj_lt = versions.LT(3, 4)
    obj = versions.Or(obj_ge, obj_lt)

    assert str(obj) == "at least v1.2 or prior to v3.4"


def test_or__str__two_items__first_is_operator():
    obj_1 = versions.Or(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.LT(3, 4)
    obj = versions.Or(obj_1, obj_2)

    assert str(obj) == "(at least v1.2 or prior to v2.0) or prior to v3.4"


def test_or__str__two_items__second_is_operator():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.And(versions.GE(2, 0), versions.LT(3, 4))
    obj = versions.Or(obj_1, obj_2)

    assert str(obj) == "at least v1.2 or (at least v2.0 and prior to v3.4)"


###############################################################################
###############################################################################
# Or class -- check_version


def test_or__check_version__both_true():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.LT(2, 0)
    obj = versions.Or(obj_1, obj_2)
    assert obj((1, 3))


def test_or__check_version__first_true():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.And(versions.GE(2, 3), versions.LT(3, 0))
    obj = versions.Or(obj_1, obj_2)
    assert obj((1, 3))


def test_or__check_version__second_true():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.And(versions.GE(2, 3), versions.LT(3, 0))
    obj = versions.Or(obj_1, obj_2)
    assert obj((2, 3))


def test_or__check_version__neither_true():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.And(versions.GE(2, 3), versions.LT(3, 0))
    obj = versions.Or(obj_1, obj_2)
    assert not obj((2, 2))


@pytest.mark.parametrize("obj_1, obj_2", _TRUNCATED_VERSIONS)
def test_or__check_version__truncated(obj_1, obj_2):
    obj = versions.Or(obj_1, obj_2)
    assert obj((1, 3, 3))


_INSUFFICIENT_VALUES_OR = (
    (versions.GE(1, 2, 2), versions.LT(2, 0)),
    (versions.GE(1, 2, 2), versions.LT(2, 0, 1)),
)


@pytest.mark.parametrize("obj_1, obj_2", _INSUFFICIENT_VALUES_OR)
def test_or__check_version__insufficient_number_of_values(obj_1, obj_2):
    obj = versions.Or(obj_1, obj_2)
    with pytest.raises(ValueError):
        obj((1, 3))


def test_or__check_version__insufficient_number_of_values__is_lazy():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.LT(2, 0, 1)
    obj = versions.Or(obj_1, obj_2)
    assert obj((1, 3))


###############################################################################
###############################################################################
# Requirement -- constructor


def test_requirementobj__init__defaults():
    obj = versions.Requirement(
        call=("echo", "foo"), search=r"(\d+)\.(\d+)", checks=versions.Any()
    )

    assert obj.name == "echo"
    assert obj.priority == 0


def test_requirementobj__init__non_defaults():
    obj = versions.Requirement(
        call=("bash", "foo"),
        search=r"(\d+)\.(\d+)",
        checks=versions.Any(),
        name="A name",
        priority=17,
    )

    assert obj.name == "A name"
    assert obj.priority == 17


###############################################################################
###############################################################################
# Requirement -- version


def _echo_version(version, dst="stdout", returncode=0):
    tmpl = "import sys; sys.%s.write(%r); sys.exit(%s);"
    return (sys.executable, "-c", tmpl % (dst, version, returncode))


_PIPES = ("stderr", "stdout")

_VERSION_CALL_RESULTS = (
    (r"v(\d+)", (3,)),
    (r"v(\d+)\.(\d+)", (3, 5)),
    (r"v(\d+)\.(\d+)\.(\d+)", (3, 5, 2)),
)


@pytest.mark.parametrize("pipe", _PIPES)
@pytest.mark.parametrize("regexp, equals", _VERSION_CALL_RESULTS)
def test_requirementobj__version__call(pipe, regexp, equals):
    call = _echo_version("v3.5.2\n", dst=pipe)
    obj = versions.Requirement(call=call, search=regexp, checks=versions.Any())
    assert obj.version() == equals


def test_requirementobj__version__version_str_not_found():
    call = _echo_version("A typical error\n")
    obj = versions.Requirement(
        call=call, search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )

    with pytest.raises(versions.VersionRequirementError):
        obj.version()


def test_requirementobj__version__command_not_found():
    obj = versions.Requirement(
        call=("xyzabcdefoo",), search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )

    with pytest.raises(
        versions.VersionRequirementError, match="No such file or directory"
    ):
        obj.version()


def test_requirementobj__version__command_not_executable():
    obj = versions.Requirement(
        call=("./README.rst",), search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )

    with pytest.raises(versions.VersionRequirementError, match="[Errno 13]"):
        obj.version()


def test_requirementobj__version__return_code_is_ignored():
    obj = versions.Requirement(
        _echo_version("v1.2.3", returncode=1),
        search=r"v(\d+)\.(\d+)",
        checks=versions.Any(),
    )
    assert obj.version() == (1, 2)


def test_requirementobj__version__func_call():
    def _return_version():
        return "This is v5.3!"

    obj = versions.Requirement(
        call=_return_version, search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )
    assert obj.version() == (5, 3)


def test_requirementobj__version__func_call_with_arguments():
    def _return_version(arg1, arg2):
        assert (arg1, arg2) == (2, "foo")
        return "This is v5.3!"

    obj = versions.Requirement(
        call=(_return_version, 2, "foo"), search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )
    assert obj.version() == (5, 3)


@pytest.mark.parametrize(
    "message", ("UnsupportedClassVersionError", "UnsupportedClassVersionError v1.2.3")
)
def test_requirementobj__version__outdated_jre__with_or_without_version_str(message):
    obj = versions.Requirement(
        call=lambda: message, search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )

    with pytest.raises(
        versions.VersionRequirementError, match="upgrade your version of Java"
    ):
        obj.version()


###############################################################################
###############################################################################
# Requirement -- executable


def test_requirementobj__executable__no_cli_args():
    obj = versions.Requirement(
        call=["samtools"], search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )

    assert obj.executable == "samtools"


def test_requirementobj__executable__with_cli_arguments():
    obj = versions.Requirement(
        call=["samtools", "--version"], search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )

    assert obj.executable == "samtools"


def test_requirementobj__executable__function():
    obj = versions.Requirement(
        call=lambda: "v1.1", search=r"v(\d+)\.(\d+)", checks=versions.Any()
    )

    assert obj.executable is None


###############################################################################
###############################################################################
# Requirement -- __call__


class CheckCounted(versions.Check):
    def __init__(self, return_value=True, expected=(1, 1)):
        self.count = 0
        self.return_value = return_value
        versions.Check.__init__(self, "counted {}", operator.eq, *expected)

    def _do_check_version(self, current, reference):
        assert current == reference
        self.count += 1
        return self.return_value


def test_requirementobj__call__check_fails__function():
    obj = versions.Requirement(
        call=lambda: "v1.0.3",
        search=r"(\d)\.(\d)",
        checks=versions.GE(1, 1),
        name="test#1",
    )

    assert not obj.check()


def test_requirementobj__call__check_fails():
    obj = versions.Requirement(
        call=_echo_version("v1.0.2"),
        search=r"(\d)\.(\d)",
        checks=versions.GE(1, 1),
        name="test#1",
    )

    assert not obj.check()


def test_requirementobj__call__check_fails__jre_outdated():
    expected = (
        "The version of the Java Runtime Environment on this\n"
        "system is too old; please check the the requirement\n"
        "for the program and upgrade your version of Java.\n"
    )

    value = "UnsupportedClassVersionError"
    obj = versions.Requirement(
        call=_echo_version(value),
        search=r"(\d)\.(\d)",
        checks=versions.GE(1, 1),
        name="test#1",
    )

    with pytest.raises(versions.VersionRequirementError, match=expected):
        obj.check()


###############################################################################
###############################################################################
# Pickling of checks

_CAN_PICKLE_VALUES = (
    versions.EQ(1, 2, 3),
    versions.GE(1, 2, 3),
    versions.LT(1, 2, 3),
    versions.Any(),
    versions.And(versions.EQ(1, 2, 3)),
    versions.Or(versions.GE(1, 2, 3)),
)


@pytest.mark.parametrize("obj", _CAN_PICKLE_VALUES)
def test_check__can_pickle(obj):
    pickle.dumps(obj)
