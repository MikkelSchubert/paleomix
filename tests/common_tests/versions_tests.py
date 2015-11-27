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
# Disable warnings on strange function names
# pylint: disable=C0103
import pickle
import operator

from nose.tools import \
    assert_is, \
    assert_is_not, \
    assert_in, \
    assert_equal, \
    assert_not_equal, \
    assert_raises

import paleomix.common.versions as versions


###############################################################################
###############################################################################
# Check class

def test_check__func_must_be_callable():
    assert_raises(TypeError, versions.Check, "FooBar", 3, 7, 5)


def test_check_str():
    obj = versions.Check("FooBar", operator.lt, 3, 7, 5)
    assert_equal(str(obj), "FooBar")


###############################################################################
###############################################################################
## Check class -- hash and comparisons

def test_check__eq_same_func_desc_and_version():
    obj_1 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    obj_2 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    assert_equal(hash(obj_1), hash(obj_2))
    assert_equal(obj_1, obj_2)


def test_check__not_eq_for_diff_func_same_desc_and_version():
    obj_1 = versions.Check("Desc {}", operator.gt, 1, 2, 3)
    obj_2 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    assert_not_equal(hash(obj_1), hash(obj_2))
    assert_not_equal(obj_1, obj_2)


def test_check__not_eq_for_diff_desc_same_func_and_version():
    obj_1 = versions.Check("Desc1 {}", operator.lt, 1, 2, 3)
    obj_2 = versions.Check("Desc2 {}", operator.lt, 1, 2, 3)
    assert_not_equal(hash(obj_1), hash(obj_2))
    assert_not_equal(obj_1, obj_2)


def test_check__not_eq_for_same_func_desc_diff_version():
    obj_1 = versions.Check("Desc {}", operator.lt, 1, 2, 3)
    obj_2 = versions.Check("Desc {}", operator.lt, 1, 3, 3)
    assert_not_equal(hash(obj_1), hash(obj_2))
    assert_not_equal(obj_1, obj_2)


###############################################################################
###############################################################################
## EQ class

def test_eq__str__one_value():
    obj = versions.EQ(1)
    assert_equal(str(obj), "v1.x")


def test_eq__str__two_values():
    obj = versions.EQ(2, 1)
    assert_equal(str(obj), "v2.1.x")


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
    assert_raises(ValueError, obj, (2,))


def test_eq__check_values__not_equal_too_few_values():
    obj = versions.EQ(2, 3)
    assert_raises(ValueError, obj, (1,))


###############################################################################
###############################################################################
## GE class

def test_ge__str__one_value():
    obj = versions.GE(1)
    assert_equal(str(obj), "at least v1.x")


def test_ge__str__two_values():
    obj = versions.GE(2, 1)
    assert_equal(str(obj), "at least v2.1.x")


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
    assert_raises(ValueError, obj, (2,))


def test_ge__check_values__not_equal_too_few_values():
    obj = versions.GE(2, 3)
    assert_raises(ValueError, obj, (1,))


###############################################################################
###############################################################################
## LT class

def test_lt__str__one_value():
    obj = versions.LT(1)
    assert_equal(str(obj), "prior to v1.x")


def test_lt__str__two_values():
    obj = versions.LT(2, 1)
    assert_equal(str(obj), "prior to v2.1.x")


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
    assert_raises(ValueError, obj, (1,))


def test_lt__check_values__not_less_than_too_few_values():
    obj = versions.LT(2, 3)
    assert_raises(ValueError, obj, (3,))


###############################################################################
###############################################################################
## Any class

def test_any__str():
    obj = versions.Any()
    assert_equal(str(obj), "any version")


def test_lt__check_values__always_true():
    obj = versions.Any()
    assert obj((1,))
    assert obj((2, 3))
    assert obj((4, 5, 6))
    assert obj((5, 6, 7, 8))


###############################################################################
###############################################################################
## And class

def test_and__init__non_check_value():
    assert_raises(ValueError, versions.And, versions.LT(2), None)


###############################################################################
###############################################################################
## And class -- str

def test_and__str__single_item():
    obj = versions.And(versions.GE(1))
    assert_equal(str(obj), "at least v1.x")


def test_and__str__two_items():
    obj_ge = versions.GE(1, 2)
    obj_lt = versions.LT(3, 4)
    obj = versions.And(obj_ge, obj_lt)

    assert_equal(str(obj), "at least v1.2.x and prior to v3.4.x")


def test_and__str__two_items__first_is_operator():
    obj_1 = versions.And(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.LT(3, 4)
    obj = versions.And(obj_1, obj_2)

    assert_equal(str(obj),
                 "(at least v1.2.x and prior to v2.0.x) and prior to v3.4.x")


def test_and__str__two_items__second_is_operator():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.Or(versions.GE(2, 0), versions.LT(3, 4))
    obj = versions.And(obj_1, obj_2)

    assert_equal(str(obj),
                 "at least v1.2.x and (at least v2.0.x or prior to v3.4.x)")


###############################################################################
###############################################################################
## And class -- check_version

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


def test_and__check_version__truncated():
    def _do_and_check_truncated(obj_1, obj_2):
        obj = versions.And(obj_1, obj_2)
        assert obj((1, 3, 3))

    yield _do_and_check_truncated, versions.GE(1, 2), versions.LT(2, 0)
    yield _do_and_check_truncated, versions.GE(1, 2, 2), versions.LT(2, 0)
    yield _do_and_check_truncated, versions.GE(1, 2), versions.LT(2, 0, 1)


def test_and__check_version__insufficient_number_of_values():
    def _do_and_check_num_values(obj_1, obj_2):
        obj = versions.And(obj_1, obj_2)
        assert_raises(ValueError, obj, (1, 3))

    yield _do_and_check_num_values, versions.GE(1, 2, 2), versions.LT(2, 0)
    yield _do_and_check_num_values, versions.GE(1, 2), versions.LT(2, 0, 1)
    yield _do_and_check_num_values, versions.GE(1, 2, 2), versions.LT(2, 0, 1)


###############################################################################
###############################################################################
## Or class

def test_or__init__non_check_value():
    assert_raises(ValueError, versions.Or, versions.LT(2), None)


###############################################################################
###############################################################################
## Or class -- str

def test_or__str__single_item():
    obj = versions.Or(versions.GE(1))
    assert_equal(str(obj), "at least v1.x")


def test_or__str__two_items():
    obj_ge = versions.GE(1, 2)
    obj_lt = versions.LT(3, 4)
    obj = versions.Or(obj_ge, obj_lt)

    assert_equal(str(obj), "at least v1.2.x or prior to v3.4.x")


def test_or__str__two_items__first_is_operator():
    obj_1 = versions.Or(versions.GE(1, 2), versions.LT(2, 0))
    obj_2 = versions.LT(3, 4)
    obj = versions.Or(obj_1, obj_2)

    assert_equal(str(obj),
                 "(at least v1.2.x or prior to v2.0.x) or prior to v3.4.x")


def test_or__str__two_items__second_is_operator():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.And(versions.GE(2, 0), versions.LT(3, 4))
    obj = versions.Or(obj_1, obj_2)

    assert_equal(str(obj),
                 "at least v1.2.x or (at least v2.0.x and prior to v3.4.x)")


###############################################################################
###############################################################################
## Or class -- check_version

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


def test_or__check_version__truncated():
    def _do_or_check_truncated(obj_1, obj_2):
        obj = versions.Or(obj_1, obj_2)
        assert obj((1, 3, 3))

    yield _do_or_check_truncated, versions.GE(1, 2), versions.LT(2, 0)
    yield _do_or_check_truncated, versions.GE(1, 2, 2), versions.LT(2, 0)
    yield _do_or_check_truncated, versions.GE(1, 2), versions.LT(2, 0, 1)


def test_or__check_version__insufficient_number_of_values():
    def _do_or_check_num_values(obj_1, obj_2):
        obj = versions.Or(obj_1, obj_2)
        assert_raises(ValueError, obj, (1, 3))

    yield _do_or_check_num_values, versions.GE(1, 2, 2), versions.LT(2, 0)
    yield _do_or_check_num_values, versions.GE(1, 2, 2), versions.LT(2, 0, 1)


def test_or__check_version__insufficient_number_of_values__is_lazy():
    obj_1 = versions.GE(1, 2)
    obj_2 = versions.LT(2, 0, 1)
    obj = versions.Or(obj_1, obj_2)
    assert obj((1, 3))


###############################################################################
###############################################################################
## RequirementObj -- constructor

def test_requirementobj__init__defaults():
    obj = versions.RequirementObj(call=("echo", "foo"),
                                  search=r"(\d+)\.(\d+)",
                                  checks=versions.Any())

    assert_equal(obj.name, "echo")
    assert_equal(obj.priority, 0)


def test_requirementobj__init__non_defaults():
    obj = versions.RequirementObj(call=("bash", "foo"),
                                  search=r"(\d+)\.(\d+)",
                                  checks=versions.Any(),
                                  name="A name",
                                  priority=17)

    assert_equal(obj.name, "A name")
    assert_equal(obj.priority, 17)


###############################################################################
###############################################################################
## RequirementObj -- version

def _echo_version(version, to="stdout", returncode=0):
    tmpl = "import sys; sys.%s.write(%r); sys.exit(%s);"
    return ("/usr/bin/python", "-c", tmpl % (to, version, returncode))
_PIPES = ("stderr", "stdout")


def test_requirementobj__version__call():
    def _do_test_version__single_digit(pipe, regexp, equals):
        call = _echo_version("v3.5.2\n", to=pipe)
        obj = versions.RequirementObj(call=call,
                                      search=regexp,
                                      checks=versions.Any())
        assert_equal(obj.version, equals)

    for pipe in _PIPES:
        yield _do_test_version__single_digit, pipe, r"v(\d+)", (3,)
        yield _do_test_version__single_digit, pipe, r"v(\d+)\.(\d+)", (3, 5)
        yield _do_test_version__single_digit, pipe, r"v(\d+)\.(\d+)\.(\d+)", \
            (3, 5, 2)


def test_requirementobj__version__version_str_not_found():
    call = _echo_version("A typical error\n")
    obj = versions.RequirementObj(call=call,
                                  search=r"v(\d+)\.(\d+)",
                                  checks=versions.Any())

    assert_raises(versions.VersionRequirementError, getattr, obj, "version")


def test_requirementobj__version__command_not_found():
    obj = versions.RequirementObj(call=("xyzabcdefoo",),
                                  search=r"v(\d+)\.(\d+)",
                                  checks=versions.Any())

    try:
        obj.version  # pylint: disable=
        assert False  # pragma: no coverage
    except versions.VersionRequirementError, error:
        # Should include OSError message
        assert_in("No such file or directory", str(error))


def test_requirementobj__version__return_code_is_ignored():
    obj = versions.RequirementObj(_echo_version("v1.2.3", returncode=1),
                                  search=r"v(\d+)\.(\d+)",
                                  checks=versions.Any())
    assert_equal(obj.version, (1, 2))


def test_requirementobj__version__func_call():
    def _return_version():
        return "This is v5.3!"

    obj = versions.RequirementObj(call=_return_version,
                                  search=r"v(\d+)\.(\d+)",
                                  checks=versions.Any())
    assert_equal(obj.version, (5, 3))


def test_requirementobj__version__func_call_with_arguments():
    def _return_version(arg1, arg2):
        assert_equal((arg1, arg2), (2, "foo"))
        return "This is v5.3!"

    obj = versions.RequirementObj(call=(_return_version, 2, "foo"),
                                  search=r"v(\d+)\.(\d+)",
                                  checks=versions.Any())
    assert_equal(obj.version, (5, 3))


def test_requirementobj__version__outdated_jre__with_or_without_version_str():
    error_msg = "upgrade your version of Java"

    def _do_test_outdated_jre(message):
        obj = versions.RequirementObj(call=lambda: message,
                                      search=r"v(\d+)\.(\d+)",
                                      checks=versions.Any())

        try:
            obj.version
            assert False  # pragma: no coverage
        except versions.VersionRequirementError, error:
            assert_in(error_msg, str(error))

    messages = [
        "UnsupportedClassVersionError",
        "UnsupportedClassVersionError v1.2.3"]

    for message in messages:
        yield _do_test_outdated_jre, message


###############################################################################
###############################################################################
## RequirementObj -- __call__

class CheckCounted(versions.Check):
    def __init__(self, return_value=True, expected=(1, 1)):
        versions.Check.__init__(self, "counted {}", operator.eq, *expected)
        object.__setattr__(self, "count", 0)
        object.__setattr__(self, "return_value", return_value)

    def _do_check_version(self, values, current):
        assert_equal(values, current)
        object.__setattr__(self, "count", self.count + 1)
        return self.return_value


def test_requirementobj__call__result_is_cached():
    counter = CheckCounted()
    obj = versions.RequirementObj(call=lambda: "v1.1.3",
                                  search=r"(\d)\.(\d)",
                                  checks=counter)

    obj()
    obj()

    assert_equal(counter.count, 1)


def test_requirementobj__call__result_is_cached_unless_forced():
    counter = CheckCounted()
    obj = versions.RequirementObj(call=lambda: "v1.1.3",
                                  search=r"(\d)\.(\d)",
                                  checks=counter)

    obj()
    obj(force=True)

    assert_equal(counter.count, 2)


def test_requirementobj__call__check_fails__function():
    expected = \
        "Version requirements not met for test#1; please refer\n" \
        "to the PALEOMIX documentation for more information.\n" \
        "\n" \
        "    Version:       v1.0.x\n" \
        "    Required:      at least v1.1.x"

    obj = versions.RequirementObj(call=lambda: "v1.0.3",
                                  search=r"(\d)\.(\d)",
                                  checks=versions.GE(1, 1),
                                  name="test#1")
    try:
        obj()
        assert False  # pragma: no coverage
    except versions.VersionRequirementError, error:
        assert_equal(str(error), expected)


def test_requirementobj__call__check_fails():
    expected = \
        "Version requirements not met for test#1; please refer\n" \
        "to the PALEOMIX documentation for more information.\n" \
        "\n" \
        "    Executable:    /usr/bin/python\n" \
        "    Call:          /usr/bin/python -c import sys; " \
        "sys.stdout.write('v1.0.2'); sys.exit(0);\n" \
        "    Version:       v1.0.x\n" \
        "    Required:      at least v1.1.x"

    obj = versions.RequirementObj(call=_echo_version("v1.0.2"),
                                  search=r"(\d)\.(\d)",
                                  checks=versions.GE(1, 1),
                                  name="test#1")
    try:
        obj()
        assert False  # pragma: no coverage
    except versions.VersionRequirementError, error:
        assert_equal(str(error), expected)


def test_requirementobj__call__check_fails__jre_outdated():
    expected = \
        "Version could not be determined for test#1:\n" \
        "\n" \
        "    Executable:    /usr/bin/python\n" \
        "    Call:          /usr/bin/python -c import sys; " \
        "sys.stdout.write('UnsupportedClassVersionError'); sys.exit(0);\n" \
        "\n" \
        "The version of the Java Runtime Environment on this\n" \
        "system is too old; please check the the requirement\n" \
        "for the program and upgrade your version of Java.\n" \
        "\n" \
        "See the documentation for more information."

    value = "UnsupportedClassVersionError"
    obj = versions.RequirementObj(call=_echo_version(value),
                                  search=r"(\d)\.(\d)",
                                  checks=versions.GE(1, 1),
                                  name="test#1")
    try:
        obj()
        assert False  # pragma: no coverage
    except versions.VersionRequirementError, error:
        assert_equal(str(error), expected)


###############################################################################
###############################################################################
## Pickling of checks

def test_check__can_pickle():
    def _do_test_can_pickle(obj):
        pickle.dumps(obj)

    yield _do_test_can_pickle, versions.EQ(1, 2, 3)
    yield _do_test_can_pickle, versions.GE(1, 2, 3)
    yield _do_test_can_pickle, versions.LT(1, 2, 3)
    yield _do_test_can_pickle, versions.Any()
    yield _do_test_can_pickle, versions.And(versions.EQ(1, 2, 3))
    yield _do_test_can_pickle, versions.Or(versions.GE(1, 2, 3))


###############################################################################
###############################################################################
## Requirement

def test_requirement__obj_is_cached_for_same_values():
    obj1 = versions.Requirement("echo", "", versions.LT(1))
    obj2 = versions.Requirement("echo", "", versions.LT(1))
    assert_is(obj1, obj2)


def test_requirement__new_obj_if_call_differ():
    obj1 = versions.Requirement("echo", "", versions.LT(1))
    obj2 = versions.Requirement("true", "", versions.LT(1))
    assert_is_not(obj1, obj2)


def test_requirement__new_obj_if_search_differ():
    obj1 = versions.Requirement("echo", r"(\d+)", versions.LT(1))
    obj2 = versions.Requirement("echo", "", versions.LT(1))
    assert_is_not(obj1, obj2)


def test_requirement__new_obj_if_checks_differ():
    obj1 = versions.Requirement("echo", "", versions.GE(1))
    obj2 = versions.Requirement("echo", "", versions.LT(1))
    assert_is_not(obj1, obj2)


def test_requirement__same_obj_if_name_differ():
    obj1 = versions.Requirement("echo", "", versions.GE(1))
    assert_equal(obj1.name, "echo")
    obj2 = versions.Requirement("echo", "", versions.GE(1), name="foo")
    assert_equal(obj2.name, "foo")
    assert_is(obj1, obj2)

    obj3 = versions.Requirement("echo", "", versions.GE(1), name="bar")
    assert_equal(obj3.name, "bar")
    assert_is(obj2, obj3)

    obj4 = versions.Requirement("echo", "", versions.GE(1))
    assert_equal(obj3.name, "bar")
    assert_is(obj3, obj4)


def test_requirement_highest_priority_overrides():
    obj1 = versions.Requirement("echo", "", versions.LT(1), priority=0)
    assert_equal(obj1.priority, 0)
    obj2 = versions.Requirement("echo", "", versions.LT(1), priority=5)
    assert_is(obj1, obj2)
    assert_equal(obj2.priority, 5)


def test_requirement_highest_priority_retained():
    obj1 = versions.Requirement("echo", "", versions.LT(1), priority=5)
    assert_equal(obj1.priority, 5)
    obj2 = versions.Requirement("echo", "", versions.LT(1), priority=0)
    assert_is(obj1, obj2)
    assert_equal(obj2.priority, 5)
