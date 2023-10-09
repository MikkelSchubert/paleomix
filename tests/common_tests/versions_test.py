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

import re
import sys
from typing import Tuple

import pytest

from paleomix.common.versions import Requirement, RequirementError

###############################################################################
###############################################################################
# Requirement -- constructor


def test_requirement__init__str():
    obj = Requirement("true")

    assert obj.name == "true"
    assert obj.call == ("true",)


def test_requirement__init__list():
    obj = Requirement(["also", "true"])

    assert obj.name == "also"
    assert obj.call == ("also", "true")


def test_requirement__init__self_call__str():
    obj = Requirement("%(PYTHON)s")

    assert obj.name == "%(PYTHON)s"
    assert obj.call == (sys.executable,)


def test_requirement__init__self_call__list():
    obj = Requirement(["%(PYTHON)s", "/path/to/blah.py"])

    assert obj.name == "%(PYTHON)s"
    assert obj.call == (sys.executable, "/path/to/blah.py")


def test_requirement__init__defaults():
    obj = Requirement(call=("echo", "foo"), regexp=r"(\d+)\.(\d+)")

    assert obj.call == ("echo", "foo")
    assert obj.regexp == re.compile(r"(\d+)\.(\d+)")
    assert obj.name == "echo"


def test_requirement__init__non_defaults():
    obj = Requirement(
        call=("bash", "foo"),
        regexp=r"(\d+)\.(\d+)",
        name="A name",
    )

    assert obj.call == ("bash", "foo")
    assert obj.regexp == re.compile(r"(\d+)\.(\d+)")
    assert obj.name == "A name"


def test_requirement__checks_require_search():
    with pytest.raises(RequirementError, match="specifiers require a regexp str"):
        Requirement("true", specifiers=">1.2.3")


###############################################################################
###############################################################################
# Requirement -- version


def _echo_version(version: str, dst: str = "stdout", returncode: int = 0):
    tmpl = "import sys; sys.%s.write(%r); sys.exit(%s);"
    return (sys.executable, "-c", tmpl % (dst, version, returncode))


_PIPES = ("stderr", "stdout")

_VERSION_CALL_RESULTS = (
    (r"v(\d+)", "3"),
    (r"v(\d+\.\d+)", "3.5"),
    (r"v(\d+\.\d+)\.(\d+)", "3.5.2"),
)


@pytest.mark.parametrize("pipe", _PIPES)
@pytest.mark.parametrize("regexp, equals", _VERSION_CALL_RESULTS)
def test_requirement__version__call(pipe: str, regexp: str, equals: Tuple[int, ...]):
    call = _echo_version("v3.5.2\n", dst=pipe)
    obj = Requirement(call=call, regexp=regexp)
    assert obj.version() == equals


def test_requirement__version__version_str_not_found():
    call = _echo_version("A typical error\n")
    obj = Requirement(call=call, regexp=r"v(\d+\.\d+)")

    with pytest.raises(RequirementError):
        obj.version()


def test_requirement__version__command_not_found():
    obj = Requirement(call=("xyzabcdefoo",), regexp=r"v(\d+\.\d+)")

    with pytest.raises(RequirementError, match="No such file or directory"):
        obj.version()


def test_requirement__version__command_not_executable():
    obj = Requirement(call=("./README.rst",), regexp=r"v(\d+\.\d+)")

    with pytest.raises(RequirementError, match="[Errno 13]"):
        obj.version()


def test_requirement__version__return_code_is_ignored():
    obj = Requirement(_echo_version("v1.2.3", returncode=1), regexp=r"v(\d+\.\d+)")
    assert obj.version() == "1.2"
    assert obj.version_str() == "v1.2"


def test_requirement__version__call_is_cached():
    call = ("echo", "v1.2.3")
    obj = Requirement(call, regexp=r"v(\d+\.\d+)")

    assert obj.call == ("echo", "v1.2.3")
    assert obj.version() == "1.2"

    obj._call = ("echo", "v3.2.1")
    assert obj.call == ("echo", "v3.2.1")
    assert obj.version() == "1.2"


def test_requirement__version__found_optional_fields_are_included():
    call = ("echo", "v1.2.5")
    obj = Requirement(call, regexp=r"v(\d+\.\d+)(\.\d+)?")

    assert obj.version() == "1.2.5"


def test_requirement__version__missing_optional_fields_are_trimmed():
    call = ("echo", "v1.2")
    obj = Requirement(call, regexp=r"v(\d+\.\d+)(\.\d+)?")

    assert obj.version() == "1.2"


def test_requirement__version__calls_without_checks_still_invoked_1():
    call = ("echo", "v1.2")
    obj = Requirement(call)

    assert not obj.version()
    assert obj.version_str() == "N/A"


def test_requirement__version__calls_without_checks_still_invoked_2(tmp_path):
    call = (str(tmp_path / "echo"), "v1.2")
    obj = Requirement(call)

    with pytest.raises(RequirementError, match="Exception was raised"):
        obj.version()


###############################################################################
###############################################################################
# Requirement -- executable


def test_requirement__executable__no_cli_args():
    obj = Requirement(call=["samtools"], regexp=r"v(\d+\.\d+)")

    assert obj.executable == "samtools"


def test_requirement__executable__with_cli_arguments():
    obj = Requirement(call=["samtools", "--version"], regexp=r"v(\d+\.\d+)")

    assert obj.executable == "samtools"


###############################################################################
###############################################################################
# Requirement -- check


def test_requirement__call__check_succeeds():
    obj = Requirement(
        call=_echo_version("v1.0.2"),
        regexp=r"(\d\.\d)",
        specifiers=">=1.0",
        name="test#1",
    )

    assert obj.check()


def test_requirement__call__check_fails():
    obj = Requirement(
        call=_echo_version("v1.0.2"),
        regexp=r"(\d\.\d)",
        specifiers=">=1.1",
        name="test#1",
    )

    assert not obj.check()


def test_requirement__call__check_succeeds_if_no_checks():
    obj = Requirement(call=_echo_version("v1.0.2"))

    assert obj.check()


def test_requirement__call__check_fails__jre_outdated():
    expected = (
        "The version of the Java Runtime Environment on this\n"
        "system is too old; please check the the requirement\n"
        "for the program and upgrade your version of Java.\n"
    )

    value = "UnsupportedClassVersionError"
    obj = Requirement(
        call=_echo_version(value),
        regexp=r"(\d)\.(\d)",
        specifiers=">=1.1",
        name="test#1",
    )

    with pytest.raises(RequirementError, match=expected):
        obj.check()
