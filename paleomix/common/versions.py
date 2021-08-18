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
import re
import subprocess
import sys
import typing
from shlex import quote
from typing import Iterable, NoReturn, Optional, Tuple, Union

from packaging.specifiers import SpecifierSet


class RequirementError(Exception):
    """Raised if version requirements are not met, or if a version could not be
    determined for a requirement check.
    """


class Requirement:
    """Version checks for executables required by PALEOMIX pipelines.

    A requirement consists of three parts, namely a command that is to be invoked, a
    regexp that extracts an 'X.Y.Z' version string from the output of the command, and a
    standard version specification (e.g. '>=X.Y.Z') to determine if the extracted
    version is acceptable.

    For example, to check that the Java version is v1.7 or later:
        obj = Requirement(call=("java", "-version"),
                        regexp='java version "(\\d+\\.\\d+)',
                        specifiers=">=1.7",
                        name="Java Runtime Environment")

        assert obj.check(), "version requirements not met"
    """

    def __init__(
        self,
        call: Union[str, Iterable[str]],
        regexp: Optional[str] = None,
        specifiers: Optional[str] = None,
        name: Optional[str] = None,
    ):
        self._call = (call,) if isinstance(call, str) else tuple(call)
        self.name = str(name or self._call[0])
        self.regexp = re.compile(regexp) if regexp else None
        self.specifiers = SpecifierSet(specifiers or "")
        self._has_cached_version = False
        self._cached_version = ""

        # Checks will always fail without a regexp string
        if specifiers and not regexp:
            raise RequirementError("specifiers require a regexp str")

    @property
    def call(self) -> Tuple[str, ...]:
        call = self._call
        if call[0] == "%(PYTHON)s":
            return (sys.executable,) + call[1:]

        return call

    def version(self, force: bool = False) -> str:
        """The version determined for the application / library. If the version
        could not be determined, a RequirementError is raised.
        """
        if force or not self._has_cached_version:
            self._cached_version = self._determine_version()
            self._has_cached_version = True

        if isinstance(self._cached_version, Exception):
            raise self._cached_version

        return self._cached_version

    def version_str(self, force: bool = False) -> str:
        version = self.version(force)
        if not version:
            return "N/A"

        return f"v{version}"

    @property
    def executable(self) -> str:
        return self.call[0]

    def check(self, force: bool = False) -> bool:
        version = self.version(force)
        if not self.specifiers:
            return True

        return version in self.specifiers

    def _determine_version(self) -> Union[str, Exception]:
        try:
            output = subprocess.run(
                self.call,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                # Merge STDERR with STDOUT output
                stderr=subprocess.STDOUT,
                encoding="utf-8",
                errors="replace",
            ).stdout
        except OSError as error:
            return self._raise_failure(error)

        # Raise an exception if the JRE is outdated
        if "UnsupportedClassVersionError" in output:
            return self._raise_failure(output)

        return self._parse_version_string(output)

    def _parse_version_string(self, output: str) -> Union[str, Exception]:
        if self.regexp:
            match = self.regexp.search(output)
            if not match:
                return self._raise_failure(output)

            fields = list(match.groups())
            # Support for optional minor versions using optional fields
            while fields and fields[-1] is None:
                fields.pop(-1)

            return ".".join(field.strip(".") for field in fields)

        return ""

    def _raise_failure(self, output: Union[str, Exception]) -> RequirementError:
        """Raises a RequirementError when a version check failed; if the
        output indicates that the JRE is outdated (i.e. the output contains
        "UnsupportedClassVersionError") a special message is given.
        """
        lines = [
            "Version could not be determined for:",
            "Command  = {}".format(" ".join(map(quote, self.call))),
            "",
        ]

        origin = None
        if isinstance(output, Exception):
            origin = output
            lines.append("Exception was raised: %r" % (output,))
        elif "UnsupportedClassVersionError" in output:
            # Raised if the JRE is too old compared to the JAR
            lines.extend(
                [
                    "The version of the Java Runtime Environment on this",
                    "system is too old; please check the the requirement",
                    "for the program and upgrade your version of Java.",
                    "",
                    "See the documentation for more information.",
                ]
            )
        else:
            lines.extend(
                [
                    "Program may be broken or a version not supported by the",
                    "pipeline; please refer to the PALEOMIX documentation.",
                    "",
                    "Requirements:   %s" % (self.specifiers,),
                    "Search string:  %r" % (self.regexp),
                    "",
                    "%s Command output %s" % ("-" * 22, "-" * 22),
                    output,
                ]
            )

        raise RequirementError("\n".join(lines)) from origin

    def __eq__(self, other: typing.Any) -> bool:
        if not isinstance(other, Requirement):
            return NotImplemented

        return self._to_tuple() == other._to_tuple()

    def __hash__(self):
        return hash(self._to_tuple())

    def _to_tuple(self):
        return (self._call, self.name, self.regexp, self.specifiers)
