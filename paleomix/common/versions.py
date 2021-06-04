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
"""Version checks for apps or libraries required by PALEOMIX pipelines.

The module contains to sets of classes: Requirement and Check. The Requirement class
implements the determation of the current version for a given tool, while the Check
(sub)classes implements various comparison to be carried against the detected version
(less than, greater than or equal, etc.).

For example, to check that the Java version is v1.7 or later:
    obj = Requirement(call=("java", "-version"),
                      search='java version "(\\d+).(\\d+)',
                      checks=GE(1, 7),
                      name="Java Runtime Environment")

    assert obj.check(), "version requirements not met"
"""
import operator
import re
import subprocess

from shlex import quote

from paleomix.common.utilities import TotallyOrdered, safe_coerce_to_tuple, try_cast


class VersionRequirementError(Exception):
    """Raised if version requirements are not met, or if a version could not be
    determined for a requirement check.
    """


class Requirement:
    """Represents a version requirement."""

    def __init__(self, call, search=None, checks=None, name=None, priority=0):
        """See function 'Requrement' for a description of parameters."""
        self._call = safe_coerce_to_tuple(call)
        self.name = str(name or self._call[0])
        self.search = search
        self.checks = Any() if checks is None else checks
        self.priority = int(priority)
        self._cached_version = None

        # Checks will always fail without a search string
        if not (self.search or isinstance(self.checks, Any)):
            raise ValueError(self.checks)

    def version(self, force=False):
        """The version determined for the application / library. If the version
        could not be determined, a VersionRequirementError is raised,
        describing the cause of the problem.
        """
        if force or self._cached_version is None:
            try:
                output = self._run(self._call)
            except OSError as error:
                self._raise_failure(error)

            if isinstance(output, bytes):
                output = output.decode("utf-8", "replace")

            # Raise an exception if the JRE is outdated, even if the
            # version could be determined (likely a false positive match).
            self._check_for_outdated_jre(output)

            if self.search:
                match = re.search(self.search, output)
                if not match:
                    self._raise_failure(output)

                self._cached_version = tuple(
                    0 if value is None else try_cast(value, int)
                    for value in match.groups()
                )
            else:
                self._cached_version = ()

        return self._cached_version

    def version_str(self, force=False):
        return _pprint_version(self.version(force))

    @property
    def executable(self):
        """Returns the executable invoked during version determination; if no
        executable is invoked, None is returned.
        """
        if not callable(self._call[0]):
            return self._call[0]

    def check(self, force=False):
        return self.checks(self.version(force))

    def _check_for_outdated_jre(self, output):
        """Checks for the error raised if the JRE is unable to run a JAR file.
        This happens if the JAR was built with a never version of the JRE, e.g.
        if Picard was built with a v1.7 JRE, and then run with a v1.6 JRE.
        """
        # This exception is raised if the JRE is incompatible with the JAR
        if "UnsupportedClassVersionError" in output:
            self._raise_failure(output)

    def _raise_failure(self, output):
        """Raises a VersionRequirementError when a version check failed; if the
        output indicates that the JRE is outdated (i.e. the output contains
        "UnsupportedClassVersionError") a special message is givenself.
        """
        lines = ["Version could not be determined:"]
        lines.extend(self._describe_call())
        lines.append("")

        if isinstance(output, OSError):
            lines.append("Exception was raised:")
            lines.append("    %s: %s" % (output.__class__.__name__, output))
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
            lines.append("Program may be broken or a version not supported by the")
            lines.append("pipeline; please refer to the PALEOMIX documentation.\n")
            lines.append("Required:       %s" % (self.checks,))
            lines.append("Search string:  %r\n" % (self.search))
            lines.append("%s Command output %s" % ("-" * 22, "-" * 22))
            lines.append(output)

        raise VersionRequirementError("\n".join(lines))

    def _describe_call(self):
        """Returns lines describing the current system call, if any."""
        if not callable(self._call[0]):
            yield "Command  = %s" % (" ".join(map(quote, self._call)),)

    @staticmethod
    def _run(call):
        if callable(call[0]):
            return call[0](*call[1:])

        proc = subprocess.Popen(
            call,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            # Merge STDERR with STDOUT output
            stderr=subprocess.STDOUT,
        )

        return proc.communicate()[0]

    def __eq__(self, other):
        if not isinstance(other, Requirement):
            return NotImplemented

        return (
            self._call == other._call
            and self.name == other.name
            and self.priority == other.priority
            and self.search == other.search
            and self.checks == other.checks
        )

    def __hash__(self):
        return hash((self._call, self.name, self.priority, self.search, self.checks))


class Check(TotallyOrdered):
    """Abstract base-class for version checks.

    Callable with a tuple of version fields (typically integers), and returns
    either True or False, depending on whether or not the specified check
    passed.

    The contructor takes a string describing the check ('description'), a
    function with the signature 'func(version, values)', where version is the
    version determined for a app/library, and where values are the values
    passed to the Check constructor.
    """

    def __init__(self, description, func, *values):
        if not callable(func):
            raise TypeError("func must be callable, not %r" % (func,))

        values = tuple(values)
        self._func = func
        self._values = values
        self._description = description

    def __str__(self):
        return self._description

    def __eq__(self, other):
        if isinstance(other, Check):
            return (
                self._description == other._description
                and self._func == other._func
                and self._values == other._values
            )

        return NotImplemented

    def __hash__(self):
        return hash((self._description, self._func, self._values))

    def __call__(self, current):
        """Takes a tuple of version fields (e.g. (1, 7)) and returns True if
        this version matches the check.
        """
        return self._do_check_version(current, self._values)

    def _do_check_version(self, current, reference):
        """Invokes the actual check; may be overridden in subclasses."""
        return self._func(current, reference)


class CheckVersion(Check):
    """Base class for comparisons involving versions; requires that the version
    checks has at least as many fields as specified for the Check object. If
    the version checked has more fields, these are truncated away.
    """

    def __init__(self, description, func, *version):
        description = description.format(_pprint_version(version))
        Check.__init__(self, description, func, *version)

    def _do_check_version(self, current, reference):
        if len(current) < len(reference):
            raise ValueError(
                "Expects at least %i fields, not %i: %r"
                % (len(reference), len(current), current)
            )

        return Check._do_check_version(self, current[: len(reference)], reference)


class EQ(CheckVersion):
    """Checks that a version is Equal to this version; note that version fields
    are truncated to the number of fields specified for this Check. As a
    consequence, EQ(1, 5) is true for (1, 5), (1, 5, 7), (1, 5, 7, 1), etc. See
    'Check' for more information.
    """

    def __init__(self, *version):
        CheckVersion.__init__(self, "{0}", operator.eq, *version)


class GE(CheckVersion):
    """Checks that a version is Greater-than or Equal to this version; note
    that version fields are truncated to the number of fields specified for
    this Check. See 'Check'.
    """

    def __init__(self, *version):
        CheckVersion.__init__(self, "at least {0}", operator.ge, *version)


class LT(CheckVersion):
    """Checks that a version is Less Than this version; note that version
    fields are truncated to the number of fields specified for this Check.
    See 'Check'.
    """

    def __init__(self, *version):
        CheckVersion.__init__(self, "prior to {0}", operator.lt, *version)


class Any(CheckVersion):
    """Dummy check; is always true."""

    def __init__(self):
        CheckVersion.__init__(self, "any version", _func_any)


class Operator(Check):
    """Base class for logical operations on Checks; and, or, etc."""

    def __init__(self, keyword, func, *checks):
        """Arguments:
        keyword -- Keyword to join description of checks by.
        func -- Function implementing the logical operation; is called as
                func(*checks). See the 'func' argument for Check.__init__.
        checks -- Zero or more Checks.
        """
        descriptions = []
        for check in checks:
            if isinstance(check, Operator):
                descriptions.append("(%s)" % (check,))
            elif isinstance(check, Check):
                descriptions.append("%s" % (check,))
            else:
                raise ValueError("%r is not of type Check" % (check,))

        description = keyword.join(descriptions)
        Check.__init__(self, description, func, *checks)


class And(Operator):
    """Carries out 'and' on a set of checks; always true for no Checks"""

    def __init__(self, *checks):
        Operator.__init__(self, " and ", _func_and, *checks)


class Or(Operator):
    """Carries out 'or' on a set of checks; always false for no Checks"""

    def __init__(self, *checks):
        Operator.__init__(self, " or ", _func_or, *checks)


###############################################################################
###############################################################################
# Check functions; must be available for pickle


def _func_any(_current, _checks):
    """Implementation of Any."""
    return True


def _func_and(current, checks):
    """Implementation of And."""
    return all(check(current) for check in checks)


def _func_or(current, checks):
    """Implementation of Or."""
    return any(check(current) for check in checks)


###############################################################################
###############################################################################
# Utility functions


def _pprint_version(value):
    """Pretty-print version tuple; takes a tuple of field numbers / values,
    and returns it as a string joined by dots with a 'v' prepended.
    """
    if not value:
        return "N/A"

    return "v%s" % (".".join(map(str, value)),)
