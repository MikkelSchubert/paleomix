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
# pylint: disable=W0223
#
"""Version checks for apps or libraries required by PALEOMIX pipelines.

The module contains to sets of classes: RequirementObj and Check. The
RequirementObj class implements the determation of the current version for a
given tool, while the Check (sub)classes implements various comparison
to be carried against the detected version (less than, greater than or equal,
etc.).

To reduce the overhead of detmining versions (which mostly involve invoking
external programs), the RequirementObj caches results. Additionally, to avoid
redundant calls, RequirementObjs are created using the 'Requirement' function
which caches RequirementObjs.

For example, to check that the Java version is v1.7 or later:
    obj = Requirement(call=("java", "-version"),
                      search='java version "(\\d+).(\\d+)',
                      checks=GE(1, 7),
                      name="Java Runtime Environment")
    try:
        obj()
    except VersionRequirementError:
        pass  # requirements not met, or failure to determine version
"""
import re
import operator
import collections

from paleomix.common.utilities import \
    Immutable, \
    TotallyOrdered, \
    safe_coerce_to_tuple, \
    try_cast
from paleomix.common.fileutils import \
    which_executable

import paleomix.common.procs as procs


# Cache used to store the output of cmd-line / function calls
_CALL_CACHE = {}
# Cache used to store Requirement object
_REQUIREMENT_CACHE = {}


class VersionRequirementError(StandardError):
    """Raised if version requirements are not met, or if a version could not be
    determined for a requirement check.
    """


def Requirement(call, search, checks, name=None, priority=0):
    # Ignore function naming scheme
    # pylint: disable=C0103
    """Returns a singleton Requirement object, based on the parameters,
    which may be used to check that version requirements are met for a
    given program/utility/module, etc.

    Parameters:
      call   -- A string, or a tuple containing strings for a system call,
                or a tuple containing a function at the first position, and
                a set of positional parameters. In the case of system calls,
                stdout and stderr are returned as a single string, in the case
                of a function call, the return value is expected to be a str.
      search -- A regular expression (string or re object), used to search
                the output of the "call". Groups are assumed to represent
                version numbers.
      checks -- A callable that implements the interface described in the
                Check class.
      name   -- Descriptive name for the executable/module/etc. If not
                specified, first value in 'call' will be used; if multiple
                otherwise identical checks are made, the last name that
                does not equal the first value of 'call' will be used.
      priority -- Order in which requirements are checked; if multiple
                  otherwise identical checks are made with different priority,
                  the highest priority takes precedence.

    Implementation detail: To reduce the need for performing calls or system-
    calls multiple times, caches are implemented using the call object as keys.
    Thus the same calls should be passed in a manner which allow equality
    between the same calls to be established.
    """
    call = safe_coerce_to_tuple(call)
    key = (call, search, checks)

    try:
        requirement = _REQUIREMENT_CACHE[key]

        # Highest priority takes precedence
        requirement.priority = max(requirement.priority, priority)
        # Last explicitly specified name takes precedence
        requirement.name = name or requirement.name
    except KeyError:
        requirement = RequirementObj(*key, name=name, priority=priority)
        _REQUIREMENT_CACHE[key] = requirement

    return requirement


class RequirementObj(object):
    """Represents a version requirement."""

    def __init__(self, call, search, checks, name=None, priority=0):
        """See function 'Requrement' for a description of parameters.
        """
        self._call = safe_coerce_to_tuple(call)
        self._done = None
        self.name = str(name or self._call[0])
        self.priority = int(priority)
        self.checks = checks
        self._rege = re.compile(search)
        self._version = None

    @property
    def version(self):
        """The version determined for the application / library. If the version
        could not be determined, a VersionRequirementError is raised,
        describing the cause of the problem.
        """
        if self._version is None:
            output = _do_call(self._call)
            # Raise an exception if the JRE is outdated, even if the
            # version could be determined (likely a false positive match).
            self._check_for_outdated_jre(output)

            match = self._rege.search(output)
            if not match:
                self._raise_failure(output)

            self._version = tuple(0 if value is None else try_cast(value, int)
                                  for value in match.groups())

        return self._version

    @property
    def executable(self):
        """Returns the executable invoked during version determination; if no
        executable is invoked, None is returned.
        """
        if not isinstance(self._call[0], collections.Callable):
            return self._call[0]

    def __call__(self, force=False):
        if force or self._done is None:
            if not self.checks(self.version):
                lines = ["Version requirements not met for %s; please refer\n"
                         "to the PALEOMIX documentation for more information."
                         "\n" % (self.name,)]
                lines.extend(self._describe_call())

                version = _pprint_version(self.version)
                lines.append("    Version:       %s" % version)
                lines.append("    Required:      %s" % self.checks)

                raise VersionRequirementError("\n".join(lines))

            self._done = True

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
        lines = ["Version could not be determined for %s:" % (self.name,)]
        lines.append("")
        lines.extend(self._describe_call())
        lines.append("")

        # Raised if the JRE is too old compared to the JAR
        if "UnsupportedClassVersionError" in output:
            lines.extend([
                "The version of the Java Runtime Environment on this",
                "system is too old; please check the the requirement",
                "for the program and upgrade your version of Java.",
                "",
                "See the documentation for more information.",
            ])
        else:
            lines.append("Program may be broken or a version not supported by the")
            lines.append("pipeline; please refer to the PALEOMIX documentation.\n")
            lines.append("    Required:       %s" % (self.checks,))
            lines.append("    Search string:  %r\n" % (self._rege.pattern))
            lines.append("%s Command output %s" % ("-" * 22, "-" * 22))
            lines.append(output)

        raise VersionRequirementError("\n".join(lines))

    def _describe_call(self):
        """Yields string describing the current system call, if any.
        """
        if self.executable:
            exec_path = which_executable(self.executable) or self.executable
            yield "    Executable:    %s" % (exec_path,)

        if not isinstance(self._call[0], collections.Callable):
            yield "    Call:          %s" % (" ".join(self._call),)


class Check(Immutable, TotallyOrdered):
    # Ignore "missing" members; required due to use of Immutable
    # pylint: disable=E1101
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
            raise TypeError('func must be callable, not %r' % (func,))

        values = tuple(values)
        Immutable.__init__(self,
                           _func=func,
                           _values=values,
                           _description=description,
                           _objs=(description, func, values))

    def __str__(self):
        return self._description

    def __lt__(self, other):
        if isinstance(other, Check):
            return self._objs < other._objs  # pylint: disable=W0212
        return NotImplemented  # pragma: no coverage

    def __hash__(self):
        return hash(self._objs)

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
            raise ValueError("Expects at least %i fields, not %i: %r"
                             % (len(reference), len(current), current))

        return Check._do_check_version(self,
                                       current[:len(reference)],
                                       reference)


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

def _run(call):
    """Carries out a system call and returns STDOUT and STDERR as a combined
    string. If an OSError is raied (e.g. due to missing executables), the
    resulting message is returned as a string.
    """
    try:
        proc = procs.open_proc(call,
                               stdout=procs.PIPE,
                               # Merge STDERR with STDOUT output
                               stderr=procs.STDOUT)

        return proc.communicate()[0]
    except (OSError, procs.CalledProcessError), error:
        return str(error)


def _do_call(call):
    """Performs a call; the result is cached, and returned upon subsequent
    calls with the same signature (either a function call or system call).
    """
    try:
        return _CALL_CACHE[call]
    except KeyError:
        if callable(call[0]):
            result = call[0](*call[1:])
        else:
            result = _run(call)
        _CALL_CACHE[call] = result
        return result


def _pprint_version(value):
    """Pretty-print version tuple; takes a tuple of field numbers / values,
    and returns it as a string joined by dots with a 'v' prepended.
    """
    return "v%s.x" % (".".join(map(str, value)),)
