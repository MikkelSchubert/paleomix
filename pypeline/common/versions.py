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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
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
import collections

from pypeline.common.utilities import \
    safe_coerce_to_tuple, \
    try_cast
from pypeline.common.fileutils import \
    which_executable


# Cache used to store the output of cmd-line / function calls
_CALL_CACHE = {}
# Cache used to store Requirement object
_REQUIREMENT_CACHE = {}


class VersionRequirementError(RuntimeError):
    pass


def Requirement(call, search, checks, name = None):
    """Returns a singleton Requirement object, based on the parameters,
    which may be used to check that version requirements are met for a
    given program/utility/module, etc.

    Parameters:
      name   -- Descriptive name for the executable/module/etc. If not
                specified, first value in call will be used.
      call   -- A string, or a tuple containing strings for a system call,
                or a tuple containing a function at the first position, and
                a set of positional parameters. In the case of system calls,
                stdout and stderr are returned as a single string.
      search -- A regular expression (string or re object), used to search
                the output of the "call". Groups are assumed to represent
                version numbers.
      checks -- A callable that implements the interface described in the
                Check class.

    Implementation detail: To reduce the need for performing calls or system-
    calls multiple times, caches are implemented using the call object as keys.
    Thus the same calls should be passed in a manner which allow equality between
    the same calls to be established.
    """
    call = safe_coerce_to_tuple(call)
    key = (call, search, checks, name)

    try:
        requirement = _REQUIREMENT_CACHE[key]
    except KeyError:
        requirement = RequirementObj(*key)
        _REQUIREMENT_CACHE[key] = requirement

    return requirement


class RequirementObj:
    def __init__(self, call, search, checks, name=None):
        self._done = None
        self.name = name or call[0]
        self._call = call
        self._reqs = checks
        self._rege = re.compile(search)
        self._version = None

    @property
    def version(self):
        if self._version is None:
            output = _do_call(self._call)
            match = self._rege.search(output)
            if not match:
                self._raise_failure(output)

            self._version = tuple(try_cast(value, int)
                                  for value in match.groups())

        return self._version

    def __call__(self, force=False):
        if force or self._done is None:
            self._reqs(self.name, self.version, self.executable)
            self._done = True

    @property
    def executable(self):
        if not isinstance(self._call[0], collections.Callable):
            return self._call[0]

    def _raise_failure(self, output):
        lines = ["Version could not be determined for %r:"
                 % (self.name,)]

        if self.executable:
            exec_path = which_executable(self.executable)
            lines.append("    Executable:     %s" % (exec_path,))

        # Raised if the JRE is too old compared to the JAR
        if "java.lang.UnsupportedClassVersionError" in output:
            lines.extend([
                "",
                "The version of the Java Runtime Environment on this",
                "system is too old; please check the the requirement",
                "for the program and upgrade your version of Java.",
                "",
                "See the WIKI for more information:",
                "https://github.com/MikkelSchubert/paleomix/wiki/"
                "Troubleshooting#JRE_outdated"
            ])
        else:
            lines.append("    Required: %s" % (self._reqs.description,))
            lines.append("    Regexp:   %r" % (self._rege.pattern))
            lines.append("    Command output:\n %r" % (output))

        raise VersionRequirementError("\n".join(lines))


class Check:
    def __init__(self, name, description, *version):
        self._version    = tuple(version)
        self._objs       = (str(name), self._version)
        self.description = description.format(_pprint(self._version))


    def __hash__(self):
        return hash(self._objs)


    def __cmp__(self, other):
        if isinstance(other, Check):
            return cmp(self._objs, other._objs) # pylint: disable=W0212
        return cmp(self.__class__, other.__class__)


    def __call__(self, name, value, executable = None):
        if not self.check_version(value):
            version = _pprint(value)

            lines = ["Version requirement not met for %r:" % (name,)]
            if executable:
                exec_path = which_executable(executable)
                lines.append("    Path:     %s" % (exec_path,))

            lines.append("    Version:  %s" % (version,))
            lines.append("    Required: %s" % (self.description,))

            raise VersionRequirementError("\n".join(lines))

    def check_version(self, value):
        raise NotImplementedError("'check_version' not implemented")


class EQ(Check):
    def __init__(self, *version):
        Check.__init__(self, "EQ", "equals {0}", *version)


    def check_version(self, value):
        return (value == self._version)


class GE(Check):
    def __init__(self, *version):
        Check.__init__(self, "GE", "at least {0}", *version)


    def check_version(self, value):
        return value >= self._version


class LT(Check):
    def __init__(self, *version):
        Check.__init__(self, "LE", "prior to {0}", *version)


    def check_version(self, value):
        return value < self._version


class And(Check):
    def __init__(self, *checks):
        self._checks = checks
        description = " and ".join("(%s)" % (check.description,) for check in checks)
        Check.__init__(self, "And", description, *checks)


    def check_version(self, value):
        return all(check.check_version(value) for check in self._checks)


class Or(Check):
    def __init__(self, *checks):
        self._checks = checks
        description = " or ".join("(%s)" % (check.description,) for check in checks)
        Check.__init__(self, "Or", description, *checks)


    def check_version(self, value):
        return any(check.check_version(value) for check in self._checks)


def _run(call):
    try:
        proc = subprocess.Popen(call,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.STDOUT)

        return proc.communicate()[0]
    except (OSError, subprocess.CalledProcessError), error:
        return str(error)


def _do_call(call):
    try:
        return _CALL_CACHE[call]
    except KeyError:
        if isinstance(call[0], collections.Callable):
            result = call[0](*call[1:])
        else:
            result = _run(call)
        _CALL_CACHE[call] = result
        return result


def _pprint(value):
    return "v" + ".".join(map(str, value))
