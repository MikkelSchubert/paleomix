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
import types
import subprocess

from pypeline.common.utilities import safe_coerce_to_tuple, try_cast


# Cache used to store the output of cmd-line / function calls
_CALL_CACHE = {}
# Cache used to store Requirement object
_REQUIREMENT_CACHE = {}


class VersionRequirementError(RuntimeError):
    pass


def Requirement(call, search, checks, pprint = None, name = None):
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
      pprint -- A function that takes a tuple of version fields and returns
                a string, or a format string that may be used to convert such
                a tuple. If not specified, the version fields will be converted
                to strings and joined by '.'s.
      checks -- A callable that carries out any required version-checks. Is
                called as "checks(value, pprint). Should raise an exception
                (e.g. VersionRequirementError) in the case of requirements
                not being met.

    Implementation detail: To reduce the need for performing calls or system-
    calls multiple times, caches are implemented using the call object as keys.
    Thus the same calls should be passed in a manner which allow equality between
    the same calls to be established.
    """
    call = safe_coerce_to_tuple(call)
    key = (call, search, pprint, checks, name)

    try:
        requirement = _REQUIREMENT_CACHE[key]
    except KeyError:
        requirement = RequirementObj(*key)
        _REQUIREMENT_CACHE[key] = requirement

    return requirement



class RequirementObj:
    def __init__(self, call, search, pprint, checks, name = None):
        self._done = None
        self.name  = name or call[0]
        self._call = call
        self._reqs = checks
        self._rege = re.compile(search)
        self._ppr  = pprint
        self._version = None


    @property
    def version(self):
        if self._version is None:
            output = _do_call(self._call)
            match = self._rege.search(output)
            if not match:
                raise VersionRequirementError("Could not determine version of '%s', searching for %s: %s" \
                                              % (self.name, repr(self._rege.pattern), repr(output)))

            self._version = tuple(try_cast(value, int) for value in match.groups())
        return self._version


    def __call__(self, force = False):
        if force or self._done is None:
            def _pprint(value):
                if not self._ppr:
                    return ".".join(map(str, value))
                elif callable(self._ppr):
                    return self._ppr(value)
                return self._ppr.format(*value)

            self._reqs(self.version, _pprint)
            self._done = True


class _Check:
    def __init__(self, name, *version):
        self._version = tuple(version)
        self._desc    = (str(name), self._version)

    def __hash__(self):
        return hash(self._desc)

    def __cmp__(self, other):
        if isinstance(other, _Check):
            return cmp(self._desc, other._desc)

        return cmp(self.__class__, other.__class__)

    def __str__(self):
        return "%s%s" % self._desc

    def __repr__(self):
        return str(self)


class EQ(_Check):
    def __init__(self, *version):
        _Check.__init__(self, "EQ", *version)

    def __call__(self, value, pprint):
        if value != self._version:
            raise VersionRequirementError("Version must be %s, found %s" \
                                          % (pprint(self._version), pprint(value)))

class GE(_Check):
    def __init__(self, *version):
        _Check.__init__(self, "GE", *version)

    def __call__(self, value, pprint):
        if not value >= self._version:
            raise VersionRequirementError("Version must be at least %s, found %s" \
                                          % (pprint(self._version), pprint(value)))


class LT(_Check):
    def __init__(self, *version):
        _Check.__init__(self, "LE", *version)

    def __call__(self, value, pprint):
        if not value < self._version:
            raise VersionRequirementError("Version must be below %s, found %s" \
                                          % (pprint(self._version), pprint(value)))


class And(_Check):
    def __init__(self, *checks):
        self._checks = checks
        _Check.__init__(self, "Or", *checks)

    def __call__(self, value, pprint):
        for check in self._checks:
            check(value, pprint)


def _run(call):
    proc = subprocess.Popen(call, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdoutdata, stderrdata = proc.communicate()

    return stdoutdata + "\n" + stderrdata


def _do_call(call):
    try:
        return _CALL_CACHE[call]
    except KeyError:
        if callable(call[0]):
            result = call[0](*call[1:])
        else:
            result = _run(call)
        _CALL_CACHE[call] = result
        return result
