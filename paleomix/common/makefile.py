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
"""Generalized methods for parsing/validating "makefiles" in YAML format.

The following example will use a imagined makefile for 'vcf_filter', which
takes a set of input files, each of which is assigned an output file, and
each of which may have a set of filters (in addition to a set of defaults):

Example makefile in YAML format:
-------------------------------------------------------------------------------
|  Defaults:
|   --min-mapq-bias:         1e-4
|   --min-end-distance-bias: 1e-4
|
|  VCF_Files:
|    "path/to/file1.vcf":
|      Output_File: "path/to/output1.vcf"
|      Options:
|        --min-strand-bias: 1e-4
|        --min-baseq-bias:  1e-4
|    "path/to/file2.vcf":
|      Output_File: "path/to/output2.vcf"
-------------------------------------------------------------------------------

Such a makefile can be parsed into a dictionary using YAML, but to help us
ensure that the makefile fits the expected layout (as above), we need to
specify the structure of the makefile.

Firstly, note that the options are specified twice, so we will make a re-usable
specification for those. In this case, this can accomplished like so:
-------------------------------------------------------------------------------
|  _SPECIFICATION_OF_OPTIONS = {
|    StringStartsWith("--") : Or(IsInt, IsFloat),
|  }
-------------------------------------------------------------------------------

or as so:
-------------------------------------------------------------------------------
|  _SPECIFICATION_OF_OPTIONS = IsDictOf(StringStartsWith("--"),
|                                       Or(IsInt, IsFloat))
-------------------------------------------------------------------------------

In either case, we require that the options be a dictionary with string keys
that start with "--", and that the values are either floats or integers. In
this case the two methods are equivalent, but normally the first method would
be preferred for more complex structures, while the second method is required
if different sub-structures are possible. For example, to require EITHER a
list of integers, or a dict of strings -> integers would have to be specified
as so:
-------------------------------------------------------------------------------
|  Or(IsListOf(IsInt), IsDictOf(IsStr, IsInt))
-------------------------------------------------------------------------------

Note that specification objects that do not take any parameters (IsInt, etc.)
do not need to be instantiated. Thus one can use both 'IsInt' or 'IsInt()',
whereas 'IsListOf', 'IsDictOf', etc. needs to be instantiated. This is purely
for convinience.


Having specified the expected structure of the options, we can specify the
remaining structure of the makefile:
-------------------------------------------------------------------------------
|  _MAKEFILE_SPECIFICATION = {
|    Defaults"  : _SPECIFICATION_OF_OPTIONS,
|
|    "VCF_Files" : {
|      Str : {
|        "Output_File" : IsStr,
|        "Options" : _SPECIFICATION_OF_OPTIONS,
|      }
|    }
|  }
-------------------------------------------------------------------------------

Finally, we can specify default values. Defaults can be specified for almost
all specification objects (excepting specifications for keys in dictionaries,
sub-specification for logical operators, and a couple of others). Let's suppose
that we always want a min/max depth set, even if the user did not include them
in the defaults:
-------------------------------------------------------------------------------
|  _SPECIFICATION_OF_OPTIONS = {
|    StringStartsWith("--") : Or(IsInt, IsFloat),
|    "--min-depth"          : IsInt(default = 8),
|    "--max-depth"          : IsInt(default = 100),
|  }
-------------------------------------------------------------------------------

These values would then be set, unless they were already set. Note that named
keys are given precedence above specification objects, when validating
key/value pairs. In other words, given this specification, the key
"--min-depth" is ALWAYS valid (even if it would fail StringStartsWith("--"),
and the value is ONLY checked against IsInt(default = 8).

Bringing all this together, we could then parse the a file containing the YAML
code shown above as follows:
-------------------------------------------------------------------------------
|  makefile = read_makefile("/path/to/makefile.yaml",
|                           _MAKEFILE_SPECIFICATION)
-------------------------------------------------------------------------------

which would yield the following dictionary:
-------------------------------------------------------------------------------
|  {'Defaults': {'--max-depth': 100,
|                '--min-depth': 8,
|                '--min-end-distance-bias': 0.001,
|                '--min-mapq-bias': 0.001},
|   'VCF_Files': {'path/to/file1.vcf':
|                   {'Options': {'--max-depth': 100,
|                    '--min-baseq-bias': 0.005,
|                    '--min-depth': 8,
|                    '--min-strand-bias': 0.005},
|                    'Output_File': 'path/to/output1.vcf'},
|                 'path/to/file2.vcf':
|                    {'Output_File': 'path/to/output1.vcf'}}},
-------------------------------------------------------------------------------

Unfortunately, the defaults are being applied to BOTH "Options" sub-trees,
which makes it impossible to tell which values are supposed to be over-ridden
for the files. To prevent this from happening, we can specify that defaults
should NOT be applied, by using the WithoutDefaults wrapper object:
-------------------------------------------------------------------------------
|  _MAKEFILE_SPECIFICATION = {
|    Defaults"  : _SPECIFICATION_OF_OPTIONS,
|
|    "VCF_Files" : {
|      Str : {
|        "Output_File" : IsStr,
|        "Options" : WithoutDefaults(_SPECIFICATION_OF_OPTIONS),
|      }
|    }
|  }
-------------------------------------------------------------------------------

Which yields the following structure following processing:
-------------------------------------------------------------------------------
|  {'Defaults': {'--max-depth': 100,
|                '--min-depth': 8,
|                '--min-end-distance-bias': 0.001,
|                '--min-mapq-bias': 0.001},
|   'VCF_Files': {'path/to/file1.vcf': {'Options': {'--min-baseq-bias': 0.005,
|                                                   '--min-strand-bias': 0.005},
|                                       'Output_File': 'path/to/output1.vcf'},
|                 'path/to/file2.vcf': {'Output_File': 'path/to/output2.vcf'}}}
-------------------------------------------------------------------------------

If the file contents does not match the expected structure, a MakefileError is
raised which describes the problem. For example, suppose that an "Output_File"
value has accidentically been left blank ('IsStr' requires a NON-EMPTY string):
-------------------------------------------------------------------------------
|  Makefile requirement not met at ...:
|    Expected value: a non-empty string
|    Observed value: ''
-------------------------------------------------------------------------------
"""
import copy
import logging

import paleomix.yaml

from paleomix.common.fileutils import fspath
from paleomix.common.utilities import group_by_pred


class MakefileError(RuntimeError):
    """Raised if a makefile is unreadable, or does not meet specifications."""


def read_makefile(filename, specification):
    """Reads and parses a makefile using the given specification."""
    try:
        with open(fspath(filename)) as handle:
            data = paleomix.yaml.safe_load(handle)
    except paleomix.yaml.YAMLError as error:
        raise MakefileError(error)

    return process_makefile(data, specification)


def process_makefile(data, specification, path=(), apply_defaults=True):
    """Validates a makefile and applies defaults to missing keys.

    Note that that default values are deep-copied before being set.
    """
    if isinstance(specification, WithoutDefaults):
        specification = specification.specification
        data = process_makefile(data, specification, path, apply_defaults=False)
    elif isinstance(specification, PreProcessMakefile):
        data, specification = specification(path, data)
        data = process_makefile(data, specification, path, apply_defaults)
    elif _is_spec(specification):
        _instantiate_spec(specification)(path, data)
    elif isinstance(data, (dict, type(None))) and isinstance(specification, dict):
        # A limitation of YAML is that empty subtrees are equal to None;
        # this check ensures that empty subtrees to be handled properly
        if data is None:
            data = {}

        _process_default_values(data, specification, path, apply_defaults)

        for cur_key in data:
            ref_key = _get_matching_spec_or_value(
                cur_key, specification, path + (cur_key,)
            )
            data[cur_key] = process_makefile(
                data[cur_key], specification[ref_key], path + (cur_key,), apply_defaults
            )
    elif isinstance(data, (list, type(None))) and isinstance(specification, list):
        if not all(_is_spec(spec) for spec in specification):
            raise TypeError(
                "Lists contains non-specification objects (%r): %r"
                % (_path_to_str(path), specification)
            )
        elif data is None:  # See comment above
            data = []

        specification = IsListOf(*specification)
        _instantiate_spec(specification)(path, data)
    elif not isinstance(specification, (dict, list)):
        raise TypeError(
            "Unexpected type in makefile specification at %r: %r!"
            % (_path_to_str(path), specification)
        )
    else:
        raise MakefileError(
            "Inconsistency between makefile specification and "
            "current makefile at %s:\n    Expected %s, "
            "found %s %r!"
            % (
                _path_to_str(path),
                type(specification).__name__,
                type(data).__name__,
                data,
            )
        )

    return data


###############################################################################
###############################################################################
# Unique 'value' used to specify that a MakefileSpec lacks a default value.
DEFAULT_NOT_SET = object()
# Unique 'value' used to specify that the user MUST supply a value
REQUIRED_VALUE = object()


class WithoutDefaults:
    """Wrapper object, that tells 'process_makefile' not to apply
    default values for the wrapped specification. See module docs
    for example usage.
    """

    def __init__(self, specification):
        self.specification = specification


class PreProcessMakefile:
    """Allows pre-processing of a part of a makefile prior to validation; when
    encountered, the object is called with the current value, and is expected
    to return a tuple containing (value, specification), which are then used
    subsequently. This allows transformation of fields for backwards
    compatibility.
    """

    def __call__(self, path, value):
        """Must return (value, specification) tuple."""
        raise NotImplementedError  # pragma: no coverage


class MakefileSpec:
    """Base-class for specifications, from which ALL specification
    objects are expected to derive. Sub-classes must implement the
    'meets_spec' function, which must return True or False depending
    on whether or not the given value meets the specification.
    """

    def __init__(self, description, default=DEFAULT_NOT_SET):
        """description -- A string describing the specification.
           default     -- A default value, or DEFAULT_NOT_SET if not used. If a
                          value is set, it is copied before being applied."""

        self.description = description
        self.default = default
        if (default not in (DEFAULT_NOT_SET, REQUIRED_VALUE)) and not self.meets_spec(
            default
        ):
            raise ValueError(
                (
                    "Default value does not meet requirements:\n"
                    "  Expected value: %s\n"
                    "  Observed value: %r\n"
                )
                % (description, default)
            )

    def __call__(self, path, value):
        if not self.meets_spec(value):
            raise MakefileError(
                (
                    "Makefile requirement not met at %r:\n"
                    "  Expected value: %s\n"
                    "  Observed value: %r\n"
                    "  Observed type:  %s"
                )
                % (_path_to_str(path), self.description, value, type(value).__name__)
            )

    def meets_spec(self, _value):
        """Return True if value meets the specification, False otherwise."""
        raise NotImplementedError


###############################################################################
###############################################################################
# Tests for basic types


class IsInt(MakefileSpec):
    """Require that the value is either an Int or a Long."""

    def __init__(self, description="an integer", default=DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, int) and not isinstance(value, bool)


class IsUnsignedInt(IsInt):
    """Require that the value is either an Int or a Long, and >= 0."""

    def __init__(self, description="an unsigned integer", default=DEFAULT_NOT_SET):
        IsInt.__init__(self, description, default)

    def meets_spec(self, value):
        return IsInt.meets_spec(self, value) and value >= 0


class IsFloat(MakefileSpec):
    """Require that the value is a float (does not cover integer types)."""

    def __init__(self, description="a float", default=DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, float)


class IsBoolean(MakefileSpec):
    """Require that the value is a boolean (True/False)."""

    def __init__(self, description="a boolean", default=DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, bool)


class IsStr(MakefileSpec):
    """Require that the value is a non-empty string."""

    def __init__(self, description=None, default=DEFAULT_NOT_SET, min_len=1):
        if description is None:
            if min_len == 0:
                description = "a string"
            elif min_len == 1:
                description = "a non-empty string"
            elif min_len >= 2:
                description = "a string at least %s characters long" % (min_len,)
            else:
                raise ValueError("min_len must be non-negative")

        self._min_len = min_len
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, str) and len(value) >= self._min_len


class IsNone(MakefileSpec):
    """Require that the value is None, typically signifying that
    the value was not set in the makefile."""

    def __init__(self, description="null or not set", default=DEFAULT_NOT_SET):
        if default is not DEFAULT_NOT_SET:
            raise NotImplementedError("IsNone does not support default values")
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return value is None


class ValueMissing(MakefileSpec):
    """Used to signify empty substructures in the makefile specification."""

    def __init__(self, description="no values"):
        MakefileSpec.__init__(self, description, DEFAULT_NOT_SET)

    def meets_spec(self, _value):
        return False


class DeprecatedOption(MakefileSpec):
    """Used to signify substructures that will eventually be removed."""

    def __init__(self, spec):
        self._spec = spec
        if not isinstance(spec, MakefileSpec):
            raise ValueError(spec)

        MakefileSpec.__init__(self, spec.description, spec.default)

    def __call__(self, path, value):
        self._spec(path, value)

        log = logging.getLogger(__name__)
        log.warning(
            "option has been deprecated and will be removed in the future: %s"
            % (_path_to_str(path),)
        )

    def meets_spec(self, value):
        return self._spec.meets_spec(value)


class RemovedOption(MakefileSpec):
    """Used to signify substructures that have been removed, and are hence ignored."""

    def __init__(self, description="removed settings"):
        MakefileSpec.__init__(self, description, DEFAULT_NOT_SET)

    def __call__(self, path, _value):
        log = logging.getLogger(__name__)
        log.warning(
            "option has been removed and no longer has any effect: %s"
            % (_path_to_str(path),)
        )

    def meets_spec(self, _value):
        return True


###############################################################################
###############################################################################
# BinaryOperators


class _BinaryOperator(MakefileSpec):
    """Base class for binary operations; takes a operation function which is
    assumed to take parameters (lvalue, rvalue), a rvalue to use when calling
    the function, and a description in the form 'operator {rvalue}' which is
    used to generate a human readable description of the specification.

    If list_kword is specified, the rvalue is assumed to be a sequence, and
    _list_values is used to convert it to human readable form.
    """

    def __init__(self, description, default, opfunc, rvalue, key=None, list_kword=None):
        self._operator = opfunc
        self._keyfunc = key
        self._rvalue = rvalue

        rvalue_repr = _list_values(rvalue, list_kword) if list_kword else rvalue
        description = description.format(rvalue=rvalue_repr)
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        if self._keyfunc is not None:
            value = self._keyfunc(value)
        return self._operator(value, self._rvalue)


class ValueIn(_BinaryOperator):
    def __init__(
        self,
        rvalues,
        key=None,
        description="value in {rvalue}",
        default=DEFAULT_NOT_SET,
    ):
        description = description.format(rvalue=_list_values(rvalues, "or"))
        _BinaryOperator.__init__(
            self,
            description=description,
            default=default,
            opfunc=self._in_operator,
            rvalue=rvalues,
            key=key,
        )

    def _in_operator(self, lvalue, rvalues):
        """Implements 'in' operator."""
        return _is_hashable(lvalue) and lvalue in rvalues


class ValuesIntersect(_BinaryOperator):
    def __init__(self, rvalues, key=None, description=None, default=DEFAULT_NOT_SET):
        if not description:
            description = "one or more of %s" % (_list_values(rvalues, "and"),)

        _BinaryOperator.__init__(
            self,
            description=description,
            default=default,
            opfunc=self._operator,
            rvalue=rvalues,
            key=key,
        )

    def _operator(self, lvalue, rvalues):
        try:
            return not isinstance(lvalue, dict) and bool(
                frozenset(lvalue).intersection(rvalues)
            )
        except TypeError:
            return False


class ValuesSubsetOf(_BinaryOperator):
    def __init__(self, rvalues, key=None, description=None, default=DEFAULT_NOT_SET):
        description = description or "subset of %s" % (_list_values(rvalues, "and"),)
        _BinaryOperator.__init__(
            self,
            description=description,
            default=default,
            opfunc=self._operator,
            rvalue=rvalues,
            key=key,
        )

    def _operator(self, lvalue, rvalues):
        try:
            return not isinstance(lvalue, dict) and bool(
                frozenset(lvalue).issubset(rvalues)
            )
        except TypeError:
            return False


###############################################################################
###############################################################################
# Logical operators


class _MultipleSpecs(MakefileSpec):
    """Base-class for logical operators for one or more specifications."""

    def __init__(
        self, specs, kwargs, name, prefix="", postfix="", join_by=" ", fmt="%s"
    ):
        self._specs = [_instantiate_spec(spec) for spec in specs]
        if not self._specs:
            raise ValueError("No specification given to %r" % (name.title(),))
        elif not all((spc.default is DEFAULT_NOT_SET) for spc in self._specs):
            raise ValueError(
                "Default values cannot be set in specs given to logical operators"
            )

        description = [(fmt % (spec.description,)) for spec in self._specs]
        description = "%s%s%s" % (prefix, join_by.join(description), postfix)
        default_value = kwargs.get("default", DEFAULT_NOT_SET)
        MakefileSpec.__init__(self, description, default_value)


class And(_MultipleSpecs):
    """Takes one or more specification objects, and requires that values meets
    all of these specifications. A default value may be set for the 'And'
    specification, but not for the specifications given to the 'And' object.
    """

    def __init__(self, *specs, **kwargs):
        _MultipleSpecs.__init__(self, specs, kwargs, "And", join_by=" and ", fmt="(%s)")

    def meets_spec(self, value):
        return all(spec.meets_spec(value) for spec in self._specs)


class Or(_MultipleSpecs):
    """Takes one or more specification objects, and requires that values meets
    at least one these specifications. A default value may be set for the 'Or'
    specification, but not for the specifications given to the 'Or' object.
    """

    def __init__(self, *specs, **kwargs):
        _MultipleSpecs.__init__(self, specs, kwargs, "Or", join_by=" or ", fmt="(%s)")

    def meets_spec(self, value):
        return any(spec.meets_spec(value) for spec in self._specs)


class Not(_MultipleSpecs):
    """Takes a single specification object, and requires that values do NOT
    meet this specification. A default value may be set for the 'Not'
    specification, but not for the specifications given to the 'Not' object.
    """

    def __init__(self, spec, **kwargs):
        _MultipleSpecs.__init__(self, [spec], kwargs, "Not", prefix="not ", fmt="(%s)")

    def meets_spec(self, value):
        return not self._specs[0].meets_spec(value)


###############################################################################
###############################################################################
# String operators
#
# In addition to providing string-specific operators (is uppercase, ends/starts
# with), "in" and set operators are provided which do case-insensitive
# comparsions. For case-sensitive operations, use the Value* specifications.


class StringIn(_BinaryOperator):
    """Require that values are found in a set of values. For strings, the
    comparison is done in a case-insensitive. For case-sensitive comparisons,
    see 'ValueIn'.
    """

    def __init__(
        self,
        rvalues,
        key=None,
        description="one of {rvalue}, case-insentive",
        default=DEFAULT_NOT_SET,
    ):
        description = description.format(rvalue=_list_values(rvalues, "or"))
        rvalues = frozenset(map(_safe_coerce_to_lowercase, rvalues))

        _BinaryOperator.__init__(
            self, description, default, self._string_in_operator, rvalues
        )

    @classmethod
    def _string_in_operator(cls, lvalue, rvalues):
        """Implements case-insensitive 'in' operator."""
        if not _is_hashable(lvalue):
            return False

        return _safe_coerce_to_lowercase(lvalue) in rvalues


class StringStartsWith(IsStr):
    """Require that the value is a string with given prefix."""

    def __init__(self, prefix, default=DEFAULT_NOT_SET):
        assert prefix and isinstance(prefix, str)
        self._prefix = prefix
        description = "a string with prefix %r" % (prefix,)
        IsStr.__init__(self, description, default)

    def meets_spec(self, value):
        return super(StringStartsWith, self).meets_spec(value) and value.startswith(
            self._prefix
        )


class StringEndsWith(IsStr):
    """Require that the value is a string with given postfix."""

    def __init__(self, postfix, default=DEFAULT_NOT_SET):
        assert postfix and isinstance(postfix, str)
        self._postfix = postfix
        description = "a string with postfix %r" % (postfix,)
        IsStr.__init__(self, description, default)

    def meets_spec(self, value):
        return super(StringEndsWith, self).meets_spec(value) and value.endswith(
            self._postfix
        )


###############################################################################
###############################################################################
# Tests for collections


class IsListOf(_MultipleSpecs):
    """Require that the value is a list, the contents of which matches one or
    more of the provided specifications; if no default value (ie. a non-empty
    list) is required, then using the following syntax is preferred:
      [IsType1, IsType2, ...]
    This is equivalent to the following:
      IsListOf(IsType1, IsType2, ...)
    """

    def __init__(self, *specs, **kwargs):
        _MultipleSpecs.__init__(
            self,
            specs,
            kwargs,
            "IsListOf",
            prefix="[",
            postfix=", ...]",
            join_by=" or ",
            fmt="(%s)",
        )

    def meets_spec(self, value):
        if not isinstance(value, list):
            return False

        return all(
            any(spec.meets_spec(lstvalue) for spec in self._specs) for lstvalue in value
        )


class IsDictOf(MakefileSpec):
    """Require that the value is a list, the keys/values of which matches
    the specifications provided for keys/values; if no default value (ie. a
    dictioanry) is required, then using the following syntax is preferred:
      {IsType1: IsType2}

    This is equivalent to the following:
      IsDictOf(IsType1, IsType2)
    but also allows multiple type-pairs to be specified.
    """

    def __init__(self, key_spec, value_spec, default=DEFAULT_NOT_SET):
        self._key_spec = _instantiate_spec(key_spec)
        self._value_spec = _instantiate_spec(value_spec)
        if self._key_spec.default is not DEFAULT_NOT_SET:
            raise ValueError("Default values cannot be set in key-specs")
        elif self._value_spec.default is not DEFAULT_NOT_SET:
            raise ValueError("Default values cannot be set in value-specs")

        description = "{(%s) : (%s)}" % (
            self._key_spec.description,
            self._value_spec.description,
        )
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        if not isinstance(value, dict):
            return False

        for (key, value) in value.items():
            if not (
                self._key_spec.meets_spec(key) and self._value_spec.meets_spec(value)
            ):
                return False

        return True


###############################################################################
###############################################################################
# Helper functions


def _is_hashable(value):
    try:
        hash(value)
        return True
    except TypeError:
        return False


def _is_spec(spec):
    """Returns true if 'spec' is a specification instance or class."""
    if isinstance(spec, MakefileSpec):
        return True
    elif isinstance(spec, type) and issubclass(spec, MakefileSpec):
        return True
    return False


def _instantiate_spec(spec):
    """Takes a specification instance or class, and returns an instance."""
    if isinstance(spec, MakefileSpec):
        return spec
    elif isinstance(spec, type) and issubclass(spec, MakefileSpec):
        return spec()
    else:
        raise TypeError("Specifications must derive from 'MakefileSpec'")


def _safe_coerce_to_lowercase(value):
    """Returns strings as lowercase, and any other types of value unchanged."""
    if isinstance(value, str):
        return value.lower()
    return value


def _list_values(values, sep):
    """Returns list of values as '[values[0], values[1], ..., sep values[-1]]':

    $ _list_values([1, 2, 3], "and")
    "[1, 2, and 3]"
    """
    values = list(map(repr, values))
    if len(values) > 2:
        values = (", ".join(values[:-1]) + ",", values[-1])
    if len(values) == 2:
        values = (" ".join((values[0], sep, values[1])),)

    return values[0]


def _get_summary_spec(specs_or_keys):
    """Returns a specification object that may be used to describe a set of
    requirements. This is used if a key or value does not match the possible
    specs, thereby describing the set of allowed values.
    """
    specs, keys = group_by_pred(_is_spec, specs_or_keys)
    if specs and keys:
        return Or(ValueIn(keys, description="key in {rvalue}"), *specs)
    elif specs:
        return Or(*specs)
    elif keys:
        return ValueIn(keys, description="key in {rvalue}")
    return ValueMissing()


def _get_matching_spec_or_value(value, specs, path):
    """Returns the specification object or value that matches the observed
    value; specs may be a list of specification objects and/or constant values
    allowed by the makefile. If no matching specification or value is found,
    an MakefileError is raised.
    """
    if value in specs:
        return value

    for spec in specs:
        if _is_spec(spec) and _instantiate_spec(spec).meets_spec(value):
            return spec

    # No matching key or spec; create combined spec to raise error message
    _get_summary_spec(specs)(path, value)
    assert False  # pragma: no coverage


def _process_default_values(data, specification, path, apply_defaults):
    """Checks a subtree against a specification, verifies that required values
    have been set, and (optionally) sets values for keys where defaults have
    been specified.
    """

    for cur_key in specification:
        if (not _is_spec(cur_key)) and (cur_key not in data):
            default_value = specification[cur_key]
            default_value_from_spec = False

            while isinstance(default_value, PreProcessMakefile):
                data, default_value = default_value(path, data)

            if _is_spec(default_value):
                default_value = _instantiate_spec(default_value)
                if default_value.default is DEFAULT_NOT_SET:
                    continue
                elif default_value.default is REQUIRED_VALUE:
                    raise MakefileError(
                        "A value MUST be supplified for %r"
                        % (_path_to_str(path + (cur_key,)))
                    )
                default_value = default_value.default
                default_value_from_spec = True

            if apply_defaults and not isinstance(
                default_value, (PreProcessMakefile, WithoutDefaults)
            ):
                if isinstance(default_value, dict):
                    # Setting of values in the dict will be accomplished
                    # in subsequent calls to _process_default_values
                    default_value = {}
                elif isinstance(default_value, list):
                    # Lists of specs defaults to empty lists
                    if not default_value_from_spec:
                        default_value = []

                # Prevent clobbering of values when re-using sub-specs
                data[cur_key] = copy.deepcopy(default_value)


def _path_to_str(path):
    """Converts a path (tuple of strings) to a printable string."""
    return " :: ".join(str(field) for field in path)


CLI_PARAMETERS = Or(IsListOf(IsStr, IsInt, IsFloat), Or(IsStr, IsInt, IsFloat, IsNone))
