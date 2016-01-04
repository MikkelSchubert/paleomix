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
|  {'Makefile': {'Defaults': {'--max-depth': 100,
|                             '--min-depth': 8,
|                             '--min-end-distance-bias': 0.001,
|                             '--min-mapq-bias': 0.001},
|                'VCF_Files': {'path/to/file1.vcf':
|                                {'Options': {'--max-depth': 100,
|                                 '--min-baseq-bias': 0.005,
|                                 '--min-depth': 8,
|                                 '--min-strand-bias': 0.005},
|                                 'Output_File': 'path/to/output1.vcf'},
|                              'path/to/file2.vcf':
|                                 {'Output_File': 'path/to/output1.vcf'}}},
|   'Statistics': {'Filename': 'makefile.yaml',
|                  'Hash': 'c0138fd4ffcbbc0dff2c82c6e7595ec38b01f532',
|                  'MTime': '2013-08-13 10:22:46.000000 '}}
-------------------------------------------------------------------------------

Note the actual contents of the makefile is found in the sub-dictionary
"Makefile", while the sub-dictionary "Statistics" contains various information
about the file itself.

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
|  {'Makefile': {'Defaults': {'--max-depth': 100,
|                             '--min-depth': 8,
|                             '--min-end-distance-bias': 0.001,
|                             '--min-mapq-bias': 0.001},
|                'VCF_Files': {'path/to/file1.vcf':
|                                {'Options': {'--min-baseq-bias': 0.005,
|                                             '--min-strand-bias': 0.005},
|                                 'Output_File': 'path/to/output1.vcf'},
|                              'path/to/file2.vcf':
|                                {'Output_File': 'path/to/output2.vcf'}}},
|   'Statistics': {'Filename': 'makefile.yaml',
|                  'Hash': 'c0138fd4ffcbbc0dff2c82c6e7595ec38b01f532',
|                  'MTime': '2013-08-13 10:22:46.000000 '}}
-------------------------------------------------------------------------------

If the file contents does not match the expected structure, a MakefileError is
raised which describes the problem. For example, suppose that an "Output_File"
value has accidentically been left blank ('IsStr' requires a NON-EMPTY string):
-------------------------------------------------------------------------------
|  Makefile requirement not met at ...:
|    Expected value(s): a non-empty string
|    Observed value(s): ''
-------------------------------------------------------------------------------
"""
import os
import copy
import types
import hashlib
import datetime
import operator

import paleomix.yaml
from paleomix.common.utilities import group_by_pred


class MakefileError(RuntimeError):
    """Raised if a makefile is unreadable, or does not meet specifications."""


def read_makefile(filename, specification):
    """Reads and parses a makefile using the given specification.

    Returns a dictionary of the form
      {
         "Makefile": <parsed makefile>,
         "Statistics": {
            "Filename": <filename>,
            "Hash": <SHA1 hash of makefile>,
            "MTime": <Modification time of makefile>,
         }
      }
    """
    try:
        with open(filename) as makefile:
            string = makefile.read()
            data = paleomix.yaml.safe_load(string)
    except paleomix.yaml.error.YAMLError, error:
        raise MakefileError(error)

    mtime = os.path.getmtime(os.path.realpath(filename))
    mtime_str = datetime.datetime.fromtimestamp(mtime).strftime("%F %T")
    return {"Makefile": process_makefile(data, specification),
            "Statistics": {"Filename": filename,
                           "Hash": hashlib.sha1(string).hexdigest(),
                           "MTime": mtime_str}}


def process_makefile(data, specification, path=("root",), apply_defaults=True):
    """Validates a makefile and applies defaults to missing keys.

    Note that that default values are deep-copied before being set.
    """
    if isinstance(specification, WithoutDefaults):
        specification = specification.specification
        data = process_makefile(data, specification, path,
                                apply_defaults=False)
    elif isinstance(specification, PreProcessMakefile):
        data, specification = specification(path, data)
        data = process_makefile(data, specification, path, apply_defaults)
    elif _is_spec(specification):
        _instantiate_spec(specification)(path, data)
    elif isinstance(data, (dict, types.NoneType)) \
            and isinstance(specification, dict):
        # A limitation of YAML is that empty subtrees are equal to None;
        # this check ensures that empty subtrees to be handled properly
        if data is None:
            data = {}

        _process_default_values(data, specification, path, apply_defaults)

        for cur_key in data:
            ref_key = _get_matching_spec_or_value(cur_key,
                                                  specification,
                                                  path + (cur_key,))
            data[cur_key] = process_makefile(data[cur_key],
                                             specification[ref_key],
                                             path + (cur_key,),
                                             apply_defaults)
    elif isinstance(data, (list, types.NoneType)) \
            and isinstance(specification, list):
        if not all(_is_spec(spec) for spec in specification):
            raise TypeError("Lists contains non-specification objects (%r): %r"
                            % (_path_to_str(path), specification))
        elif data is None:  # See comment above
            data = []

        specification = IsListOf(*specification)
        _instantiate_spec(specification)(path, data)
    elif not isinstance(specification, (dict, list)):
        raise TypeError("Unexpected type in makefile specification at %r: %r!"
                        % (_path_to_str(path), specification))
    else:
        raise MakefileError("Inconsistency between makefile specification and "
                            "current makefile at %s:\n    Expected %s, "
                            "found %s %r!" % (_path_to_str(path),
                                           type(specification).__name__,
                                           type(data).__name__,
                                           data))

    return data


###############################################################################
###############################################################################
# Unique 'value' used to specify that a MakefileSpec lacks a default value.
DEFAULT_NOT_SET = object()
# Unique 'value' used to specify that the user MUST supply a value
REQUIRED_VALUE = object()


class WithoutDefaults(object):
    """Wrapper object, that tells 'process_makefile' not to apply
    default values for the wrapped specification. See module docs
    for example usage.
    """

    def __init__(self, specification):
        self.specification = specification


class PreProcessMakefile(object):
    """Allows pre-processing of a part of a makefile prior to validation; when
    encountered, the object is called with the current value, and is expected
    to return a tuple containing (value, specification), which are then used
    subsequently. This allows transformation of fields for backwards
    compatibility.
    """

    def __call__(self, path, value):
        """Must return (value, specification) tuple."""
        raise NotImplementedError


class MakefileSpec(object):
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
        if (default not in (DEFAULT_NOT_SET, REQUIRED_VALUE)) \
                and not self.meets_spec(default):
            raise ValueError(("Default value does not meet requirements:\n"
                              "  Expected value(s): %s\n"
                              "  Observed value(s): %r\n")
                             % (description, default))

    def __call__(self, path, value):
        if not self.meets_spec(value):
            raise MakefileError(("Makefile requirement not met at %r:\n"
                                 "  Expected value(s): %s\n"
                                 "  Observed value(s): %r\n"
                                 "  Observed type:     %s")
                                % (_path_to_str(path), self.description,
                                   value, type(value).__name__))

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
        return isinstance(value, (types.IntType, types.LongType)) \
            and not isinstance(value, types.BooleanType)


class IsUnsignedInt(IsInt):
    """Require that the value is either an Int or a Long, and >= 0."""

    def __init__(self, description="an unsigned integer",
                 default=DEFAULT_NOT_SET):
        IsInt.__init__(self, description, default)

    def meets_spec(self, value):
        return IsInt.meets_spec(self, value) & (value >= 0)


class IsFloat(MakefileSpec):
    """Require that the value is a float (does not cover integer types)."""

    def __init__(self, description="a float", default=DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, types.FloatType)


class IsBoolean(MakefileSpec):
    """Require that the value is a boolean (True/False)."""

    def __init__(self, description="a boolean", default=DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, types.BooleanType)


class IsStr(MakefileSpec):
    """Require that the value is a non-empty string."""

    def __init__(self, description="a non-empty string",
                 default=DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, types.StringTypes) and value


class IsNone(MakefileSpec):
    """Require that the value is None, typically signifying that
    the value was not set in the makefile."""

    def __init__(self, description="None or not set", default=DEFAULT_NOT_SET):
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

    def __init__(self, description, default, opfunc, rvalue, key=None,
                 list_kword=None):
        self._operator = opfunc
        self._keyfunc = key
        self._rvalue = rvalue

        repr_func = repr
        if list_kword is not None:
            repr_func = lambda value: _list_values(value, list_kword)
        description = description.format(rvalue=repr_func(rvalue))
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        if self._keyfunc is not None:
            value = self._keyfunc(value)
        return self._operator(value, self._rvalue)


def _create_binary_operator(operator_func, description, list_kword=None):
    """Creates and returns a BinaryOperator class based on the given
    operator_func function, which is assumed to be a function taking two
    arguments (lvalue, rvalue) and returning a boolean value.
    """

    class _BinaryOperatorImpl(_BinaryOperator):
        """Implements a binary operator specfication."""

        def __init__(self, rvalue, key=None, description=description,
                     default=DEFAULT_NOT_SET):
            _BinaryOperator.__init__(self, description, default, operator_func,
                                     rvalue, key, list_kword)
    return _BinaryOperatorImpl


def _create_set_operator(operator_func, description):
    """Creates and returns a BinaryOperator designed to operate on sets of
    values. Thus, values in the makefile are expected to be either lists or
    strings (for case sensitive operations).
    """

    def _operator(lvalue, rvalue):
        """Operator function for set based operations."""
        if not isinstance(lvalue, (types.ListType,) + types.StringTypes):
            return False

        return bool(operator_func(frozenset(lvalue), rvalue))

    description = "%s {rvalue}" % (description,)
    return _create_binary_operator(_operator, description, "or")


ValueLT = _create_binary_operator(operator.lt, "value < {rvalue}")
ValueLE = _create_binary_operator(operator.le, "value <= {rvalue}")
ValueEQ = _create_binary_operator(operator.eq, "value = {rvalue}")
ValueGE = _create_binary_operator(operator.ge, "value >= {rvalue}")
ValueGT = _create_binary_operator(operator.gt, "value > {rvalue}")
ValueIn = _create_binary_operator(lambda lvalue, rvalue: lvalue in rvalue,
                                  "value in {rvalue}", "or")
ValuesIntersect = _create_set_operator(frozenset.intersection, "contains")
ValuesSubsetOf = _create_set_operator(frozenset.issubset, "subset of")


###############################################################################
###############################################################################
# Logical operators

class _MultipleSpecs(MakefileSpec):  # pylint: disable=W0223
    """Base-class for logical operators for one or more specifications."""

    def __init__(self, specs, kwargs, name, prefix="", postfix="",
                 join_by=" ", fmt="%s"):
        self._specs = [_instantiate_spec(spec) for spec in specs]
        if not self._specs:
            raise ValueError("No specification given to %r" % (name.title(),))
        elif not all((spc.default is DEFAULT_NOT_SET) for spc in self._specs):
            raise ValueError("Default values cannot be set in specs given to "
                             "logical operators")

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
        _MultipleSpecs.__init__(self, specs, kwargs, "And",
                                join_by=" and ", fmt="(%s)")

    def meets_spec(self, value):
        return all(spec.meets_spec(value) for spec in self._specs)


class Or(_MultipleSpecs):
    """Takes one or more specification objects, and requires that values meets
    at least one these specifications. A default value may be set for the 'Or'
    specification, but not for the specifications given to the 'Or' object.
    """

    def __init__(self, *specs, **kwargs):
        _MultipleSpecs.__init__(self, specs, kwargs, "Or",
                                join_by=" or ", fmt="(%s)")

    def meets_spec(self, value):
        return any(spec.meets_spec(value) for spec in self._specs)


class Xor(_MultipleSpecs):
    """Takes two specification objects, and requires that values meets ONE and
    ONLY ONE of these specifications. A default value may be set for the 'Xor'
    specification, but not for the specifications given to the 'Xor' object.
    """

    def __init__(self, *specs, **kwargs):
        if len(specs) != 2:
            raise ValueError("'Xor' takes exactly 2 specifications, not %i"
                             % (len(specs),))

        _MultipleSpecs.__init__(self, specs, kwargs, "Xor",
                                join_by=" xor ", fmt="(%s)")

    def meets_spec(self, value):
        return operator.xor(*(spec.meets_spec(value) for spec in self._specs))


class Not(_MultipleSpecs):
    """Takes a single specification object, and requires that values do NOT
    meet this specification. A default value may be set for the 'Not'
    specification, but not for the specifications given to the 'Not' object.
    """

    def __init__(self, spec, **kwargs):
        _MultipleSpecs.__init__(self, [spec], kwargs, "Not",
                                prefix="not ", fmt="(%s)")

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

    def __init__(self, rvalues, key=None,
                 description="one of {rvalue}, case-insentive",
                 default=DEFAULT_NOT_SET):
        description = description.format(rvalue=_list_values(rvalues, "or"))
        rvalues = frozenset(map(_safe_coerce_to_lowercase, rvalues))

        _BinaryOperator.__init__(self, description, default,
                                 self._string_in_operator, rvalues)

    @classmethod
    def _string_in_operator(cls, lvalue, rvalues):
        """Implements case-insensitive 'in' operator."""
        return _safe_coerce_to_lowercase(lvalue) in rvalues


class _StrSetOperator(_BinaryOperator):
    """Base class for set operations involving case-insensitive strings."""

    def __init__(self, description, default, opfunc, rvalues, key=None):
        rvalues = frozenset(map(_safe_coerce_to_lowercase, rvalues))
        _BinaryOperator.__init__(self, description, default, opfunc, rvalues,
                                 key)

    def meets_spec(self, value):
        if not isinstance(value, (types.ListType,) + types.StringTypes):
            return False

        lvalues = frozenset(map(_safe_coerce_to_lowercase, value))
        return _BinaryOperator.meets_spec(self, lvalues)


class StringsIntersect(_StrSetOperator):
    """Require that a set of values overlap with a pre-defined set of values
    (as set in the constructor). For strings, values are compared in a case-
    insensitive manner. For case-sensitive comparisons, see 'ValuesIntersect'.
    """

    def __init__(self, rvalue, key=None,
                 description="contains {rvalue}, case-insentive",
                 default=DEFAULT_NOT_SET):
        description = description.format(rvalue=_list_values(rvalue, "and/or"))

        _StrSetOperator.__init__(self, description, default,
                                 self._string_intersection_operator,
                                 frozenset(rvalue), key)

    @classmethod
    def _string_intersection_operator(cls, lvalue, rvalues):
        """Implements case-insensitive 'intersect' operator."""
        return bool(frozenset(lvalue).intersection(rvalues))


class StringsSubsetOf(_StrSetOperator):
    """Require that a set of values are a subset of a pre-defined set of values
    (as set in the constructor). For strings, values are compared in a case-
    insensitive manner. For case-sensitive comparisons, see 'ValuesSubsetOf'.

    Note that empty sets are always considered to be a subset of the
    pre-defined set.
    """

    def __init__(self, rvalue, key=None,
                 description="subset of {rvalue}, case-insentive",
                 default=DEFAULT_NOT_SET):
        description = description.format(rvalue=_list_values(rvalue, "and"))

        _StrSetOperator.__init__(self, description, default,
                                 self._operator_func, frozenset(rvalue), key)

    @classmethod
    def _operator_func(cls, lvalue, rvalue):
        """Operator implementation."""
        return bool(frozenset(lvalue).issubset(rvalue))


class StringIsUppercase(IsStr):
    """Require that the value is a uppercase, non-empty string."""

    def __init__(self, default=DEFAULT_NOT_SET):
        IsStr.__init__(self, "an uppercase non-empty string", default)

    def meets_spec(self, value):
        return IsStr().meets_spec(value) and value.isupper()


class StringStartsWith(IsStr):
    """Require that the value is a string with given prefix."""

    def __init__(self, prefix, default=DEFAULT_NOT_SET):
        assert prefix and isinstance(prefix, types.StringTypes)
        self._prefix = prefix
        description = "a string with the prefix %r" % (prefix,)
        IsStr.__init__(self, description, default)

    def meets_spec(self, value):
        return IsStr.meets_spec(self, value) \
            and value.startswith(self._prefix)


class StringEndsWith(IsStr):
    """Require that the value is a string with given postfix."""

    def __init__(self, postfix, default=DEFAULT_NOT_SET):
        assert postfix and isinstance(postfix, types.StringTypes)
        self._postfix = postfix
        description = "a string with the postfix %r" % (postfix,)
        IsStr.__init__(self, description, default)

    def meets_spec(self, value):
        return IsStr.meets_spec(self, value) and value.endswith(self._postfix)


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
        _MultipleSpecs.__init__(self, specs, kwargs, "IsListOf",
                                prefix="[", postfix=", ...]",
                                join_by=" or ", fmt="(%s)")

    def meets_spec(self, value):
        if not isinstance(value, types.ListType):
            return False

        return all(any(spec.meets_spec(lstvalue) for spec in self._specs)
                   for lstvalue in value)


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

        description = "{(%s) : (%s)}" \
            % (self._key_spec.description, self._value_spec.description)
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        if not isinstance(value, types.DictType):
            return False

        for (key, value) in value.iteritems():
            if not (self._key_spec.meets_spec(key)
                    and self._value_spec.meets_spec(value)):
                return False

        return True


###############################################################################
###############################################################################
# Helper functions

def _is_spec(spec):
    """Returns true if 'spec' is a specification instance or class."""
    if isinstance(spec, MakefileSpec):
        return True
    elif isinstance(spec, types.TypeType) and issubclass(spec, MakefileSpec):
        return True
    return False


def _instantiate_spec(spec):
    """Takes a specification instance or class, and returns an instance."""
    if isinstance(spec, MakefileSpec):
        return spec
    elif isinstance(spec, types.TypeType) and issubclass(spec, MakefileSpec):
        return spec()
    else:
        raise TypeError("Specifications must derive from 'MakefileSpec'")


def _safe_coerce_to_lowercase(value):
    """Returns strings as lowercase, and any other types of value unchanged."""
    if isinstance(value, types.StringTypes):
        return value.lower()
    return value


def _list_values(values, sep):
    """Returns list of values as '[values[0], values[1], ..., sep values[-1]]':

    $ _list_values([1, 2, 3], "and")
    "[1, 2, and 3]"
    """
    values = map(repr, values)
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
                    raise MakefileError("A value MUST be supplified for %r"
                                        % (_path_to_str(path + (cur_key,))))
                default_value = default_value.default
                default_value_from_spec = True

            if apply_defaults \
                    and not isinstance(default_value, (PreProcessMakefile,
                                                       WithoutDefaults)):
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
    return ":".join(str(field) for field in path)


CLI_PARAMETERS = Or(IsListOf(IsStr, IsInt, IsFloat),
                    Or(IsStr, IsInt, IsFloat, IsNone))
