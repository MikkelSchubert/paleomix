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
"""Generalized methods for parsing/validating "makefiles" in YAML format.

The following example will use a imagined makefile for 'vcf_filter', which
takes a set of input files, each of which is assigned an output file, and
each of which may have a set of filters (in addition to a set of defaults):

The following example can also be found in pypeline/examples/pypeline/makefile.

Example makefile in YAML format:
--------------------------------------------------------------------------------
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
--------------------------------------------------------------------------------

Such a makefile can be parsed into a dictionary using YAML, but to help us
ensure that the makefile fits the expected layout (as above), we need to specify
the structure of the makefile.

Firstly, note that the options are specified twice, so we will start by making a
re-usable specification for those. In this case, this can accomplished like so:
--------------------------------------------------------------------------------
|  _SPECIFICATION_OF_OPTIONS = {
|    StringStartsWith("--") : Or(IsInt, IsFloat),
|  }
--------------------------------------------------------------------------------

or as so:
--------------------------------------------------------------------------------
|  _SPECIFICATION_OF_OPTIONS = IsDictOf(StringStartsWith("--"),
|                                       Or(IsInt, IsFloat))
--------------------------------------------------------------------------------

In either case, we require that the options be a dictionary with string keys
that start with "--", and that the values are either floats or integers. In this
case the two methods are equivalent, but normally the first method would be
preferred for more complex structures, while the second method is required if a
different sub-structures are possible. For example, to require EITHER a list of
integers, or a dict of strings -> integers would have to be specified as so:
--------------------------------------------------------------------------------
|  Or(IsListOf(IsInt), IsDictOf(IsStr, IsInt))
--------------------------------------------------------------------------------

Note that specification objects that do not take any parameters (IsInt, etc.) do
not need to be instantiated. Thus one can use both 'IsInt' or 'IsInt()', whereas
'IsListOf', 'IsDictOf', etc. needs to be instantiated. This is purely for
convinience.


Having specified the expected structure of the options, we can specify the
remaining structure of the makefile:
--------------------------------------------------------------------------------
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
--------------------------------------------------------------------------------

Finally, we can specify default values. Defaults can be specified for almost all
specification objects (excepting specifications for keys in dictionaries, sub-
specification for logical operators, and a couple of others). Let's suppose that
we always want a min/max depth set, even if the user did not include them in the
defaults:
--------------------------------------------------------------------------------
|  _SPECIFICATION_OF_OPTIONS = {
|    StringStartsWith("--") : Or(IsInt, IsFloat),
|    "--min-depth"          : IsInt(default = 8),
|    "--max-depth"          : IsInt(default = 100),
|  }
--------------------------------------------------------------------------------

These values would then be set, unless they were already set. Note that named
keys are given precedence above specification objects, when validating key/value
pairs. In other words, given this specification, the key "--min-depth" is ALWAYS
valid (even if it would fail StringStartsWith("--"), and the value is ONLY
checked against IsInt(default = 8).

Bringing all this together, we could then parse the a file containing the YAML
code shown above as follows:
--------------------------------------------------------------------------------
|  makefile = read_makefile("/path/to/makefile.yaml",
|                           _MAKEFILE_SPECIFICATION)
--------------------------------------------------------------------------------

which would yield the following dictionary:
--------------------------------------------------------------------------------
|  {'Makefile': {'Defaults': {'--max-depth': 100,
|                             '--min-depth': 8,
|                             '--min-end-distance-bias': 0.001,
|                             '--min-mapq-bias': 0.001},
|                'VCF_Files': {'path/to/file1.vcf': {'Options': {'--max-depth': 100,
|                                                                '--min-baseq-bias': 0.005,
|                                                                '--min-depth': 8,
|                                                                '--min-strand-bias': 0.005},
|                                                    'Output_File': 'path/to/output1.vcf'},
|                              'path/to/file2.vcf': {'Output_File': 'path/to/output1.vcf'}}},
|   'Statistics': {'Filename': 'makefile.yaml',
|                  'Hash': 'c0138fd4ffcbbc0dff2c82c6e7595ec38b01f532',
|                  'MTime': '2013-08-13 10:22:46.000000 '}}
--------------------------------------------------------------------------------

Note the actual contents of the makefile is found in the sub-dictionary
"Makefile", while the sub-dictionary "Statistics" contains various information
about the file itself.

Unfortunately, the defaults are being applied to BOTH "Options" sub-trees, which
makes it impossible to tell which values are supposed to be over-ridden for the
files. To prevent this from happening, we can specify that defaults should NOT
be applied, by using the WithoutDefaults wrapper object:
--------------------------------------------------------------------------------
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
--------------------------------------------------------------------------------

Which yields the following structure following processing:
--------------------------------------------------------------------------------
|  {'Makefile': {'Defaults': {'--max-depth': 100,
|                             '--min-depth': 8,
|                             '--min-end-distance-bias': 0.001,
|                             '--min-mapq-bias': 0.001},
|                'VCF_Files': {'path/to/file1.vcf': {'Options': {'--min-baseq-bias': 0.005,
|                                                                '--min-strand-bias': 0.005},
|                                                    'Output_File': 'path/to/output1.vcf'},
|                              'path/to/file2.vcf': {'Output_File': 'path/to/output2.vcf'}}},
|   'Statistics': {'Filename': 'makefile.yaml',
|                  'Hash': 'c0138fd4ffcbbc0dff2c82c6e7595ec38b01f532',
|                  'MTime': '2013-08-13 10:22:46.000000 '}}
--------------------------------------------------------------------------------

If the file contents does not match the expected structure, a MakefileError is
raised which describes the problem. For example, suppose that an "Output_File"
value has accidentically been left blank ('IsStr' requires a NON-EMPTY string):
--------------------------------------------------------------------------------
|  Makefile requirement not met at 'root:VCF_Files:path/to/file1.vcf:Output_File':
|    Expected value(s): a non-empty string
|    Observed value(s): ''
--------------------------------------------------------------------------------
"""
import yaml

import os
import copy
import types
import hashlib
import datetime
import operator

from pypeline.common.utilities import group_by_pred



class MakefileError(RuntimeError):
    pass


def read_makefile(filename, specification):
    try:
        with open(filename) as makefile:
            string = makefile.read()
            data = yaml.safe_load(string)
    except RuntimeError, error:
        raise MakefileError(error)

    mtime = os.path.getmtime(os.path.realpath(filename))
    return {"Makefile"   : process_makefile(data, specification),
            "Statistics" : {"Filename" : filename,
                            "Hash"     : hashlib.sha1(string).hexdigest(),
                            "MTime"    : datetime.datetime.fromtimestamp(mtime).strftime("%F %T.%f %z")}}


def process_makefile(data, specification, path = ("root",), apply_defaults = True):
    """Validates a makefile and applies defaults to missing keys.

    Note that that default values are deep-copied before being set."""
    if isinstance(specification, WithoutDefaults):
        specification = specification.specification
        process_makefile(data, specification, path, apply_defaults = False)
    elif _is_spec(specification):
        _instantiate_spec(specification)(path, data)
    elif isinstance(data, dict) and isinstance(specification, dict):
        if apply_defaults:
            _set_default_values(data, specification)

        for cur_key in data:
            ref_key = _get_matching_spec_or_value(cur_key, specification, path + (cur_key,))
            process_makefile(data[cur_key], specification[ref_key], path + (cur_key,), apply_defaults)
    elif not isinstance(specification, dict):
        raise TypeError("Unexpected type in makefile specification at %r: %r!" \
                        % (_path_to_str(path), specification))
    else:
        raise MakefileError(("Inconsistency between makefile format and current makefile at %s\n"
                             "Expected dict, found %r!") \
                            % (_path_to_str(path), data))

    return data




################################################################################
################################################################################
# Unique 'value' used to specify that a MakefileSpec does not have a default value.
DEFAULT_NOT_SET = object()


class WithoutDefaults:
    """Wrapper object, that tells 'process_makefile' not to apply
    default values for the wrapped specification. See module docs
    for example usage."""
    def __init__(self, specification):
        self.specification = specification


class MakefileSpec:
    """Base-class for specifications, from which ALL specification
    objects are expected to derive. Sub-classes must implement the
    'meets_spec' function, which must return True or False depending
    on whether or not the given value meets the specification."""

    def __init__(self, description, default = DEFAULT_NOT_SET):
        """description -- A string describing the specification.
           default     -- A default value, or DEFAULT_NOT_SET if not used. If a
                          value is set, it is copied before being applied."""

        self.description = description
        self.default     = default
        if (default is not DEFAULT_NOT_SET) and not self.meets_spec(default):
            raise ValueError(("Default value does not meet requirements:\n"
                              "  Expected value(s): %s\n"
                              "  Observed value(s): %r\n")
                              % (description, default))

    def __call__(self, path, value):
        if not self.meets_spec(value):
            raise MakefileError(("Makefile requirement not met at %r:\n"
                                 "  Expected value(s): %s\n"
                                 "  Observed value(s): %r\n")
                                 % (_path_to_str(path), self.description, value))

    def meets_spec(self, _value):
        """Return True if value meets the specification, False otherwise."""
        raise NotImplementedError("'meets_spec' must be implemented in subclasses")




################################################################################
################################################################################
## Tests for basic types

class IsInt(MakefileSpec):
    """Require that the value is either an Int or a Long."""

    def __init__(self, description = "an integer", default = DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, (types.IntType, types.LongType)) \
          and not isinstance(value, types.BooleanType)


class IsUnsignedInt(IsInt):
    """Require that the value is either an Int or a Long, and >= 0."""

    def __init__(self, description = "an unsigned integer", default = DEFAULT_NOT_SET):
        IsInt.__init__(self, description, default)

    def meets_spec(self, value):
        return IsInt.meets_spec(self, value) & (value >= 0)


class IsFloat(MakefileSpec):
    """Require that the value is a float (does not cover integer types)."""

    def __init__(self, description = "a float", default = DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, types.FloatType)


class IsBoolean(MakefileSpec):
    """Require that the value is a boolean (True/False)."""

    def __init__(self, description = "a boolean", default = DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, types.BooleanType)


class IsStr(MakefileSpec):
    """Require that the value is a non-empty string."""

    def __init__(self, description = "a non-empty string", default = DEFAULT_NOT_SET):
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return isinstance(value, types.StringTypes) and value


class IsNone(MakefileSpec):
    """Require that the value is None, typically signifying that
    the value was not set in the makefile."""

    def __init__(self, description = "None or not set", default = DEFAULT_NOT_SET):
        if default is not DEFAULT_NOT_SET:
            raise NotImplementedError("IsNone does not implement a default value")
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        return value is None




################################################################################
################################################################################
## BinaryOperators

class _BinaryOperator(MakefileSpec):
    def __init__(self, description, default, opfunc, rvalue, key = None):
        self._operator = opfunc
        self._keyfunc  = key
        self._rvalue   = rvalue
        description = description.format(rvalue = repr(rvalue))
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        if self._keyfunc is not None:
            value = self._keyfunc(value)
        return self._operator(value, self._rvalue)


class _BinarySetOperator(_BinaryOperator):
    def __init__(self, description, default, opfunc, rvalue, key = None):
        _BinaryOperator.__init__(self, description, default, opfunc, rvalue, key)

    def meets_spec(self, value):
        if not isinstance(value, (types.ListType,) + types.StringTypes):
            return False
        return _BinaryOperator.meets_spec(self, value)



class ValueLT(_BinaryOperator):
    """Less-than operator. If a key-function is set, this function is applied to
    the value before checking the specification (key(value) < rvalue)."""

    def __init__(self, rvalue, key = None, description = "value < {rvalue}", default = DEFAULT_NOT_SET):
        _BinaryOperator.__init__(self, description, default, operator.lt, rvalue, key)


class ValueLE(_BinaryOperator):
    """Less-than-or-equal-to operator. If a key-function is set, this function is
    applied to the value before checking the specification (key(value) <= rvalue)."""

    def __init__(self, rvalue, key = None, description = "value <= {rvalue}", default = DEFAULT_NOT_SET):
        _BinaryOperator.__init__(self, description, default, operator.le, rvalue, key)


class ValueGE(_BinaryOperator):
    """Greater-than-or-equal-to operator. If a key-function is set, this function is
    applied to the value before checking the specification (key(value) >= rvalue)."""

    def __init__(self, rvalue, key = None, description = "value >= {rvalue}", default = DEFAULT_NOT_SET):
        _BinaryOperator.__init__(self, description, default, operator.ge, rvalue, key)


class ValueGT(_BinaryOperator):
    """Greater-than operator. If a key-function is set, this function is applied
    to the value before checking the specification (key(value) > rvalue)."""

    def __init__(self, rvalue, key = None, description = "value > {rvalue}", default = DEFAULT_NOT_SET):
        _BinaryOperator.__init__(self, description, default, operator.gt, rvalue, key)


class ValueIn(_BinaryOperator):
    """Require that values are found in a set of values. Strings are compared in a
    case-sensitive manner. For case-insensitive comparisons, see 'StringIn'."""

    def __init__(self, rvalue, key = None, description = "value in {rvalue}", default = DEFAULT_NOT_SET):
        description = description.format(rvalue = _list_values(rvalue, "or"))
        in_operator = lambda lvalue, rvalue: (lvalue in rvalue)
        _BinaryOperator.__init__(self, description, default, in_operator, rvalue, key)


class ValuesIntersect(_BinarySetOperator):
    """Require that a set of values overlap with a pre-defined set of values (as set in
    the constructor). For strings, values are compared in a non case-sensitive manner.
    For case-insensitive comparisons, see 'StringsIntersect'."""

    def __init__(self, rvalue, key = None, description = "contains {rvalue}", default = DEFAULT_NOT_SET):
        description        = description.format(rvalue = _list_values(rvalue, "and/or"))
        intersect_operator = lambda lvalue, rvalue: bool(frozenset(lvalue).intersection(rvalue))
        _BinarySetOperator.__init__(self, description, default, intersect_operator, frozenset(rvalue), key)


class ValuesSubsetOf(_BinarySetOperator):
    """Require that a set of values are a subset of a pre-defined set of values (as set in
    the constructor). For strings, values are compared in a case-sensitive manner. For
    case-insensitive comparisons, see 'StringsSubsetOf'.

    Note that empty sets are always considered to be a subset of the pre-defined set."""

    def __init__(self, rvalue, key = None, description = "subset of {rvalue}", default = DEFAULT_NOT_SET):
        description     = description.format(rvalue = _list_values(rvalue, "and"))
        subset_operator = lambda lvalue, rvalue: bool(frozenset(lvalue).issubset(rvalue))
        _BinarySetOperator.__init__(self, description, default, subset_operator, frozenset(rvalue), key)


class ValueMissing(MakefileSpec):
    """Used to signify empty substructures in the makefile specification."""

    def __init__(self, description = "no values"):
        MakefileSpec.__init__(self, description, DEFAULT_NOT_SET)

    def meets_spec(self, _value):
        return False


################################################################################
################################################################################
## Logical operators

class _MultipleSpecs(MakefileSpec): # pylint: disable=W0223
    def __init__(self, specs, kwargs, name, prefix = "", postfix = "", join_by = " ", fmt = "%s"):
        self._specs = _instantiate_specs(specs)
        if not self._specs:
            raise ValueError("No specification given to %r" % (name.title(),))
        elif not all((spec.default is DEFAULT_NOT_SET) for spec in self._specs):
            raise ValueError("Default values cannot be set in specs given to logical operators")

        description = []
        for spec in self._specs:
            description.append(fmt % (spec.description,))
        description = "%s%s%s" % (prefix, join_by.join(description), postfix)
        MakefileSpec.__init__(self, description, kwargs.get("default", DEFAULT_NOT_SET))


class And(_MultipleSpecs):
    """Takes one or more specification objects, and requires that values meets all
    of these specifications. A default value may be set for the 'And' specification,
    but not for the specifications given to the 'And' object."""

    def __init__(self, *specs, **kwargs):
        _MultipleSpecs.__init__(self, specs, kwargs, "And", join_by = " and ", fmt = "(%s)")

    def meets_spec(self, value):
        return all(spec.meets_spec(value) for spec in self._specs)


class Or(_MultipleSpecs):
    """Takes one or more specification objects, and requires that values meets at least
    one these specifications. A default value may be set for the 'Or' specification,
    but not for the specifications given to the 'Or' object."""

    def __init__(self, *specs, **kwargs):
        _MultipleSpecs.__init__(self, specs, kwargs, "Or", join_by = " or ", fmt = "(%s)")

    def meets_spec(self, value):
        return any(spec.meets_spec(value) for spec in self._specs)


class Xor(_MultipleSpecs):
    """Takes two specification objects, and requires that values meets ONE and ONLY ONE
    of these specifications. A default value may be set for the 'Xor' specification,
    but not for the specifications given to the 'Xor' object."""

    def __init__(self, *specs, **kwargs):
        if len(specs) != 2:
            raise ValueError("'Xor' takes exactly 2 specification, not %i" % (len(specs),))

        _MultipleSpecs.__init__(self, specs, kwargs, "Xor", join_by = " xor ", fmt = "(%s)")

    def meets_spec(self, value):
        return sum(spec.meets_spec(value) for spec in self._specs) == 1


class Not(_MultipleSpecs):
    """Takes a single specification object, and requires that values do NOT meet this
    specification. A default value may be set for the 'Not' specification, but not for
    the specifications given to the 'Not' object."""

    def __init__(self, spec, **kwargs):
        _MultipleSpecs.__init__(self, [spec], kwargs, "Not", prefix = "not ", fmt = "(%s)")

    def meets_spec(self, value):
        return not self._specs[0].meets_spec(value)




################################################################################
################################################################################
## String operators
##
## In addition to providing string-specific operators (is uppercase, ends/starts
## with), "in" and set operators are provided which do case-insensitive
## comparsions. For case-sensitive operations, use the Value* specifications.


class StringIn(_BinaryOperator):
    """Require that values are found in a set of values. For strings, the comparison
    is done in a case-insensitive. For case-sensitive comparisons, see 'ValueIn'."""

    def __init__(self, rvalues, key = None, description = "one of {rvalue}, case-insentive", default = DEFAULT_NOT_SET):
        description = description.format(rvalue = _list_values(rvalues, "or"))
        rvalues     = frozenset(map(_safe_coerce_to_lowercase, rvalues))
        opfunc      = lambda lvalue, rvalue: (_safe_coerce_to_lowercase(lvalue) in rvalue)
        _BinaryOperator.__init__(self, description, default, opfunc, rvalues)


class _StrSetOperator(_BinaryOperator):
    def __init__(self, description, default, opfunc, rvalues, key = None):
        rvalues     = frozenset(map(_safe_coerce_to_lowercase, rvalues))
        _BinaryOperator.__init__(self, description, default, opfunc, rvalues, key)

    def meets_spec(self, value):
        if not isinstance(value, (types.ListType,) + types.StringTypes):
            return False

        lvalues = frozenset(map(_safe_coerce_to_lowercase, value))
        return _BinaryOperator.meets_spec(self, lvalues)


class StringsIntersect(_StrSetOperator):
    """Require that a set of values overlap with a pre-defined set of values (as set in
    the constructor). For strings, values are compared in a case-insensitive manner.
    For case-sensitive comparisons, see 'ValuesIntersect'."""

    def __init__(self, rvalue, key = None,
                 description = "contains {rvalue}, case-insentive",
                 default = DEFAULT_NOT_SET):
        description        = description.format(rvalue = _list_values(rvalue, "and/or"))
        intersect_operator = lambda lvalue, rvalue: bool(frozenset(lvalue).intersection(rvalue))
        _StrSetOperator.__init__(self, description, default, intersect_operator, frozenset(rvalue), key)


class StringsSubsetOf(_StrSetOperator):
    """Require that a set of values are a subset of a pre-defined set of values (as set in
    the constructor). For strings, values are compared in a case-insensitive manner. For
    case-sensitive comparisons, see 'ValuesSubsetOf'.

    Note that empty sets are always considered to be a subset of the pre-defined set."""

    def __init__(self, rvalue, key = None,
                 description = "subset of {rvalue}, case-insentive",
                 default = DEFAULT_NOT_SET):
        description     = description.format(rvalue = _list_values(rvalue, "and"))
        subset_operator = lambda lvalue, rvalue: bool(frozenset(lvalue).issubset(rvalue))
        _StrSetOperator.__init__(self, description, default, subset_operator, frozenset(rvalue), key)


class StringIsUppercase(IsStr):
    """Require that the value is a uppercase, non-empty string."""

    def __init__(self, default = DEFAULT_NOT_SET):
        IsStr.__init__(self, "an uppercase non-empty string", default)

    def meets_spec(self, value):
        return IsStr().meets_spec(value) and value.isupper()


class StringStartsWith(IsStr):
    """Require that the value is a string with given prefix."""

    def __init__(self, prefix, default = DEFAULT_NOT_SET):
        assert prefix and isinstance(prefix, types.StringTypes)
        self._prefix = prefix
        IsStr.__init__(self, "a string with the prefix %r" % (prefix,), default)

    def meets_spec(self, value):
        return IsStr.meets_spec(self, value) \
          and value.startswith(self._prefix)


class StringEndsWith(IsStr):
    """Require that the value is a string with given postfix."""

    def __init__(self, postfix, default = DEFAULT_NOT_SET):
        assert postfix and isinstance(postfix, types.StringTypes)
        self._postfix = postfix
        IsStr.__init__(self, "a string with the postfix %r" % (postfix,), default)

    def meets_spec(self, value):
        return IsStr.meets_spec(self, value) \
          and value.endswith(self._postfix)




################################################################################
################################################################################
## Tests for collections

class IsListOf(_MultipleSpecs):
    """Require that the value is a list, the contents of which matches one or
    more of the provided specifications."""

    def __init__(self, *specs, **kwargs):
        _MultipleSpecs.__init__(self, specs, kwargs, "IsListOf",
                                prefix = "[", postfix = ", ...]",
                                join_by = " or ", fmt = "(%s)")

    def meets_spec(self, value):
        if not isinstance(value, types.ListType):
            return False
        return all(any(spec.meets_spec(lstvalue) for spec in self._specs) for lstvalue in value)


class IsDictOf(MakefileSpec):
    """Require that the value is a list, the keys/values of which matches
    the specifications provided for keys/values."""

    def __init__(self, key_spec, value_spec, default = DEFAULT_NOT_SET):
        self._key_spec   = _instantiate_spec(key_spec)
        self._value_spec = _instantiate_spec(value_spec)
        if self._key_spec.default is not DEFAULT_NOT_SET:
            raise ValueError("Default values cannot be set in key-specs for IsDictOf")
        elif self._value_spec.default is not DEFAULT_NOT_SET:
            raise ValueError("Default values cannot be set in value-specs for IsDictOf")

        description = "{(%s) : (%s)}" \
          % (self._key_spec.description, self._value_spec.description)
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value):
        if not isinstance(value, types.DictType):
            return False

        for (key, value) in value.iteritems():
            if not (self._key_spec.meets_spec(key) and self._value_spec.meets_spec(value)):
                return False

        return True




################################################################################
################################################################################
## Helper functions

def _is_spec(spec):
    if isinstance(spec, MakefileSpec):
        return True
    elif isinstance(spec, types.ClassType) and issubclass(spec, MakefileSpec):
        return True
    return False


def _instantiate_spec(spec):
    if isinstance(spec, MakefileSpec):
        return spec
    elif isinstance(spec, types.ClassType) and issubclass(spec, MakefileSpec):
        return spec()
    else:
        raise TypeError("Specifications must derive from 'MakefileSpec'")


def _instantiate_specs(specs):
    return map(_instantiate_spec, specs)


def _safe_coerce_to_lowercase(value):
    if isinstance(value, types.StringTypes):
        return value.lower()
    return value


def _list_values(values, sep):
    values = map(repr, values)
    if len(values) > 2:
        values = (", ".join(values[:-1]) + ",", values[-1])
    if len(values) == 2:
        values = (" ".join((values[0], sep, values[1])),)

    return values[0]


def _get_summary_spec(specs_or_keys):
    specs, keys = group_by_pred(_is_spec, specs_or_keys)
    if specs and keys:
        return Or(ValueIn(keys, description = "key in {rvalue}"), *specs)
    elif specs:
        return Or(*specs)
    elif keys:
        return ValueIn(keys, description = "key in {rvalue}")
    return ValueMissing()


def _get_matching_spec_or_value(value, specs, path):
    if value in specs:
        return value

    for spec in specs:
        if _is_spec(spec):
            try:
                _instantiate_spec(spec)(path, value)
                return spec
            except MakefileError:
                pass
    # No matching key or spec; create combined spec for error message
    _get_summary_spec(specs)(path, value)
    assert False # pragma: no coverage


def _set_default_values(data, specification):
    for cur_key in specification:
        if (not _is_spec(cur_key)) and (cur_key not in data):
            default_value = specification[cur_key]
            if _is_spec(default_value):
                default_value = _instantiate_spec(default_value)
                if (default_value.default is DEFAULT_NOT_SET):
                    continue
                default_value = default_value.default

            if not isinstance(default_value, WithoutDefaults):
                # Setting of values in the dict will be accomplished
                # in subsequent calls to _set_default_values
                if isinstance(default_value, dict):
                    default_value = {}

                # Prevent clobbering of values when re-using sub-specs
                data[cur_key] = copy.deepcopy(default_value)


def _path_to_str(path):
    return ":".join(str(field) for field in path)


CLI_PARAMETERS = Or(IsListOf(IsStr, IsInt, IsFloat),
                    Or(IsStr, IsInt, IsFloat, IsNone))
