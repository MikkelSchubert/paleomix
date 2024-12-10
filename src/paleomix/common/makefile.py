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
for convenience.


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
value has accidentally been left blank ('IsStr' requires a NON-EMPTY string):
-------------------------------------------------------------------------------
|  Makefile requirement not met at ...:
|    Expected value: a non-empty string
|    Observed value: ''
-------------------------------------------------------------------------------
"""

from __future__ import annotations

import copy
import logging
import re
from collections.abc import Collection
from typing import Dict, Hashable, List, Sequence, Tuple, Type, Union

from typing_extensions import TypeGuard

from paleomix.common import yaml
from paleomix.common.fileutils import PathTypes, fspath

BasicType = Union[str, int, float, bool, None]
SpecType = Union["_SpecBase", "Type[_SpecBase]"]
SpecTree = Union[
    BasicType,
    SpecType,
    Dict[Union[str, SpecType], "SpecTree"],
    List["SpecTree"],
]
SpecPath = Tuple[BasicType, ...]


class MakefileError(RuntimeError):
    """Raised if a makefile is unreadable, or does not meet specifications."""


def read_makefile(filename: PathTypes, specification: SpecTree) -> object:
    """Reads and parses a makefile using the given specification."""
    try:
        with open(fspath(filename)) as handle:
            data = yaml.safe_load(handle)
    except (yaml.YAMLError, yaml.DuplicateKeyError) as error:
        raise MakefileError(error) from None

    return process_makefile(data, specification)


def process_makefile(
    data: object,
    specification: SpecTree,
    path: SpecPath = (),
    apply_defaults: bool = True,
) -> object:
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
        data = _instantiate_spec(specification)(path, data)
    elif isinstance(data, (dict, type(None))) and isinstance(specification, dict):
        # A limitation of YAML is that empty subtrees are equal to None;
        # this check ensures that empty subtrees to be handled properly
        if data is None:
            data = {}

        _process_default_values(data, specification, path, apply_defaults)

        for cur_key in data:
            ref_key = _get_matching_spec_or_value(
                cur_key, specification, (*path, cur_key)
            )
            data[cur_key] = process_makefile(
                data[cur_key], specification[ref_key], (*path, cur_key), apply_defaults
            )
    elif isinstance(data, (list, type(None))) and isinstance(specification, list):
        if not _is_spec_list(specification):
            raise TypeError(
                "Lists contains non-specification objects ({!r}): {!r}".format(
                    _path_to_str(path), specification
                )
            )
        elif data is None:  # See comment above
            data = []

        specification = IsListOf(*specification)
        _instantiate_spec(specification)(path, data)
    elif not isinstance(specification, (dict, list)):
        raise TypeError(
            "Unexpected type in makefile specification at {!r}: {!r}!".format(
                _path_to_str(path), specification
            )
        )
    else:
        raise MakefileError(
            "Inconsistency between makefile specification and "
            "current makefile at {}:\n    Expected {}, "
            "found {} {!r}!".format(
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


class _SpecBase:
    pass


class WithoutDefaults(_SpecBase):
    """Wrapper object, that tells 'process_makefile' not to apply
    default values for the wrapped specification. See module docs
    for example usage.
    """

    def __init__(self, specification: SpecTree) -> None:
        self.specification = specification


class PreProcessMakefile(_SpecBase):
    """Allows pre-processing of a part of a makefile prior to validation; when
    encountered, the object is called with the current value, and is expected
    to return a tuple containing (value, specification), which are then used
    subsequently. This allows transformation of fields for backwards
    compatibility.
    """

    def __call__(self, path: SpecPath, value: object) -> tuple[object, SpecTree]:
        """Must return (value, specification) tuple."""
        raise NotImplementedError  # pragma: no coverage


class MakefileSpec(_SpecBase):
    """Base-class for specifications, from which ALL specification
    objects are expected to derive. Sub-classes must implement the
    'meets_spec' function, which must return True or False depending
    on whether or not the given value meets the specification.
    """

    def __init__(
        self,
        description: str = "N/A",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        """description -- A string describing the specification.
        default     -- A default value, or DEFAULT_NOT_SET if not used. If a
                       value is set, it is copied before being applied."""

        self.description = description
        self.default = default
        if (default not in (DEFAULT_NOT_SET, REQUIRED_VALUE)) and not self.meets_spec(
            default
        ):
            raise ValueError(
                "Default value does not meet requirements:\n"
                f"  Expected value: {description}\n"
                f"  Observed value: {default!r}\n"
            )

    def __call__(self, path: SpecPath, value: object) -> object:
        if not self.meets_spec(value):
            raise MakefileError(
                f"Makefile requirement not met at {_path_to_str(path)!r}:\n"
                f"  Expected value: {self.description}\n"
                f"  Observed value: {value!r}\n"
                f"  Observed type:  {type(value).__name__}"
            )

        return value

    def meets_spec(self, _value: object) -> bool:
        """Return True if value meets the specification, False otherwise."""
        raise NotImplementedError


###############################################################################
###############################################################################
# Tests for basic types


class IsAny(MakefileSpec):
    """Any value is allowed; if validation is too complex to implement using specs"""

    def __init__(
        self,
        description: str = "any value",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, _value: object) -> bool:
        return True


class IsInt(MakefileSpec):
    """Require that the value is either an Int or a Long."""

    def __init__(
        self,
        description: str = "an integer",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value: object) -> TypeGuard[int]:
        return isinstance(value, int) and not isinstance(value, bool)


class IsUnsignedInt(IsInt):
    """Require that the value is either an Int or a Long, and >= 0."""

    def __init__(
        self,
        description: str = "an unsigned integer",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        IsInt.__init__(self, description, default)

    def meets_spec(self, value: object) -> bool:
        return super().meets_spec(value) and value >= 0


class IsFloat(MakefileSpec):
    """Require that the value is a float (does not cover integer types)."""

    def __init__(
        self,
        description: str = "a float",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value: bool) -> TypeGuard[float]:
        return isinstance(value, float)


class IsBoolean(MakefileSpec):
    """Require that the value is a boolean (True/False)."""

    def __init__(
        self,
        description: str = "a boolean",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value: object) -> TypeGuard[bool]:
        return isinstance(value, bool)


class IsStr(MakefileSpec):
    """Require that the value is a non-empty string."""

    def __init__(
        self,
        description: str | None = None,
        default: object = DEFAULT_NOT_SET,
        min_len: int = 1,
    ) -> None:
        if description is None:
            if min_len == 0:
                description = "a string"
            elif min_len == 1:
                description = "a non-empty string"
            elif min_len >= 2:
                description = f"a string at least {min_len} characters long"
            else:
                raise ValueError("min_len must be non-negative")

        self._min_len = min_len
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value: object) -> TypeGuard[str]:
        return isinstance(value, str) and len(value) >= self._min_len


class IsNone(MakefileSpec):
    """Require that the value is None, typically signifying that
    the value was not set in the makefile."""

    def __init__(
        self,
        description: str = "null or not set",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        if default is not DEFAULT_NOT_SET:
            raise NotImplementedError("IsNone does not support default values")
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value: object) -> TypeGuard[None]:
        return value is None


class ValueMissing(MakefileSpec):
    """Used to signify empty substructures in the makefile specification."""

    def __init__(self, description: str = "no values") -> None:
        MakefileSpec.__init__(self, description, DEFAULT_NOT_SET)

    def meets_spec(self, _value: object) -> bool:
        return False


class DeprecatedOption(MakefileSpec):
    """Used to signify substructures that will eventually be removed."""

    def __init__(self, spec: MakefileSpec) -> None:
        self._spec = spec
        if not isinstance(spec, MakefileSpec):
            raise ValueError(spec)

        MakefileSpec.__init__(self, spec.description, spec.default)

    def __call__(self, path: SpecPath, value: object) -> object:
        self._spec(path, value)

        log = logging.getLogger(__name__)
        log.warning(
            "option has been deprecated and will be removed in the future: %s",
            _path_to_str(path),
        )

        return value

    def meets_spec(self, value: object) -> bool:
        return self._spec.meets_spec(value)


class RemovedOption(MakefileSpec):
    """Used to signify substructures that have been removed, and are hence ignored."""

    def __init__(self, description: str = "removed settings") -> None:
        MakefileSpec.__init__(self, description, DEFAULT_NOT_SET)

    def __call__(self, path: SpecPath, value: object) -> object:
        log = logging.getLogger(__name__)
        log.warning(
            "option has been removed and no longer has any effect: %s",
            _path_to_str(path),
        )

        return value

    def meets_spec(self, _value: object) -> bool:
        return True


###############################################################################
###############################################################################
# BinaryOperators


class ValueIn(MakefileSpec):
    def __init__(
        self,
        rvalues: Sequence[Hashable],
        description: str = "value in {rvalue}",
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        self._rvalues = tuple(rvalues)
        self._normalized_values: dict[str, str] = {}
        for value in self._rvalues:
            if isinstance(value, str):
                assert value.lower() not in self._normalized_values, value
                self._normalized_values[value.lower()] = value

        if isinstance(default, str):
            default = self._normalized_values.get(default, default)

        MakefileSpec.__init__(
            self,
            description=description.format(rvalue=_list_values(self._rvalues, "or")),
            default=default,
        )

    def meets_spec(self, value: object) -> bool:
        return _is_hashable(value) and value in self._rvalues

    def __call__(self, path: SpecPath, value: object) -> object:
        if isinstance(value, str):
            value = self._normalized_values.get(value.lower(), value)

        return super().__call__(path, value)


class ValuesIntersect(MakefileSpec):
    def __init__(
        self,
        rvalues: Sequence[Hashable],
        description: str | None = None,
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        self._rvalues = rvalues
        if not description:
            description = "one or more of {}".format(_list_values(rvalues, "and"))

        MakefileSpec.__init__(self, description=description, default=default)

    def meets_spec(self, lvalue: object) -> bool:
        try:
            return not isinstance(lvalue, dict) and bool(
                frozenset(lvalue).intersection(self._rvalues)
            )
        except TypeError:
            return False


class ValuesSubsetOf(MakefileSpec):
    def __init__(
        self,
        rvalues: Sequence[Hashable],
        description: str | None = None,
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        self._rvalues = rvalues
        description = description or "subset of {}".format(_list_values(rvalues, "and"))

        MakefileSpec.__init__(self, description=description, default=default)

    def meets_spec(self, lvalue: object) -> bool:
        try:
            return not isinstance(lvalue, dict) and frozenset(lvalue).issubset(
                self._rvalues
            )
        except TypeError:
            return False


###############################################################################
###############################################################################
# Logical operators


class _MultipleSpecs(MakefileSpec):
    """Base-class for logical operators for one or more specifications."""

    def __init__(
        self,
        specs: Sequence[SpecType],
        default: object,
        name: str,
        prefix: str = "",
        postfix: str = "",
        join_by: str = " ",
        fmt: str = "%s",
    ) -> None:
        self._specs = [_instantiate_spec(spec) for spec in specs]
        if not self._specs:
            raise ValueError(f"No specification given to {name.title()!r}")
        elif not all((spc.default is DEFAULT_NOT_SET) for spc in self._specs):
            raise ValueError(
                "Default values cannot be set in specs given to logical operators"
            )

        description = [(fmt % (spec.description,)) for spec in self._specs]
        description = f"{prefix}{join_by.join(description)}{postfix}"
        MakefileSpec.__init__(self, description, default)


class And(_MultipleSpecs):
    """Takes one or more specification objects, and requires that values meets
    all of these specifications. A default value may be set for the 'And'
    specification, but not for the specifications given to the 'And' object.
    """

    def __init__(self, *specs: SpecType, default: object = DEFAULT_NOT_SET) -> None:
        _MultipleSpecs.__init__(
            self, specs, default, "And", join_by=" and ", fmt="(%s)"
        )

    def meets_spec(self, value: object) -> bool:
        return all(spec.meets_spec(value) for spec in self._specs)


class Or(_MultipleSpecs):
    """Takes one or more specification objects, and requires that values meets
    at least one these specifications. A default value may be set for the 'Or'
    specification, but not for the specifications given to the 'Or' object.
    """

    def __init__(self, *specs: SpecType, default: object = DEFAULT_NOT_SET) -> None:
        _MultipleSpecs.__init__(self, specs, default, "Or", join_by=" or ", fmt="(%s)")

    def meets_spec(self, value: object) -> bool:
        return any(spec.meets_spec(value) for spec in self._specs)


class Not(_MultipleSpecs):
    """Takes a single specification object, and requires that values do NOT
    meet this specification. A default value may be set for the 'Not'
    specification, but not for the specifications given to the 'Not' object.
    """

    def __init__(self, spec: SpecType, default: object = DEFAULT_NOT_SET) -> None:
        _MultipleSpecs.__init__(self, [spec], default, "Not", prefix="not ", fmt="(%s)")

    def meets_spec(self, value: object) -> bool:
        return not self._specs[0].meets_spec(value)


###############################################################################
###############################################################################
# String operators
#


class StringStartsWith(IsStr):
    """Require that the value is a string with given prefix."""

    def __init__(self, prefix: str, default: object = DEFAULT_NOT_SET) -> None:
        self._prefix = prefix
        IsStr.__init__(self, f"a string with prefix {prefix!r}", default)

    def meets_spec(self, value: object) -> TypeGuard[str]:
        return super().meets_spec(value) and value.startswith(self._prefix)


class StringEndsWith(IsStr):
    """Require that the value is a string with given postfix."""

    def __init__(self, postfix: str, default: object = DEFAULT_NOT_SET) -> None:
        self._postfix = postfix
        IsStr.__init__(self, f"a string with postfix {postfix!r}", default)

    def meets_spec(self, value: object) -> TypeGuard[str]:
        return super().meets_spec(value) and value.endswith(self._postfix)


class FASTQPath(IsStr):
    """String path with optional, case-insensitive {pair} key to specify where the
    mate 1/2 identifier is located."""

    _PAIR_KEY = re.compile("{pair}", re.IGNORECASE)

    def __init__(
        self,
        paired_end: bool | None = False,
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        self._paired_end = paired_end
        if paired_end is None:
            desc = "a path with an optional {pair} key"
        elif paired_end:
            desc = "a path with a {pair} key"
        else:
            desc = "a path without a {pair} key"

        IsStr.__init__(self, desc, default)

    def meets_spec(self, value: object) -> TypeGuard[str]:
        return super().meets_spec(value) and (
            self._paired_end is None or self.is_paired_end(value) == self._paired_end
        )

    @classmethod
    def is_paired_end(cls, path: str) -> bool:
        """Returns true if a path contains a {pair} component."""
        return bool(cls._PAIR_KEY.search(path))

    @classmethod
    def format(cls, path: str, pair: int) -> str:
        assert pair in (1, 2), pair

        result, count = cls._PAIR_KEY.subn(str(pair), path)
        assert count

        return result


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

    def __init__(self, *specs: SpecType, default: object = DEFAULT_NOT_SET) -> None:
        _MultipleSpecs.__init__(
            self,
            specs,
            default,
            "IsListOf",
            prefix="[",
            postfix=", ...]",
            join_by=" or ",
            fmt="(%s)",
        )

    def meets_spec(self, value: object) -> bool:
        if not isinstance(value, list):
            return False

        return all(any(spec.meets_spec(it) for spec in self._specs) for it in value)


class IsDictOf(MakefileSpec):
    """Require that the value is a list, the keys/values of which matches
    the specifications provided for keys/values; if no default value (ie. a
    dictionary) is required, then using the following syntax is preferred:
      {IsType1: IsType2}

    This is equivalent to the following:
      IsDictOf(IsType1, IsType2)
    but also allows multiple type-pairs to be specified.
    """

    def __init__(
        self,
        key_spec: SpecType,
        value_spec: SpecType,
        default: object = DEFAULT_NOT_SET,
    ) -> None:
        self._key_spec = _instantiate_spec(key_spec)
        self._value_spec = _instantiate_spec(value_spec)
        if self._key_spec.default is not DEFAULT_NOT_SET:
            raise ValueError("Default values cannot be set in key-specs")
        elif self._value_spec.default is not DEFAULT_NOT_SET:
            raise ValueError("Default values cannot be set in value-specs")

        description = "{{({}) : ({})}}".format(
            self._key_spec.description,
            self._value_spec.description,
        )
        MakefileSpec.__init__(self, description, default)

    def meets_spec(self, value: object) -> TypeGuard[dict]:
        if not isinstance(value, dict):
            return False

        return all(
            (self._key_spec.meets_spec(key) and self._value_spec.meets_spec(value))
            for key, value in value.items()
        )


###############################################################################
###############################################################################
# Helper functions


def _is_hashable(value: object) -> bool:
    try:
        hash(value)
    except TypeError:
        return False
    else:
        return True


def _is_spec(spec: object) -> TypeGuard[SpecType]:
    """Returns true if 'spec' is a specification instance or class."""
    return isinstance(spec, MakefileSpec) or (
        isinstance(spec, type) and issubclass(spec, MakefileSpec)
    )


def _is_spec_list(spec: list[object]) -> TypeGuard[list[SpecType]]:
    return all(_is_spec(it) for it in spec)


def _instantiate_spec(spec: SpecTree) -> MakefileSpec:
    """Takes a specification instance or class, and returns an instance."""
    if isinstance(spec, MakefileSpec):
        return spec
    elif isinstance(spec, type) and issubclass(spec, MakefileSpec):
        return spec()
    else:
        raise TypeError("Specifications must derive from 'MakefileSpec'")


def _list_values(raw: Sequence[Hashable], sep: str) -> str:
    """Returns list of values as 'values[0], values[1], ..., sep values[-1]':

    $ _list_values([1, 2, 3], "and")
    "1, 2, and 3"
    """
    *head, tail = list(map(repr, raw))
    if len(head) > 1:
        return "{}, {} {}".format(", ".join(head), sep, tail)
    elif len(head) == 1:
        return "{} {} {}".format(*head, sep, tail)
    else:
        return tail


def _get_summary_spec(specs_or_keys: Collection[SpecType | str]) -> MakefileSpec:
    """Returns a specification object that may be used to describe a set of
    requirements. This is used if a key or value does not match the possible
    specs, thereby describing the set of allowed values.
    """
    specs = [it for it in specs_or_keys if not isinstance(it, str)]
    keys = [it for it in specs_or_keys if isinstance(it, str)]
    if specs and keys:
        return Or(ValueIn(keys, description="key in {rvalue}"), *specs)
    elif specs:
        return Or(*specs)
    elif keys:
        return ValueIn(keys, description="key in {rvalue}")
    return ValueMissing()


def _get_matching_spec_or_value(
    value: object,
    specs: Collection[SpecType | str],
    path: SpecPath,
) -> object:
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

    raise AssertionError("combined spec-test should have failed")  # pragma: no coverage


def _process_default_values(
    data: object,
    specification: SpecTree,
    path: SpecPath,
    apply_defaults: bool,
) -> None:
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
                        "A value MUST be supplied for {!r}".format(
                            _path_to_str((*path, cur_key))
                        )
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
                elif isinstance(default_value, list) and not default_value_from_spec:
                    # Lists of specs defaults to empty lists
                    default_value = []

                # Prevent clobbering of values when re-using sub-specs
                data[cur_key] = copy.deepcopy(default_value)


def _path_to_str(path: SpecPath) -> str:
    """Converts a path (tuple of strings) to a printable string."""
    return " :: ".join(str(field) for field in path)
