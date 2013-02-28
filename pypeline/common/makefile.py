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
import yaml

import os
import copy
import types
import hashlib
import datetime

from pypeline.common.utilities import safe_coerce_to_tuple


class MakefileError(RuntimeError):
    pass


def read_makefile(filename, defaults, validation):
    try:
        with open(filename) as makefile:
            string = makefile.read()
            data = yaml.safe_load(string)
    except Exception, e:
        raise MakefileError(e)

    final = copy.deepcopy(defaults)
    final = apply_defaults(data, final, ("root",))
    validate_makefile(final, validation, ("root",))

    mtime = os.path.getmtime(os.path.realpath(filename))
    return {"Makefile"   : final,
            "Statistics" : {"Filename" : filename,
                          "Hash"     : hashlib.sha1(string).hexdigest(),
                          "MTime"    : datetime.datetime.fromtimestamp(mtime).strftime("%F %T.%f %z")}}


def apply_defaults(current, defaults, path):
    if not isinstance(current, dict) or not isinstance(defaults, dict):
        raise MakefileError("Invalid structure in makefile at path %s" % ":".join(path))

    for key in current:
        if isinstance(current[key], dict) and isinstance(defaults.get(key), dict):
            apply_defaults(current[key], defaults[key], path + (key,))
        else:
            defaults[key] = current[key]

    return defaults



def validate_makefile(current, reference, path):
    if callable(reference):
        reference(path, current)
    elif isinstance(current, dict) and isinstance(reference, dict):
        for cur_key in current:
            if cur_key not in reference:
                for ref_key in reference:
                    if callable(ref_key):
                        ref_key(path, cur_key)
                        validate_makefile(current[cur_key], reference[ref_key], path + (cur_key,))
                        break
                else:
                    raise MakefileError("Unknown key in makefile at '%s': %s" \
                                    % (":".join(path), cur_key))
            else:
                validate_makefile(current[cur_key], reference[cur_key], path + (cur_key,))
    else:
        raise MakefileError("Inconsistency between makefile format and current makefile at %s\nExpected dict, found %s!" \
                            % (":".join(path), current.__class__.__name__))


def And(*funcs):
    def _And(path, value):
        for func in funcs:
            func(path, value)
    return _And


def Or(*funcs):
    names = [func.__name__ for func in funcs]
    name = "(%s)" % (" or ".join(names))

    def _Or(path, value):
        error = RuntimeError("No function given to Or() at '%s'!" % ":".join(path))
        for func in funcs:
            try:
                func(path, value)
                return
            except MakefileError, error:
                pass

        raise MakefileError("Value %s at '%s' in makefile does not meet requirements: %s" \
                            % (repr(value), ":".join(path), name))
    _Or.__name__ = name
    return _Or


def OneOf(*args, **kwargs):
    assert "case_sensitive" not in kwargs or len(kwargs) == 1, "No keywords but 'case_sensitive' are allowed: " + repr(kwargs)
    if kwargs.get("case_sensitive", True):
        def _OneOf_CaseSensitive(path, value):
            if value not in args:
                raise MakefileError("Value for '%s' must be one of '%s', not '%s'!" \
            % (":".join(path), "', '".join(map(repr, args)), repr(value)))
        return _OneOf_CaseSensitive

    args = map(_safe_coerce_to_lowercase, args)
    def _OneOf_CaseInsensitive(path, value):
        if _safe_coerce_to_lowercase(value) not in args:
            raise MakefileError("Value for '%s' must be one of %s (case insensitive), not %s!" \
                                % (":".join(path), ", ".join(map(repr, args)), repr(value)))
    return _OneOf_CaseInsensitive


def AnyOf(*args, **kwargs):
    keywords = set(("case_sensitive", "min_items", "max_items"))
    if set(kwargs) - keywords:
        raise RuntimeError("Invalid keyword(s) given to AnyOf: %s" \
                           ", ".join(map(repr, set(kwargs) - keywords)))

    key_func = _this
    if not kwargs.get("case_sensitive", True):
        args = map(_safe_coerce_to_lowercase, args)
        key_func = _safe_coerce_to_lowercase

    min_items = int(kwargs.get("min_items", 0))
    max_items = int(kwargs.get("max_items", len(args)))

    def _AnyOf(path, value):
        hits, values = 0, safe_coerce_to_tuple(value)
        for value in values:
            if key_func(value) not in args:
                raise MakefileError("Value for '%s' must be among %s, not %s!" \
                                    % (":".join(path), ", ".join(map(repr, args)), repr(value)))

        if not min_items <= len(values) <= max_items:
            raise MakefileError("Expected %s to %s values, found %i!" \
                                % (repr(min_items), repr(max_items), len(values)))
    return _AnyOf



def IsListOf(*args):
    test_func = Or(*args)
    def _IsListOf(path, value):
        if not isinstance(value, list):
            raise MakefileError("Value for '%s' must be a list, not '%s'!" \
                                % (":".join(path), repr(value)))

            for item in value:
                test_func(item)
    return _IsListOf


def IsInt(path, value):
    if isinstance(value, bool) or not isinstance(value, int):
        raise MakefileError("Value for '%s' must be an integer, not '%s'!" \
                            % (":".join(path), repr(value)))

def IsFloat(path, value):
    if not isinstance(value, float):
        raise MakefileError("Value for '%s' must be an float, not '%s'!" \
                            % (":".join(path), repr(value)))


def IsInRange(min, max):
    assert min < max
    def _IsInRange(path, value):
        if not min <= value < max:
            raise MakefileError("Value for '%s' must be in range %s <= value < %s, not '%s'!" \
                                % (":".join(path), min, max, value))
    return _IsInRange


def IsBoolean(path, value):
    if not isinstance(value, bool):
        raise MakefileError("Value for '%s' must be a boolean, not %s!" \
                            % (":".join(path), repr(value)))


def IsStr(path, value):
    if not isinstance(value, types.StringTypes):
        raise MakefileError("Value for '%s' must be a string, not '%s'!" \
                            % (":".join(path), repr(value)))

def IsNone(path, value):
    if value is not None:
        raise MakefileError("Value for '%s' must be unspecified/None, not '%s'!" \
                            % (":".join(path), repr(value)))


def _safe_coerce_to_lowercase(value):
    if isinstance(value, types.StringTypes):
        return value.lower()
    return value


def _this(value):
    return value
