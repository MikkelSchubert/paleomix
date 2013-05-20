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
from __future__ import print_function

import os
import sys
import types
import subprocess


def _is_cls(obj, *cls_names):
    return obj.__class__.__name__ in cls_names


def _describe_cls(atomiccmd):
    if _is_cls(atomiccmd, "ParallelCmds"):
        return "Atomic commands executed in parallel"
    elif _is_cls(atomiccmd, "SequentialCmds"):
        return "Atomic commands executed sequentially"
    else:
        return "Atomic command"


def _collect_stats(atomiccmd, stats):
    assert atomiccmd not in stats["id"]

    if _is_cls(atomiccmd, "AtomicCmd"):
        cmd_id = stats["id"][atomiccmd] = len(stats["id"])
        pipe   = atomiccmd._files.get("IN_STDIN")
        if _is_cls(pipe, "AtomicCmd"):
            stats["pipe"][pipe] = atomiccmd
    elif _is_cls(atomiccmd, "ParallelCmds", "SequentialCmds"):
        for subcmd in atomiccmd._commands:
            _collect_stats(subcmd, stats)

    return stats


def _build_status(atomiccmd, stats, indent, lines):
    prefix = " " * indent + "Status  = "
    if atomiccmd._proc:
        if atomiccmd.ready():
            lines.append(prefix + "Exited with return-code %i" \
                         % tuple(atomiccmd.join()))
        else:
            lines.append(prefix + "Running ...")


def _build_stdin(atomiccmd, stats, indent, lines):
    pipe   = atomiccmd._files.get("IN_STDIN")
    prefix = " " * indent + "STDIN   = "
    if pipe and pipe in stats["id"]:
        lines.append("%s<%02i>" % (prefix, stats["id"][pipe],))
    elif isinstance(pipe, types.StringTypes):
        lines.append("%s'%s'" % (prefix, pipe))
    elif pipe:
        lines.append("%s<PIPE>" % (prefix,))


def _build_stdout(atomiccmd, stats, indent, lines):
    prefix = " " * indent + "STDOUT  = "
    if atomiccmd in stats["pipe"]:
        pipe = stats["pipe"].get(atomiccmd)
        lines.append("%s<%02i>" % (prefix, stats["id"][pipe],))
    elif atomiccmd._files.get("OUT_STDOUT"):
        filename = atomiccmd._files.get("OUT_STDOUT")
        if filename is not subprocess.PIPE:
            lines.append("%s'%s'" % (prefix, filename))
        else:
            lines.append("%s<PIPE>" % (prefix,))
    elif "OUT_STDOUT" in atomiccmd._pipes:
        filename = atomiccmd._pipes["OUT_STDOUT"]
        lines.append("%s'%s'" % (prefix, filename))


def _build_stderr(atomiccmd, stats, indent, lines):
    prefix = " " * indent + "STDERR  = "
    if atomiccmd._files.get("OUT_STDERR"):
        filename = atomiccmd._files.get("OUT_STDERR")
        lines.append("%s'%s'" % (prefix, filename,))
    elif "OUT_STDERR" in atomiccmd._pipes:
        filename = atomiccmd._pipes["OUT_STDERR"]
        lines.append("%s'%s'" % (prefix, filename))


def _build_cwd(atomiccmd, stats, indent, lines):
    prefix = " " * indent + "CWD     = "
    if atomiccmd._set_cwd:
        cwd = atomiccmd._temp or "<Temporary folder>"
        lines.append("%s'%s'" % (prefix, cwd,))
    else:
        lines.append("%s'%s'" % (prefix, os.getcwd()))


def _pformat(atomiccmd, stats, indent, lines, include_prefix = True):
    s_prefix     = ""
    if include_prefix:
        s_prefix = " " * indent + "- "
        if _is_cls(atomiccmd, "AtomicCmd"):
            cmd_id = stats["id"][atomiccmd]
            s_prefix += "<%02i> " % (cmd_id,)
    s_prefix_len = len(s_prefix)

    if _is_cls(atomiccmd, "AtomicCmd"):
        temp = None if atomiccmd._set_cwd else atomiccmd._temp
        c_prefix = s_prefix + "Command = "
        for line in _pformat_list(atomiccmd._generate_call(temp)).split("\n"):
            lines.append("%s%s" % (c_prefix, line))
            c_prefix = " " * len(c_prefix)

        _build_status(atomiccmd, stats, s_prefix_len, lines)
        _build_stdin(atomiccmd,  stats, s_prefix_len, lines)
        _build_stdout(atomiccmd, stats, s_prefix_len, lines)
        _build_stderr(atomiccmd, stats, s_prefix_len, lines)
        _build_cwd(atomiccmd,    stats, s_prefix_len, lines)
    elif _is_cls(atomiccmd, "ParallelCmds", "SequentialCmds"):
        lines.append("%s%s:" % (s_prefix, _describe_cls(atomiccmd)))
        for subcmd in atomiccmd._commands:
            _pformat(subcmd, stats, s_prefix_len + 2, lines)
    else:
        assert False


def _pformat_list(lst, width = 80):
    """Return a printable representation of a list, where line-breaks
    are inserted between items to minimize the number of lines with a
    width greater than 'width'. Very long items may cause this maximum
    to be exceeded."""
    result = [[]]
    current_width = 0
    for item in map(repr, lst):
        if current_width + len(item) > width:
            result.append([item])
            current_width = len(item)
        else:
            result[-1].append(item)
            current_width += len(item) + 2

    return "[%s]" % (",\n ".join(", ".join(line) for line in result))


def pformat(atomiccmd):
    """Returns a human readable description of an Atomic Cmd or Atomic Set
    of commands. This is currently equivalent to str(cmd_obj)."""
    lines = []
    stats = _collect_stats(atomiccmd, {"id" : {}, "pipe" : {}})
    _pformat(atomiccmd, stats, 0, lines, False)
    return "<%s>" % "\n".join(lines)


def pprint(atomiccmd, out = sys.stdout):
    """Prints a human readable description of an Atomic Cmd or Atomic Set
    of commands. This is currently equivalent to print(str(cmd_obj), ...)."""
    print(pformat(atomiccmd), file = out)

