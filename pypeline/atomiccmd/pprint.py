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
# pylint: disable=W0212

from __future__ import print_function

import os
import sys
import types
import subprocess


def _is_cls(obj, *cls_names):
    return obj.__class__.__name__ in cls_names

def _get_pipe_name(files, pipe):
    if pipe in files:
        return pipe.split("_")[-1] + " "
    return pipe.split("_")[-1] + "*"

def _get_pipe_file(files, pipe):
    pipe_filename = files.get(pipe)
    if pipe_filename:
        return pipe_filename
    return files.get("TEMP_%s" % (pipe,))

def _describe_cls(atomiccmd):
    if _is_cls(atomiccmd, "ParallelCmds"):
        return "Parallel commands"
    elif _is_cls(atomiccmd, "SequentialCmds"):
        return "Sequential commands"
    assert False # pragma: no coverage


def _collect_stats(atomiccmd, stats):
    assert atomiccmd not in stats["id"]

    if _is_cls(atomiccmd, "AtomicCmd"):
        stats["id"][atomiccmd] = len(stats["id"])
        pipe   = _get_pipe_file(atomiccmd._files, "IN_STDIN")
        if _is_cls(pipe, "AtomicCmd"):
            stats["pipe"][pipe] = atomiccmd
    elif _is_cls(atomiccmd, "ParallelCmds", "SequentialCmds"):
        for subcmd in atomiccmd._commands:
            _collect_stats(subcmd, stats)
    else:
        assert False # pragma: no coverage

    return stats


def _build_status(atomiccmd, _stats, indent, lines):
    prefix = " " * indent + "Status  = "
    if atomiccmd._proc:
        if atomiccmd.ready():
            return_code = tuple(atomiccmd.join())
            if isinstance(return_code[0], types.StringTypes):
                lines.append(prefix + "Terminated with signal %s" % return_code)
            else:
                lines.append(prefix + "Exited with return-code %i" % return_code)
        else:
            lines.append(prefix + "Running ...")


def _build_stdin(atomiccmd, files, stats, indent, lines):
    pipe_name = _get_pipe_name(files, "IN_STDIN")
    pipe   = _get_pipe_file(files, "IN_STDIN")
    prefix = "%s%s  = " % (" " * indent, pipe_name)
    if pipe and pipe in stats["id"]:
        lines.append("%s<%02i>" % (prefix, stats["id"][pipe],))
    elif isinstance(pipe, types.StringTypes):
        if atomiccmd._set_cwd and (pipe_name == "STDIN*"):
            pipe = os.path.basename(pipe)
        lines.append("%s'%s'" % (prefix, pipe))
    elif pipe:
        lines.append("%s<PIPE>" % (prefix,))


def _build_out_pipe(atomiccmd, files, stats, indent, lines, pipe):
    pipe_name = _get_pipe_name(files, pipe)
    prefix = "%s%s = " % (" " * indent, pipe_name)

    if (atomiccmd in stats["pipe"]) and (pipe == "OUT_STDOUT"):
        pipe = stats["pipe"].get(atomiccmd)
        lines.append("%s<%02i>" % (prefix, stats["id"][pipe],))
        return

    filename = _get_pipe_file(files, pipe)
    if filename is not subprocess.PIPE:
        lines.append("%s'%s'" % (prefix, filename))
    else:
        lines.append("%s<PIPE>" % (prefix,))


def _build_cwd(atomiccmd, indent, lines):
    prefix = " " * indent + "CWD     = "
    if atomiccmd._temp:
        if atomiccmd._set_cwd:
            lines.append("%s'%s'" % (prefix, atomiccmd._temp,))
        else:
            lines.append("%s'%s'" % (prefix, os.getcwd()))
    elif atomiccmd._set_cwd:
        lines.append("%s'%s'" % (prefix, "${TEMP_DIR}"))


def _pformat(atomiccmd, stats, indent, lines, include_prefix = True):
    s_prefix     = ""
    if include_prefix:
        s_prefix = " " * indent + "- "
        if _is_cls(atomiccmd, "AtomicCmd"):
            cmd_id = stats["id"][atomiccmd]
            s_prefix += "<%02i> " % (cmd_id,)
    s_prefix_len = len(s_prefix)

    if _is_cls(atomiccmd, "AtomicCmd"):
        temp = "" if atomiccmd._set_cwd else (atomiccmd._temp or "${TEMP_DIR}")
        files = atomiccmd._generate_filenames(atomiccmd._files, temp)

        c_prefix = s_prefix + "Command = "
        for line in _pformat_list(atomiccmd._generate_call(temp)).split("\n"):
            lines.append("%s%s" % (c_prefix, line))
            c_prefix = " " * len(c_prefix)

        if not s_prefix_len:
            s_prefix_len += 1

        _build_status(atomiccmd,   stats, s_prefix_len, lines)
        _build_stdin(atomiccmd,    files, stats, s_prefix_len, lines)
        _build_out_pipe(atomiccmd, files, stats, s_prefix_len, lines, "OUT_STDOUT")
        _build_out_pipe(atomiccmd, files, stats, s_prefix_len, lines, "OUT_STDERR")
        _build_cwd(atomiccmd,                    s_prefix_len, lines)
    elif _is_cls(atomiccmd, "ParallelCmds", "SequentialCmds"):
        lines.append("%s%s:" % (s_prefix, _describe_cls(atomiccmd)))
        for subcmd in atomiccmd._commands:
            _pformat(subcmd, stats, s_prefix_len + 2, lines)
    else:
        assert False # pragma: no coverage


def _pformat_list(lst, width = 80):
    """Return a printable representation of a list, where line-breaks
    are inserted between items to minimize the number of lines with a
    width greater than 'width'. Very long items may cause this maximum
    to be exceeded."""
    result = [[]]
    current_width = 0
    for item in map(repr, lst):
        if current_width + len(item) + 2 > width:
            if not result[-1]:
                result[-1] = [item]
                current_width = len(item) + 2
            else:
                result.append([item])
                current_width = len(item) + 2
        else:
            result[-1].append(item)
            current_width += len(item) + 2

    return "[%s]" % (",\n ".join(", ".join(line) for line in result))


def pformat(atomiccmd):
    """Returns a human readable description of an Atomic Cmd or Atomic Set
    of commands. This is currently equivalent to str(cmd_obj)."""
    if not _is_cls(atomiccmd, "AtomicCmd", "ParallelCmds", "SequentialCmds"):
        raise TypeError("Invalid type in pformat: %r" % atomiccmd.__class__.__name__)

    lines = []
    stats = _collect_stats(atomiccmd, {"id" : {}, "pipe" : {}})
    _pformat(atomiccmd, stats, 0, lines, False)
    return "<%s>" % "\n".join(lines)


def pprint(atomiccmd, out = sys.stdout):
    """Prints a human readable description of an Atomic Cmd or Atomic Set
    of commands. This is currently equivalent to print(str(cmd_obj), ...)."""
    print(pformat(atomiccmd), file = out)

