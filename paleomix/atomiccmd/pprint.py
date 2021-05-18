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
import os
import shlex


def _is_cls(obj, *cls_names):
    return obj.__class__.__name__ in cls_names


def _describe_cls(atomiccmd):
    if _is_cls(atomiccmd, "ParallelCmds"):
        return "Parallel processes"
    elif _is_cls(atomiccmd, "SequentialCmds"):
        return "Sequential processes"
    assert False  # pragma: no coverage


def _collect_stats(atomiccmd, stats):
    assert atomiccmd not in stats["id"]

    if _is_cls(atomiccmd, "AtomicCmd"):
        stats["id"][atomiccmd] = len(stats["id"]) + 1
        if _is_cls(atomiccmd._stdin, "AtomicCmd"):
            stats["pipe"][atomiccmd._stdin] = atomiccmd
    elif _is_cls(atomiccmd, "ParallelCmds", "SequentialCmds"):
        for subcmd in atomiccmd._commands:
            _collect_stats(subcmd, stats)
    else:
        assert False  # pragma: no coverage

    return stats


def _build_status(atomiccmd, _stats, indent, lines):
    prefix = " " * indent + "Status  = "
    if atomiccmd._proc:
        if atomiccmd.ready():
            return_code = tuple(atomiccmd.join())
            if atomiccmd._terminated:
                lines.append(prefix + "Automatically terminated by PALEOMIX")
            elif isinstance(return_code[0], str):
                lines.append(prefix + "Terminated with signal %s" % return_code)
            else:
                lines.append(prefix + "Exited with return-code %i" % return_code)
        else:
            lines.append(prefix + "Running")


def _build_stdin(atomiccmd, stats, indent, lines):
    pipe = atomiccmd._stdin
    prefix = "%s%s   = " % (" " * indent, "STDIN")
    if atomiccmd._stdin in stats["id"]:
        lines.append("%sPiped from process %i" % (prefix, stats["id"][pipe]))
    elif _is_cls(pipe, "InputFile", "TempInputFile"):
        temp = "${TEMP_DIR}" if atomiccmd._temp is None else atomiccmd._temp
        path = atomiccmd._to_path(temp, pipe)
        lines.append("%s%s" % (prefix, shlex.quote(path)))
    elif pipe:  # FIXME: DEVNULL/PIPE
        lines.append("%s<PIPE>" % (prefix,))


def _build_stdout(atomiccmd, stats, indent, lines):
    prefix = "%sSTDOUT  = " % (" " * indent,)

    pipe = atomiccmd._stdout
    if atomiccmd in stats["pipe"]:
        pipe = stats["pipe"][atomiccmd]
        lines.append("%sPiped to process %i" % (prefix, stats["id"][pipe]))
    elif _is_cls(atomiccmd._stdout, "OutputFile", "TempOutputFile"):
        temp = "${TEMP_DIR}" if atomiccmd._temp is None else atomiccmd._temp
        path = atomiccmd._to_path(temp, pipe)

        lines.append("%s%s" % (prefix, shlex.quote(path)))
    elif pipe == atomiccmd.PIPE:
        lines.append("%s<PIPE>" % (prefix,))
    elif pipe == atomiccmd.DEVNULL:
        lines.append("%s<DEVNULL>" % (prefix,))


def _build_stderr(atomiccmd, stats, indent, lines):
    prefix = "%sSTDERR  = " % (" " * indent,)

    pipe = atomiccmd._stdout
    if _is_cls(atomiccmd._stderr, "OutputFile", "TempOutputFile"):
        pipe = atomiccmd._stderr
        temp = "${TEMP_DIR}" if atomiccmd._temp is None else atomiccmd._temp
        path = atomiccmd._to_path(temp, pipe)

        lines.append("%s%s" % (prefix, shlex.quote(path)))
    elif pipe == atomiccmd.PIPE:
        lines.append("%s<PIPE>" % (prefix,))
    elif pipe == atomiccmd.DEVNULL:
        lines.append("%s<DEVNULL>" % (prefix,))


def _build_cwd(atomiccmd, indent, lines):
    prefix = " " * indent + "CWD     = "
    if atomiccmd._temp:
        if atomiccmd._set_cwd:
            lines.append("%s%s" % (prefix, shlex.quote(atomiccmd._temp)))
        else:
            lines.append("%s%s" % (prefix, shlex.quote(os.getcwd())))
    elif atomiccmd._set_cwd:
        lines.append("%s%s" % (prefix, shlex.quote("${TEMP_DIR}")))


def _pformat(atomiccmd, stats, indent, lines, include_prefix=True):
    s_prefix = ""
    if include_prefix:
        s_prefix = " " * indent
        if _is_cls(atomiccmd, "AtomicCmd"):
            cmd_id = stats["id"][atomiccmd]
            lines.append(s_prefix + "Process %i:" % (cmd_id,))
            s_prefix += "  "
    s_prefix_len = len(s_prefix)

    if _is_cls(atomiccmd, "AtomicCmd"):
        print(atomiccmd._set_cwd, atomiccmd._temp)
        temp = "" if atomiccmd._set_cwd else (atomiccmd._temp or "${TEMP_DIR}")

        c_prefix = s_prefix + "Command = "
        for line in _pformat_list(atomiccmd.to_call(temp)).split("\n"):
            lines.append("%s%s" % (c_prefix, line))
            c_prefix = " " * len(c_prefix)

        _build_status(atomiccmd, stats, s_prefix_len, lines)
        _build_stdin(atomiccmd, stats, s_prefix_len, lines)
        _build_stdout(atomiccmd, stats, s_prefix_len, lines)
        _build_stderr(atomiccmd, stats, s_prefix_len, lines)
        _build_cwd(atomiccmd, s_prefix_len, lines)
    elif _is_cls(atomiccmd, "ParallelCmds", "SequentialCmds"):
        lines.append("%s%s:" % (s_prefix, _describe_cls(atomiccmd)))
        for subcmd_idx, subcmd in enumerate(atomiccmd._commands):
            if subcmd_idx:
                lines.append("")

            _pformat(subcmd, stats, s_prefix_len + 2, lines)
    else:
        assert False  # pragma: no coverage


def _pformat_list(lst, width=80):
    """Return a printable representation of a list, where line-breaks
    are inserted between items to minimize the number of lines with a
    width greater than 'width'. Very long items may cause this maximum
    to be exceeded."""
    result = [[]]
    current_width = 0
    for item in (shlex.quote(str(value)) for value in lst):
        if current_width + len(item) + 1 > width:
            if not result[-1]:
                result[-1] = [item]
            else:
                result.append([item])

            current_width = len(item) + 1
        else:
            result[-1].append(item)
            current_width += len(item) + 1

    return " \\\n    ".join(" ".join(line) for line in result)


def pformat(atomiccmd):
    """Returns a human readable description of an Atomic Cmd or Atomic Set
    of commands. This is currently equivalent to str(cmd_obj)."""
    if not _is_cls(atomiccmd, "AtomicCmd", "ParallelCmds", "SequentialCmds"):
        raise TypeError("Invalid type in pformat: %r" % atomiccmd.__class__.__name__)

    lines = []
    stats = _collect_stats(atomiccmd, {"id": {}, "pipe": {}})
    _pformat(atomiccmd, stats, 0, lines, False)

    return "\n".join(lines)
