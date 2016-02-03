#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
from __future__ import print_function

import sys


# Global setting to enable, disable, or force the use console codes when the
# print_* functions are used. By default (ON), color codes will only be used if
# the destination is a TTY.
COLORS_OFF, COLORS_ON, COLORS_FORCED = range(3)


def set_color_output(clr):
    """Set use of colors in print functions; possible values are COLORS_OFF,
    COLORS_ON, and COLORS_FORCED.
    """
    global _COLORS
    if clr not in (COLORS_OFF, COLORS_ON, COLORS_FORCED):
        raise ValueError("Invalid value in set_color_output; must be one of "
                         "COLORS_OFF, COLORS_ON, or COLORS_FORCED.")

    _COLORS = clr


def print_msg(*vargs, **kwargs):
    """Equivalent to print. Currently does not apply a color to the text"""
    print(*vargs, **kwargs)


def print_debug(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (blue)."""
    _do_print_color(*vargs, colorcode=36, **kwargs)


def print_info(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (green)."""
    _do_print_color(*vargs, colorcode=32, **kwargs)


def print_err(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (red)."""
    _do_print_color(*vargs, colorcode=31, **kwargs)


def print_warn(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (yellow)."""
    _do_print_color(*vargs, colorcode=33, **kwargs)


def print_disabled(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (gray)."""
    _do_print_color(*vargs, colorcode=30, **kwargs)


def _do_print_color(*vargs, **kwargs):
    """Utility function: Prints using shell colors."""
    colorcode = kwargs.pop("colorcode")
    destination = kwargs.pop("file", sys.stderr)

    # No colors if output is redirected (e.g. less, file, etc.)
    colors_on = _COLORS != COLORS_OFF
    colors_forced = _COLORS == COLORS_FORCED
    if colors_on and (destination.isatty() or colors_forced):
        vargs = list(vargs)
        for (index, varg) in enumerate(vargs):
            varg_lines = []
            # Newlines terminate the color-code for e.g. 'less', so ensure that
            # each line is color-coded, while preserving the list of arguments
            for line in str(varg).split("\n"):
                varg_lines.append("\033[00;%im%s\033[00m" % (colorcode, line))
            vargs[index] = "\n".join(varg_lines)

    print(*vargs, file=destination, **kwargs)

    if '\n' in kwargs.get('end', '\n'):
        destination.flush()

_COLORS = COLORS_ON
