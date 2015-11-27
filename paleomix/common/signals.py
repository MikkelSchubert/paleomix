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
import signal
import types


def from_str(value):
    if not isinstance(value, types.StringTypes):
        raise TypeError("'from_str' takes strings, not %r" % (value.__class__.__name__,))
    return FROM_STRING[value]


def to_str(value):
    if not isinstance(value, (types.IntType, types.LongType)):
        raise TypeError("'from_str' takes strings, not %r" % (value.__class__.__name__,))
    return FROM_SIGNAL[value]


def _get_signals():
    signals = {}
    for key in dir(signal):
        if key.startswith("SIG") and not key.startswith("SIG_"):
            signals[getattr(signal, key)] = key

    # The following signals may have synonyms, so specify which to use.
    # For example, SIGIOT is a synonym of SIGABTR on OSX.
    for key in ("SIGABRT", "SIGCHLD", "SIGIO"):
        value = getattr(signal, key)
        signals[value] = key

    return signals

FROM_SIGNAL = _get_signals()
FROM_STRING = dict(zip(FROM_SIGNAL.values(), FROM_SIGNAL.keys()))
