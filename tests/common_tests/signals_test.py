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
import paleomix.common.signals as signals

import nose
from nose.tools import assert_equal


def test_signal__sigterm_to_str():
    assert_equal(signals.to_str(signal.SIGTERM), "SIGTERM")


def test_signal__str_to_sigterm():
    assert_equal(signals.from_str("SIGTERM"), signal.SIGTERM)


@nose.tools.raises(KeyError)
def test_signal__to_str__unknown_signal():
    signals.to_str(1024)


@nose.tools.raises(KeyError)
def test_signal__from_str__unknown_signal():
    signals.from_str("SIGSMURF")


@nose.tools.raises(TypeError)
def test_signal__to_str__wrong_type():
    signals.to_str("SIGTERM")


@nose.tools.raises(TypeError)
def test_signal__from_str__wrong_type():
    signals.from_str(signal.SIGTERM)
