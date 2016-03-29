#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is herby granted, free of charge, to any person obtaining a copy
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
# Disable warning for missing docstring
# pylint: disable=C0111
# Disable warning caused by "invalid" function names
# pylint: disable=C0103
# Disable warning caused by touching private member variables/functions
# TODO: Remove this / fix places touching privates
# pylint: disable=W0212
import os

from flexmock import \
    flexmock

from paleomix.common.testing import \
    with_temp_folder, \
    set_file_contents

from paleomix.nodegraph import \
    NodeGraph, \
    FileStatusCache


def test_dir():
    return os.path.dirname(__file__)


def test_file(*args):
    return os.path.join(test_dir(), "data", *args)


_DESCRIPTION = "My description of a node"
_IN_FILES = frozenset((test_file("empty_file_1"),
                       test_file("empty_file_2")))
_OUT_FILES = frozenset((test_file("missing_out_file_1"),
                        test_file("missing_out_file_2")))
_EXEC_FILES = frozenset(("ls", "sh"))
_AUX_FILES = frozenset((test_file("rCRS.fasta"),
                        test_file("rCRS.fasta.fai")))
_REQUIREMENTS = frozenset((id, str))


###############################################################################
###############################################################################
# Setup timestamps for test files

def setup_module():
    timestamps = {test_file("timestamp_a_older"): 1000190760,
                  test_file("timestamp_b_older"): 1000190760,
                  test_file("timestamp_a_younger"): 1120719000,
                  test_file("timestamp_b_younger"): 1120719000}

    for filename, timestamp in timestamps.iteritems():
        # Set atime and mtime
        os.utime(filename, (timestamp, timestamp))


###############################################################################
###############################################################################
# NodeGraph: _is_done
# TODO: Avoid testing private function, mock cache

def test_nodegraph_is_done__no_output():
    cache = FileStatusCache()
    node = flexmock(output_files=())
    assert NodeGraph.is_done(node, cache)


@with_temp_folder
def test_nodegraph_is_done__output_changes(temp_folder):
    temp_file_1 = os.path.join(temp_folder, "file_1.txt")
    temp_file_2 = os.path.join(temp_folder, "file_2.txt")
    my_node = flexmock(output_files=(temp_file_1, temp_file_2))
    assert not NodeGraph.is_done(my_node, FileStatusCache())
    set_file_contents(temp_file_1, "foo")
    assert not NodeGraph.is_done(my_node, FileStatusCache())
    set_file_contents(temp_file_2, "bar")
    assert NodeGraph.is_done(my_node, FileStatusCache())


@with_temp_folder
def test_nodegraph_is_done__subnode_not_considered(temp_folder):
    temp_file = os.path.join(temp_folder, "file.txt")
    subnode = flexmock(output_files=(temp_file,))
    my_node = flexmock(output_files=(),
                       subnodes=(subnode,))
    assert NodeGraph.is_done(my_node, FileStatusCache())


def test_nodegraph_is_outdated__no_output():
    my_node = flexmock(input_files=(),
                       output_files=())
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__input_but_no_output():
    my_node = flexmock(input_files=_IN_FILES,
                       output_files=())
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__output_but_no_input():
    my_node = flexmock(input_files=(),
                       output_files=_OUT_FILES)
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__not_outdated():
    my_node = flexmock(input_files=(test_file("timestamp_a_older"),),
                       output_files=(test_file("timestamp_a_younger"),))
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__outdated():
    my_node = flexmock(input_files=(test_file("timestamp_a_younger"),),
                       output_files=(test_file("timestamp_a_older"),))
    assert NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__updates():
    my_node = flexmock(input_files=(test_file("timestamp_a_older"),),
                       output_files=(test_file("timestamp_a_younger"),))
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())
    my_node = flexmock(input_files=(test_file("timestamp_a_younger"),),
                       output_files=(test_file("timestamp_a_older"),))
    assert NodeGraph.is_outdated(my_node, FileStatusCache())
