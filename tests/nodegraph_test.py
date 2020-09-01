#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
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
import os

from unittest.mock import Mock

from paleomix.common.fileutils import fspath
from paleomix.nodegraph import NodeGraph, FileStatusCache


_TIMESTAMP_1 = 1000190760
_TIMESTAMP_2 = 1120719000


def create_test_file(utime, *args):
    filename = os.path.join(*map(fspath, args))
    with open(filename, "wb"):
        pass

    os.utime(filename, (utime, utime))

    return filename


###############################################################################
###############################################################################
# NodeGraph: _is_done
# TODO: Avoid testing private function, mock cache


def test_nodegraph_is_done__no_output():
    cache = FileStatusCache()
    node = Mock(output_files=())
    assert NodeGraph.is_done(node, cache)


def test_nodegraph_is_done__output_changes(tmp_path):
    temp_file_1 = tmp_path / "file_1.txt"
    temp_file_2 = tmp_path / "file_2.txt"
    my_node = Mock(output_files=(str(temp_file_1), str(temp_file_2)))
    assert not NodeGraph.is_done(my_node, FileStatusCache())
    temp_file_1.write_text("foo")
    assert not NodeGraph.is_done(my_node, FileStatusCache())
    temp_file_2.write_text("bar")
    assert NodeGraph.is_done(my_node, FileStatusCache())


def test_nodegraph_is_done__subnode_not_considered(tmp_path):
    temp_file = tmp_path / "file.txt"
    subnode = Mock(output_files=(str(temp_file),))
    my_node = Mock(output_files=(), subnodes=(subnode,))
    assert NodeGraph.is_done(my_node, FileStatusCache())


def test_nodegraph_is_outdated__no_output():
    my_node = Mock(input_files=(), output_files=())
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__input_but_no_output(tmp_path):
    input_file = tmp_path / "file"
    input_file.touch()

    my_node = Mock(input_files=(str(input_file),), output_files=())
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__output_but_no_input(tmp_path):
    output_file = tmp_path / "file"
    output_file.touch()

    my_node = Mock(input_files=(), output_files=(str(output_file),))
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__not_outdated(tmp_path):
    my_node = Mock(
        input_files=(create_test_file(_TIMESTAMP_1, tmp_path, "older_file"),),
        output_files=(create_test_file(_TIMESTAMP_2, tmp_path, "younger_file"),),
    )
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__outdated(tmp_path):
    my_node = Mock(
        input_files=(create_test_file(_TIMESTAMP_2, tmp_path, "younger_file"),),
        output_files=(create_test_file(_TIMESTAMP_1, tmp_path, "older_file"),),
    )
    assert NodeGraph.is_outdated(my_node, FileStatusCache())


def test_nodegraph_is_outdated__updates(tmp_path):
    older_file = create_test_file(_TIMESTAMP_1, tmp_path, "older_file")
    younger_file = create_test_file(_TIMESTAMP_2, tmp_path, "younger_file")

    my_node = Mock(input_files=(older_file,), output_files=(younger_file,),)
    assert not NodeGraph.is_outdated(my_node, FileStatusCache())
    my_node = Mock(input_files=(younger_file,), output_files=(older_file,),)
    assert NodeGraph.is_outdated(my_node, FileStatusCache())
