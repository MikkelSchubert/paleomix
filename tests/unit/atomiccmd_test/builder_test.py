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
# pylint: disable=C0103

from flexmock import flexmock
from nose.tools import assert_equal, assert_raises # pylint: disable=E0611

from tests.common.utils import monkeypatch

from pypeline.atomiccmd.builder import \
     AtomicCmdBuilder, \
     AtomicCmdBuilderError, \
     AtomicJavaCmdBuilder, \
     AtomicMPICmdBuilder



################################################################################
################################################################################
## AtomicCmdBuilder: Constructor


def test_builder__simple__call():
    builder = AtomicCmdBuilder(["ls"])
    assert_equal(builder.call, ["ls"])

def test_builder__simple__str():
    builder = AtomicCmdBuilder("ls")
    assert_equal(builder.call, ["ls"])

def test_builder__simple__iterable():
    builder = AtomicCmdBuilder(iter(["ls"]))
    assert_equal(builder.call, ["ls"])

def test_builder__complex():
    builder = AtomicCmdBuilder(["java", "jar", "/a/jar"])
    assert_equal(builder.call, ["java", "jar", "/a/jar"])

def test_builder__kwargs__empty():
    builder = AtomicCmdBuilder(["ls"])
    assert_equal(builder.kwargs, {})

def test_builder__kwargs():
    expected = {"IN_FILE"  : "/abc/def.txt",
                "OUT_FILE" : "/etc/fstab"}
    builder = AtomicCmdBuilder(["ls"], **expected)
    assert_equal(builder.kwargs, expected)

def test_builder__kwargs__set_cwd():
    builder = AtomicCmdBuilder(["ls"], set_cwd = True)
    assert_equal(builder.kwargs, {"set_cwd" : True})



################################################################################
################################################################################
## AtomicCmdBuilder: set_option


def test_builder__set_option():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt")
    assert_equal(builder.call, ["find", "-name", "*.txt"])

def test_builder__set_option__overwrite():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt", fixed = False)
    builder.set_option("-name", "*.bat")
    assert_equal(builder.call, ["find", "-name", "*.bat"])

def test_builder__set_option__overwrite_fixed():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt")
    assert_raises(AtomicCmdBuilderError, builder.set_option, "-name", "*.bat")




################################################################################
################################################################################
## AtomicCmdBuilder: add_option

def test_builder__add_option():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-name", "*.txt")
    assert_equal(builder.call, ["find", "-name", "*.txt"])

def test_builder__add_option__overwrite():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-name", "*.txt")
    builder.add_option("-or")
    builder.add_option("-name", "*.bat")
    assert_equal(builder.call, ["find", "-name", "*.txt", "-or", "-name", "*.bat"])



################################################################################
################################################################################
## AtomicCmdBuilder: add_option / set_option common tests


def test_builder__add_or_set_option__without_value():
    def _do_test_builder__add_or_set_option__without_value(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-delete")
        assert_equal(builder.call, ["find", "-delete"])
    yield _do_test_builder__add_or_set_option__without_value, AtomicCmdBuilder.add_option
    yield _do_test_builder__add_or_set_option__without_value, AtomicCmdBuilder.set_option


def test_builder__add_or_set_option__with_sep():
    def _do_test_builder__add_or_set_option__with_sep(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-size", "0", sep = "=")
        assert_equal(builder.call, ["find", "-size=0"])
    yield _do_test_builder__add_or_set_option__with_sep, AtomicCmdBuilder.add_option
    yield _do_test_builder__add_or_set_option__with_sep, AtomicCmdBuilder.set_option

def test_builder__add_or_set_option__with_non_str_value():
    def _do_test_test_builder__add_or_set_option__with_non_str_value(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-size", 0)
        assert_equal(builder.call, ["find", "-size", 0])
    yield _do_test_test_builder__add_or_set_option__with_non_str_value, AtomicCmdBuilder.add_option
    yield _do_test_test_builder__add_or_set_option__with_non_str_value, AtomicCmdBuilder.set_option

def test_builder__add_or_set_option__add_and_set():
    def _do_test_builder__add_or_set_option__add_and_set(setter_1, setter_2):
        builder = AtomicCmdBuilder("find")
        setter_1(builder, "-name", "*.txt")
        assert_raises(AtomicCmdBuilderError, setter_2, builder, "-name", "*.bat")
    yield _do_test_builder__add_or_set_option__add_and_set, AtomicCmdBuilder.set_option, AtomicCmdBuilder.add_option
    yield _do_test_builder__add_or_set_option__add_and_set, AtomicCmdBuilder.add_option, AtomicCmdBuilder.set_option

def test_builder__add_or_set_option__with_non_str_key():
    def _do_test_builder__add_or_set_option__with_non_str_key(setter):
        builder = AtomicCmdBuilder("find")
        assert_raises(TypeError, setter, builder, 7913, "True")
    yield _do_test_builder__add_or_set_option__with_non_str_key, AtomicCmdBuilder.add_option
    yield _do_test_builder__add_or_set_option__with_non_str_key, AtomicCmdBuilder.set_option

def test_builder__add_or_set_option__after_finalize():
    def _do_test_builder__add_or_set_option__after_finalize(setter):
        builder = AtomicCmdBuilder("find")
        builder.finalize()
        assert_raises(AtomicCmdBuilderError, setter, builder, "-size", "1")
    yield _do_test_builder__add_or_set_option__after_finalize, AtomicCmdBuilder.add_option
    yield _do_test_builder__add_or_set_option__after_finalize, AtomicCmdBuilder.set_option


def test_builder__add_or_set_option__empty_key():
    def _do_test_builder__add_or_set_option__empty_key(setter):
        builder = AtomicCmdBuilder("find")
        assert_raises(KeyError, setter, builder, "", "1")
    yield _do_test_builder__add_or_set_option__empty_key, AtomicCmdBuilder.add_option
    yield _do_test_builder__add_or_set_option__empty_key, AtomicCmdBuilder.set_option



################################################################################
################################################################################
## AtomicCmdBuilder: pop_option

def test_builder__pop_option():
    def _do_test_builder__pop_option(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-size", "0", fixed = False)
        builder.pop_option("-size")
        assert_equal(builder.call, ["find"])
    yield _do_test_builder__pop_option, AtomicCmdBuilder.set_option
    yield _do_test_builder__pop_option, AtomicCmdBuilder.add_option


def test_builder__pop_option__last_option():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-size", "0", fixed = False)
    builder.add_option("-size", "1", fixed = False)
    builder.pop_option("-size")
    assert_equal(builder.call, ["find", "-size", "0"])

def test_builder__pop_option__different_options():
    def _do_test_builder__pop_option(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-empty", fixed = False)
        setter(builder, "-size", "1", fixed = False)
        setter(builder, "-name", "*.txt", fixed = False)
        builder.pop_option("-size")
        assert_equal(builder.call, ["find", "-empty", "-name", "*.txt"])
    yield _do_test_builder__pop_option, AtomicCmdBuilder.set_option
    yield _do_test_builder__pop_option, AtomicCmdBuilder.add_option



def test_builder__pop_option__is_fixed():
    def _do_test_builder__pop_option__is_fixed(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-size", "0")
        assert_raises(AtomicCmdBuilderError, builder.pop_option, "-size")
    yield _do_test_builder__pop_option__is_fixed, AtomicCmdBuilder.set_option
    yield _do_test_builder__pop_option__is_fixed, AtomicCmdBuilder.add_option

def test_builder__pop_option__empty():
    builder = AtomicCmdBuilder("find")
    assert_raises(KeyError, builder.pop_option, "-size")

def test_builder__pop_option__missing_key():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-size", 0)
    assert_raises(KeyError, builder.pop_option, "-isize")

def test_builder__pop_option__with_non_str_key():
    builder = AtomicCmdBuilder("find")
    assert_raises(TypeError, builder.pop_option, 7913)

def test_builder__pop_option__with_empty_key():
    builder = AtomicCmdBuilder("find")
    assert_raises(KeyError, builder.pop_option, "")




################################################################################
################################################################################
## AtomicCmdBuilder: add_value

def test_builder__add_value():
    builder = AtomicCmdBuilder("ls")
    builder.add_value("%(IN_FILE)s")
    assert_equal(builder.call, ["ls", "%(IN_FILE)s"])


def test_builder__add_value__two_values():
    builder = AtomicCmdBuilder("ls")
    builder.add_value("%(IN_FILE)s")
    builder.add_value("%(OUT_FILE)s")
    assert_equal(builder.call, ["ls", "%(IN_FILE)s", "%(OUT_FILE)s"])




################################################################################
################################################################################
## AtomicCmdBuilder: set_kwargs


def test_builder__set_kwargs__called_once():
    expected = {"IN_PATH" : "/a/b/", "OUT_PATH" : "/dst/file"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(**expected)
    assert_equal(builder.kwargs, expected)

def test_builder__set_kwargs__called_twice():
    expected = {"IN_PATH" : "/a/b/", "OUT_PATH" : "/dst/file"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(OUT_PATH = "/dst/file")
    builder.set_kwargs(IN_PATH = "/a/b/")
    assert_equal(builder.kwargs, expected)

def test_builder__set_kwargs__atomiccmdbuilder():
    mock = flexmock(AtomicCmdBuilder("true"))
    mock.should_receive('finalize').and_return("finalized!")
    builder = AtomicCmdBuilder("ls", IN_BUILDER = mock)
    assert_equal(builder.kwargs, {"IN_BUILDER" : "finalized!"})

def test_builder__set_kwargs__after_finalize():
    expected = {"IN_PATH" : "/a/b/"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(IN_PATH = "/a/b/")
    builder.finalize()
    assert_raises(AtomicCmdBuilderError, builder.set_kwargs, OUT_PATH = "/dst/file")
    assert_equal(builder.kwargs, expected)

def test_builder__set__kwargs__overwriting():
    expected = {"IN_PATH" : "/a/b/"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(IN_PATH = "/a/b/")
    assert_raises(AtomicCmdBuilderError, builder.set_kwargs, IN_PATH = "/dst/file")
    assert_equal(builder.kwargs, expected)




################################################################################
################################################################################
## AtomicCmdBuilder: finalize

def test_builder__finalize__returns_singleton():
    builder = AtomicCmdBuilder("echo")
    assert builder.finalize() is builder.finalize()

def test_builder__finalize__calls_atomiccmd():
    was_called = []
    class _AtomicCmdMock:
        def __init__(self, *args, **kwargs):
            assert_equal(args, (["echo", "-out", "%(OUT_FILE)s", "%(IN_FILE)s"],))
            assert_equal(kwargs, {"IN_FILE" : "/in/file", "OUT_FILE" : "/out/file", "set_cwd" : True})
            was_called.append(True)

    with monkeypatch("pypeline.atomiccmd.builder.AtomicCmd", _AtomicCmdMock):
        builder = AtomicCmdBuilder("echo", set_cwd = True)
        builder.add_option("-out", "%(OUT_FILE)s")
        builder.add_value("%(IN_FILE)s")
        builder.set_kwargs(OUT_FILE = "/out/file",
                           IN_FILE  = "/in/file")

        builder.finalize()
        assert was_called



################################################################################
################################################################################
## AtomicJavaCmdBuilder

JAVA_CFG = flexmock(temp_root = "/disk/tmp")

def test_java_builder__defaults__call():
    builder = AtomicJavaCmdBuilder(JAVA_CFG, "/path/Foo.jar")
    assert_equal(builder.call, ["java", "-server", "-Xmx4g",
                                "-Djava.io.tmpdir=/disk/tmp",
                                "-XX:+UseSerialGC",
                                "-jar", "%(AUX_JAR)s"])

def test_java_builder__defaults__kwargs():
    builder = AtomicJavaCmdBuilder(JAVA_CFG, "/path/Foo.jar")
    assert_equal(builder.kwargs, {"AUX_JAR" : "/path/Foo.jar"})

def test_java_builder__multithreaded_gc():
    builder = AtomicJavaCmdBuilder(JAVA_CFG, "/path/Foo.jar", gc_threads = 3)
    assert_equal(builder.call, ["java", "-server", "-Xmx4g",
                                "-Djava.io.tmpdir=/disk/tmp",
                                "-XX:ParallelGCThreads=3",
                                "-jar", "%(AUX_JAR)s"])

def test_java_builder__multithreaded_gc__zero_or_negative_threads():
    assert_raises(ValueError, AtomicJavaCmdBuilder, JAVA_CFG, "/path/Foo.jar", gc_threads = 0)
    assert_raises(ValueError, AtomicJavaCmdBuilder, JAVA_CFG, "/path/Foo.jar", gc_threads = -1)

def test_java_builder__multithreaded_gc__non_int_threads():
    assert_raises(TypeError, AtomicJavaCmdBuilder, JAVA_CFG, "/path/Foo.jar", gc_threads = "3")

def test_java_builder__kwargs():
    builder = AtomicJavaCmdBuilder(JAVA_CFG, "/path/Foo.jar", set_cwd = True)
    assert_equal(builder.kwargs, {"AUX_JAR" : "/path/Foo.jar", "set_cwd" : True})



################################################################################
################################################################################
## AtomicMPICmdBuilder

def test_mpi_builder__defaults__str():
    builder = AtomicMPICmdBuilder("ls")
    assert_equal(builder.call, ["ls"])
    assert_equal(builder.kwargs, {})

def test_mpi_builder__multithreaded__str():
    builder = AtomicMPICmdBuilder("ls", threads = 3)
    assert_equal(builder.call, ["mpirun", "-n", 3, "ls"])
    assert_equal(builder.kwargs, {"EXEC_MAIN" : "ls"})

def test_mpi_builder__defaults__complex_cmd():
    builder = AtomicMPICmdBuilder(["python", "/foo/run.py"])
    assert_equal(builder.call, ["python", "/foo/run.py"])
    assert_equal(builder.kwargs, {})

def test_mpi_builder__multithreaded__complex_cmd():
    builder = AtomicMPICmdBuilder(["python", "/foo/run.py"], threads = 3)
    assert_equal(builder.call, ["mpirun", "-n", 3, "python", "/foo/run.py"])
    assert_equal(builder.kwargs, {"EXEC_MAIN" : "python"})

def test_mpi_builder__kwargs():
    builder = AtomicMPICmdBuilder("ls", set_cwd = True)
    assert_equal(builder.kwargs, {"set_cwd" : True})

def test_mpi_builder__threads__zero_or_negative():
    assert_raises(ValueError, AtomicMPICmdBuilder, "ls", threads =  0)
    assert_raises(ValueError, AtomicMPICmdBuilder, "ls", threads = -1)

def test_mpi_builder__threads__non_int():
    assert_raises(TypeError, AtomicMPICmdBuilder, "ls",  threads = "3")
