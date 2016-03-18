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
from flexmock import flexmock
from nose.tools import assert_equal, assert_raises

from paleomix.common.testing import Monkeypatch

from paleomix.atomiccmd.builder import \
    AtomicCmdBuilder, \
    AtomicCmdBuilderError, \
    AtomicJavaCmdBuilder, \
    AtomicMPICmdBuilder, \
    apply_options, \
    use_customizable_cli_parameters, \
    create_customizable_cli_parameters, \
    JAVA_VERSIONS


###############################################################################
###############################################################################
# AtomicCmdBuilder: Constructor


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
    expected = {"IN_FILE": "/abc/def.txt",
                "OUT_FILE": "/etc/fstab"}
    builder = AtomicCmdBuilder(["ls"], **expected)
    assert_equal(builder.kwargs, expected)


def test_builder__kwargs__set_cwd():
    builder = AtomicCmdBuilder(["ls"], set_cwd=True)
    assert_equal(builder.kwargs, {"set_cwd": True})


###############################################################################
###############################################################################
# AtomicCmdBuilder: set_option


def test_builder__set_option():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt")
    assert_equal(builder.call, ["find", "-name", "*.txt"])


def test_builder__set_option__overwrite():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt", fixed=False)
    builder.set_option("-name", "*.bat")
    assert_equal(builder.call, ["find", "-name", "*.bat"])


def test_builder__set_option__overwrite_fixed():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt")
    assert_raises(AtomicCmdBuilderError, builder.set_option, "-name", "*.bat")


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_option

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


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_option / set_option common tests

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
        setter(builder, "-size", "0", sep="=")
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


###############################################################################
###############################################################################
# AtomicCmdBuilder: pop_option

def test_builder__pop_option():
    def _do_test_builder__pop_option(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-size", "0", fixed=False)
        builder.pop_option("-size")
        assert_equal(builder.call, ["find"])
    yield _do_test_builder__pop_option, AtomicCmdBuilder.set_option
    yield _do_test_builder__pop_option, AtomicCmdBuilder.add_option


def test_builder__pop_option__last_option():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-size", "0", fixed=False)
    builder.add_option("-size", "1", fixed=False)
    builder.pop_option("-size")
    assert_equal(builder.call, ["find", "-size", "0"])


def test_builder__pop_option__different_options():
    def _do_test_builder__pop_option(setter):
        builder = AtomicCmdBuilder("find")
        setter(builder, "-empty", fixed=False)
        setter(builder, "-size", "1", fixed=False)
        setter(builder, "-name", "*.txt", fixed=False)
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


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_value

def test_builder__add_value():
    builder = AtomicCmdBuilder("ls")
    builder.add_value("%(IN_FILE)s")
    assert_equal(builder.call, ["ls", "%(IN_FILE)s"])


def test_builder__add_value__two_values():
    builder = AtomicCmdBuilder("ls")
    builder.add_value("%(IN_FILE)s")
    builder.add_value("%(OUT_FILE)s")
    assert_equal(builder.call, ["ls", "%(IN_FILE)s", "%(OUT_FILE)s"])


###############################################################################
###############################################################################
# AtomicCmdBuilder: set_kwargs

def test_builder__set_kwargs__called_once():
    expected = {"IN_PATH": "/a/b/", "OUT_PATH": "/dst/file"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(**expected)
    assert_equal(builder.kwargs, expected)


def test_builder__set_kwargs__called_twice():
    expected = {"IN_PATH": "/a/b/", "OUT_PATH": "/dst/file"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(OUT_PATH="/dst/file")
    builder.set_kwargs(IN_PATH="/a/b/")
    assert_equal(builder.kwargs, expected)


def test_builder__set_kwargs__atomiccmdbuilder():
    mock = flexmock(AtomicCmdBuilder("true"))
    mock.should_receive('finalize').and_return("finalized!")
    builder = AtomicCmdBuilder("ls", IN_BUILDER=mock)
    assert_equal(builder.kwargs, {"IN_BUILDER": "finalized!"})


def test_builder__set_kwargs__after_finalize():
    expected = {"IN_PATH": "/a/b/"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(IN_PATH="/a/b/")
    builder.finalize()
    assert_raises(AtomicCmdBuilderError, builder.set_kwargs, OUT_PATH="/dst/file")
    assert_equal(builder.kwargs, expected)


def test_builder__set__kwargs__overwriting():
    expected = {"IN_PATH": "/a/b/"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(IN_PATH="/a/b/")
    assert_raises(AtomicCmdBuilderError, builder.set_kwargs, IN_PATH="/dst/file")
    assert_equal(builder.kwargs, expected)


###############################################################################
###############################################################################
# AtomicCmdBuilder: finalize

def test_builder__finalized_call__simple_command():
    builder = AtomicCmdBuilder("echo")
    assert_equal(builder.finalized_call, ["echo"])


def test_builder__finalized_call__kwargs_are_instantiated():
    builder = AtomicCmdBuilder(("echo", "%(ARG1)s", "X=%(ARG2)s"),
                               ARG1="/foo/bar",
                               ARG2="zod")
    assert_equal(builder.finalized_call, ["echo", "/foo/bar", "X=zod"])


def test_builder__finalized_call__kwargs_are_instantiated__with_temp_dir():
    builder = AtomicCmdBuilder(("echo", "%(ARG)s", "%(TEMP_DIR)s"),
                               ARG="/foo/bar")
    assert_equal(builder.finalized_call, ["echo", "/foo/bar", "%(TEMP_DIR)"])


def test_builder__finalized_call__kwargs_are_instantiated__with_non_str_arg():
    builder = AtomicCmdBuilder(("echo", "%(ARG)s", 17),
                               ARG="/foo/bar")
    assert_equal(builder.finalized_call, ["echo", "/foo/bar", "17"])


###############################################################################
###############################################################################
# AtomicCmdBuilder: finalize

def test_builder__finalize__returns_singleton():
    builder = AtomicCmdBuilder("echo")
    assert builder.finalize() is builder.finalize()


def test_builder__finalize__calls_atomiccmd():
    was_called = []

    class _AtomicCmdMock:
        def __init__(self, *args, **kwargs):
            assert_equal(args, (["echo", "-out", "%(OUT_FILE)s", "%(IN_FILE)s"],))
            assert_equal(kwargs, {"IN_FILE": "/in/file",
                                  "OUT_FILE": "/out/file",
                                  "set_cwd": True})
            was_called.append(True)

    with Monkeypatch("paleomix.atomiccmd.builder.AtomicCmd", _AtomicCmdMock):
        builder = AtomicCmdBuilder("echo", set_cwd=True)
        builder.add_option("-out", "%(OUT_FILE)s")
        builder.add_value("%(IN_FILE)s")
        builder.set_kwargs(OUT_FILE="/out/file",
                           IN_FILE="/in/file")

        builder.finalize()
        assert was_called


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_multiple_options

def test_builder__add_multiple_options():
    values = ("file_a", "file_b")
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", values)

    assert_equal(kwargs, expected)
    assert_equal(builder.kwargs, expected)
    assert_equal(builder.call, ["ls",
                                "-i", "%(IN_FILE_01)s",
                                "-i", "%(IN_FILE_02)s"])


def test_builder__add_multiple_options_with_sep():
    values = ("file_a", "file_b")
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", values, sep="=")

    assert_equal(kwargs, expected)
    assert_equal(builder.kwargs, expected)
    assert_equal(builder.call, ["ls",
                                "-i=%(IN_FILE_01)s",
                                "-i=%(IN_FILE_02)s"])


def test_builder__add_multiple_options_with_template():
    values = ("file_a", "file_b")
    expected = {"OUT_BAM_1": "file_a", "OUT_BAM_2": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", values, template="OUT_BAM_%i")

    assert_equal(kwargs, expected)
    assert_equal(builder.kwargs, expected)
    assert_equal(builder.call, ["ls",
                                "-i", "%(OUT_BAM_1)s",
                                "-i", "%(OUT_BAM_2)s"])


def test_builder__add_multiple_options_multiple_times():
    expected = {"IN_FILE_01": "file_a",
                "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", ("file_a",))
    assert_equal(kwargs, {"IN_FILE_01": "file_a"})
    kwargs = builder.add_multiple_options("-i", ("file_b",))
    assert_equal(kwargs, {"IN_FILE_02": "file_b"})

    assert_equal(builder.kwargs, expected)
    assert_equal(builder.call, ["ls",
                                "-i", "%(IN_FILE_01)s",
                                "-i", "%(IN_FILE_02)s"])


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_multiple_values

def test_builder__add_multiple_values():
    values = ("file_a", "file_b")
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_values(values)

    assert_equal(kwargs, expected)
    assert_equal(builder.kwargs, expected)
    assert_equal(builder.call, ["ls", "%(IN_FILE_01)s", "%(IN_FILE_02)s"])


def test_builder__add_multiple_values_with_template():
    values = ("file_a", "file_b")
    expected = {"OUT_BAM_1": "file_a", "OUT_BAM_2": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_values(values, template="OUT_BAM_%i")

    assert_equal(kwargs, expected)
    assert_equal(builder.kwargs, expected)
    assert_equal(builder.call, ["ls", "%(OUT_BAM_1)s", "%(OUT_BAM_2)s"])


def test_builder__add_multiple_values_multiple_times():
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_values(("file_a",))
    assert_equal(kwargs, {"IN_FILE_01": "file_a"})
    kwargs = builder.add_multiple_values(("file_b",))
    assert_equal(kwargs, {"IN_FILE_02": "file_b"})

    assert_equal(builder.kwargs, expected)
    assert_equal(builder.call, ["ls", "%(IN_FILE_01)s", "%(IN_FILE_02)s"])


###############################################################################
###############################################################################
# AtomicJavaCmdBuilder

def test_java_builder__default__no_config():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar")
    assert_equal(builder.call, ["java",
                                "-server",
                                "-Djava.io.tmpdir=%(TEMP_DIR)s",
                                "-Djava.awt.headless=true",
                                "-XX:+UseSerialGC",
                                "-Xmx4g",
                                "-jar", "%(AUX_JAR)s"])


def test_java_builder__defaults__call():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar", temp_root="/disk/tmp")
    assert_equal(builder.call, ["java",
                                "-server",
                                "-Djava.io.tmpdir=/disk/tmp",
                                "-Djava.awt.headless=true",
                                "-XX:+UseSerialGC",
                                "-Xmx4g",
                                "-jar", "%(AUX_JAR)s"])


def test_java_builder__defaults__kwargs():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar")
    assert_equal(builder.kwargs, {"AUX_JAR": "/path/Foo.jar",
                                  "CHECK_JRE": JAVA_VERSIONS[(1, 6)]})


def test_java_builder__multithreaded_gc():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar",
                                   temp_root="/disk/tmp",
                                   gc_threads=3)
    assert_equal(builder.call, ["java",
                                "-server",
                                "-Djava.io.tmpdir=/disk/tmp",
                                "-Djava.awt.headless=true",
                                "-XX:ParallelGCThreads=3",
                                "-Xmx4g",
                                "-jar", "%(AUX_JAR)s"])


def test_java_builder__multithreaded_gc__zero_or_negative_threads():
    assert_raises(ValueError, AtomicJavaCmdBuilder, "/path/Foo.jar", gc_threads=0)
    assert_raises(ValueError, AtomicJavaCmdBuilder, "/path/Foo.jar", gc_threads=-1)


def test_java_builder__multithreaded_gc__non_int_threads():
    assert_raises(TypeError, AtomicJavaCmdBuilder, "/path/Foo.jar", gc_threads="3")


def test_java_builder__kwargs():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar", set_cwd=True)
    assert_equal(builder.kwargs, {"AUX_JAR": "/path/Foo.jar",
                                  "set_cwd": True,
                                  "CHECK_JRE": JAVA_VERSIONS[(1, 6)]})


###############################################################################
###############################################################################
# AtomicMPICmdBuilder

def test_mpi_builder__defaults__str():
    builder = AtomicMPICmdBuilder("ls")
    assert_equal(builder.call, ["ls"])
    assert_equal(builder.kwargs, {"EXEC_MPI": "mpirun"})


def test_mpi_builder__multithreaded__str():
    builder = AtomicMPICmdBuilder("ls", threads=3)
    assert_equal(builder.call, ["mpirun", "-n", 3, "ls"])
    assert_equal(builder.kwargs, {"EXEC_MAIN": "ls"})


def test_mpi_builder__defaults__complex_cmd():
    builder = AtomicMPICmdBuilder(["python", "/foo/run.py"])
    assert_equal(builder.call, ["python", "/foo/run.py"])
    assert_equal(builder.kwargs, {"EXEC_MPI": "mpirun"})


def test_mpi_builder__multithreaded__complex_cmd():
    builder = AtomicMPICmdBuilder(["python", "/foo/run.py"], threads=3)
    assert_equal(builder.call, ["mpirun", "-n", 3, "python", "/foo/run.py"])
    assert_equal(builder.kwargs, {"EXEC_MAIN": "python"})


def test_mpi_builder__kwargs():
    builder = AtomicMPICmdBuilder("ls", set_cwd=True)
    assert_equal(builder.kwargs, {"set_cwd": True, "EXEC_MPI": "mpirun"})


def test_mpi_builder__threads__zero_or_negative():
    assert_raises(ValueError, AtomicMPICmdBuilder, "ls", threads=0)
    assert_raises(ValueError, AtomicMPICmdBuilder, "ls", threads=-1)


def test_mpi_builder__threads__non_int():
    assert_raises(TypeError, AtomicMPICmdBuilder, "ls", threads="3")


###############################################################################
###############################################################################
# create_customizable_cli_parameters

def test_custom_cli__single_named_arg():
    class SingleNamedArg:
        @create_customizable_cli_parameters
        def customize(cls, argument):
            return {}

    value = "A value"
    obj = SingleNamedArg.customize(value)
    assert_equal(obj.argument, value)


def test_custom_cli__adding_new_values():
    class SingleNamedArg:
        @create_customizable_cli_parameters
        def customize(cls):
            return {"dynamic": 12345}

    obj = SingleNamedArg.customize()
    assert_equal(obj.dynamic, 12345)


def test_custom_cli__multiple_named_args():
    class SingleNamedArg:
        @create_customizable_cli_parameters
        def customize(cls, first, second):
            return {}

    obj = SingleNamedArg.customize(123, 456)
    assert_equal(obj.first,  123)
    assert_equal(obj.second, 456)


def test_custom_cli__only_customize_is_valid_function_name():
    try:
        class ClassWithMisnamedFunction:
            @create_customizable_cli_parameters
            def not_called_customize(cls, first, second):
                return {}  # pragma: no coverage

        assert False, "ValueError not raised"  # pragma: no coverage
    except ValueError:
        pass


###############################################################################
###############################################################################
# use_customizable_cli_parameters


###############################################################################
###############################################################################
# apply_options

def test_apply_options__single_option__default_pred__set_when_pred_is_true():
    mock = flexmock()
    mock.should_receive('set_option').with_args('--foo', 17).once()
    apply_options(mock, {"--foo": 17})


def test_apply_options__single_option__default_pred__ignore_when_pred_is_false():
    mock = flexmock()
    apply_options(mock, {"Other": None})


def _user_pred(key):
    return key.startswith("FOO")


def test_apply_options__single_option__user_pred__set_when_pred_is_true():
    mock = flexmock()
    mock.should_receive('set_option').with_args('FOO_BAR', 17).once()
    apply_options(mock, {"FOO_BAR": 17}, _user_pred)


def test_apply_options__single_option__user_pred__ignore_when_pred_is_false():
    mock = flexmock()
    apply_options(mock, {"BAR_FOO": 17}, _user_pred)


def test_apply_options__single_option__boolean__set_when_value_is_true():
    mock = flexmock()
    mock.should_receive('set_option').with_args('-v')
    apply_options(mock, {"-v": True})


def test_apply_options__single_option__boolean__set_when_value_is_none():
    mock = flexmock()
    mock.should_receive('set_option').with_args('-v')
    apply_options(mock, {"-v": None})


def test_apply_options__single_option__boolean__pop_when_value_is_false():
    mock = flexmock()
    mock.should_receive('pop_option').with_args('-v')
    apply_options(mock, {"-v": False})


def test_apply_options__single_option__boolean__pop_missing_throws():
    mock = flexmock()
    mock.should_receive('pop_option').with_args('-v').and_raise(KeyError('-v'))
    assert_raises(KeyError, apply_options, mock, {"-v": False})


def test_apply_options__multiple_option():
    mock = flexmock()
    mock.should_receive('add_option').with_args('--foo', 3).once()
    mock.should_receive('add_option').with_args('--foo', 17).once()
    apply_options(mock, {"--foo": [3, 17]})


def test_apply_options__boolean_and_none_is_single_value_only():
    mock = flexmock()
    assert_raises(TypeError, apply_options, mock, {"--foo": [True]})
    assert_raises(TypeError, apply_options, mock, {"--foo": [False]})
    assert_raises(TypeError, apply_options, mock, {"--foo": [None]})


def test_apply_options__unexpected_types_in_values():
    mock = flexmock()
    assert_raises(TypeError, apply_options, mock, {"--foo": object()})
    assert_raises(TypeError, apply_options, mock, {"--foo": iter([])})
    assert_raises(TypeError, apply_options, mock, {"--foo": {}})
    assert_raises(TypeError, apply_options, mock, {"--foo": set()})


def test_apply_options__non_string_types_in_keys():
    mock = flexmock()
    assert_raises(TypeError, apply_options, mock, {1: 17})
    assert_raises(TypeError, apply_options, mock, {("foo",): 17})


def test_apply_options__not_dict_like():
    mock = flexmock()
    assert_raises(TypeError, apply_options, mock, None)
    assert_raises(TypeError, apply_options, mock, [1, 2, 3])
