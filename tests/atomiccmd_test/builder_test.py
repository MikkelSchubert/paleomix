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
from unittest.mock import call, Mock, patch

import pytest


from paleomix.atomiccmd.builder import (
    AtomicCmdBuilder,
    AtomicCmdBuilderError,
    AtomicJavaCmdBuilder,
    AtomicMPICmdBuilder,
    apply_options,
)


###############################################################################
###############################################################################
# AtomicCmdBuilder: Constructor


def test_builder__simple__call():
    builder = AtomicCmdBuilder(["ls"])
    assert builder.call == ["ls"]


def test_builder__simple__str():
    builder = AtomicCmdBuilder("ls")
    assert builder.call == ["ls"]


def test_builder__simple__iterable():
    builder = AtomicCmdBuilder(iter(["ls"]))
    assert builder.call == ["ls"]


def test_builder__complex():
    builder = AtomicCmdBuilder(["java", "jar", "/a/jar"])
    assert builder.call == ["java", "jar", "/a/jar"]


def test_builder__kwargs__empty():
    builder = AtomicCmdBuilder(["ls"])
    assert builder.kwargs == {}


def test_builder__kwargs():
    expected = {"IN_FILE": "/abc/def.txt", "OUT_FILE": "/etc/fstab"}
    builder = AtomicCmdBuilder(["ls"], **expected)
    assert builder.kwargs == expected


def test_builder__kwargs__set_cwd():
    builder = AtomicCmdBuilder(["ls"], set_cwd=True)
    assert builder.kwargs == {"set_cwd": True}


###############################################################################
###############################################################################
# AtomicCmdBuilder: set_option


def test_builder__set_option():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt")
    assert builder.call == ["find", "-name", "*.txt"]


def test_builder__set_option__overwrite():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt", fixed=False)
    builder.set_option("-name", "*.bat")
    assert builder.call == ["find", "-name", "*.bat"]


def test_builder__set_option__overwrite_fixed():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt")
    with pytest.raises(AtomicCmdBuilderError):
        builder.set_option("-name", "*.bat")


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_option


def test_builder__add_option():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-name", "*.txt")
    assert builder.call == ["find", "-name", "*.txt"]


def test_builder__add_option__overwrite():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-name", "*.txt")
    builder.add_option("-or")
    builder.add_option("-name", "*.bat")
    assert builder.call == ["find", "-name", "*.txt", "-or", "-name", "*.bat"]


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_option / set_option common tests

_ADD_SET_OPTION = [AtomicCmdBuilder.add_option, AtomicCmdBuilder.set_option]


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__add_or_set_option__without_value(setter):
    builder = AtomicCmdBuilder("find")
    setter(builder, "-delete")
    assert builder.call == ["find", "-delete"]


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__add_or_set_option__with_sep(setter):
    builder = AtomicCmdBuilder("find")
    setter(builder, "-size", "0", sep="=")
    assert builder.call == ["find", "-size=0"]


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__add_or_set_option__with_non_str_value(setter):
    builder = AtomicCmdBuilder("find")
    setter(builder, "-size", 0)
    assert builder.call == ["find", "-size", 0]


def test_builder__add_or_set_option__add_and_set():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-name", "*.txt")
    with pytest.raises(AtomicCmdBuilderError):
        AtomicCmdBuilder.add_option(builder, "-name", "*.bat")


def test_builder__add_or_set_option__set_and_add():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-name", "*.txt")
    with pytest.raises(AtomicCmdBuilderError):
        AtomicCmdBuilder.set_option(builder, "-name", "*.bat")


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__add_or_set_option__with_non_str_key(setter):
    builder = AtomicCmdBuilder("find")
    with pytest.raises(TypeError):
        setter(builder, 7913, "True")


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__add_or_set_option__after_finalize(setter):
    builder = AtomicCmdBuilder("find")
    builder.finalize()
    with pytest.raises(AtomicCmdBuilderError):
        setter(builder, "-size", "1")


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__add_or_set_option__empty_key(setter):
    builder = AtomicCmdBuilder("find")
    with pytest.raises(KeyError):
        setter(builder, "", "1")


###############################################################################
###############################################################################
# AtomicCmdBuilder: pop_option


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__pop_option(setter):
    builder = AtomicCmdBuilder("find")
    setter(builder, "-size", "0", fixed=False)
    builder.pop_option("-size")
    assert builder.call == ["find"]


def test_builder__pop_option__last_option():
    builder = AtomicCmdBuilder("find")
    builder.add_option("-size", "0", fixed=False)
    builder.add_option("-size", "1", fixed=False)
    builder.pop_option("-size")
    assert builder.call == ["find", "-size", "0"]


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__pop_option__different_options(setter):
    builder = AtomicCmdBuilder("find")
    setter(builder, "-empty", fixed=False)
    setter(builder, "-size", "1", fixed=False)
    setter(builder, "-name", "*.txt", fixed=False)
    builder.pop_option("-size")
    assert builder.call == ["find", "-empty", "-name", "*.txt"]


@pytest.mark.parametrize("setter", _ADD_SET_OPTION)
def test_builder__pop_option__is_fixed(setter):
    builder = AtomicCmdBuilder("find")
    setter(builder, "-size", "0")
    with pytest.raises(AtomicCmdBuilderError):
        builder.pop_option("-size")


def test_builder__pop_option__empty():
    builder = AtomicCmdBuilder("find")
    with pytest.raises(KeyError):
        builder.pop_option("-size")


def test_builder__pop_option__missing_key():
    builder = AtomicCmdBuilder("find")
    builder.set_option("-size", 0)
    with pytest.raises(KeyError):
        builder.pop_option("-isize")


def test_builder__pop_option__with_non_str_key():
    builder = AtomicCmdBuilder("find")
    with pytest.raises(TypeError):
        builder.pop_option(7913)


def test_builder__pop_option__with_empty_key():
    builder = AtomicCmdBuilder("find")
    with pytest.raises(KeyError):
        builder.pop_option("")


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_value


def test_builder__add_value():
    builder = AtomicCmdBuilder("ls")
    builder.add_value("%(IN_FILE)s")
    assert builder.call == ["ls", "%(IN_FILE)s"]


def test_builder__add_value__two_values():
    builder = AtomicCmdBuilder("ls")
    builder.add_value("%(IN_FILE)s")
    builder.add_value("%(OUT_FILE)s")
    assert builder.call == ["ls", "%(IN_FILE)s", "%(OUT_FILE)s"]


###############################################################################
###############################################################################
# AtomicCmdBuilder: set_kwargs


def test_builder__set_kwargs__called_once():
    expected = {"IN_PATH": "/a/b/", "OUT_PATH": "/dst/file"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(**expected)
    assert builder.kwargs == expected


def test_builder__set_kwargs__called_twice():
    expected = {"IN_PATH": "/a/b/", "OUT_PATH": "/dst/file"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(OUT_PATH="/dst/file")
    builder.set_kwargs(IN_PATH="/a/b/")
    assert builder.kwargs == expected


def test_builder__set_kwargs__atomiccmdbuilder():
    mock = Mock(AtomicCmdBuilder("true"))
    mock.finalize.return_value = "finalized!"
    builder = AtomicCmdBuilder("ls", IN_BUILDER=mock)
    assert builder.kwargs == {"IN_BUILDER": "finalized!"}


def test_builder__set_kwargs__after_finalize():
    expected = {"IN_PATH": "/a/b/"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(IN_PATH="/a/b/")
    builder.finalize()
    with pytest.raises(AtomicCmdBuilderError):
        builder.set_kwargs(OUT_PATH="/dst/file")
    assert builder.kwargs == expected


def test_builder__set__kwargs__overwriting():
    expected = {"IN_PATH": "/a/b/"}
    builder = AtomicCmdBuilder("echo")
    builder.set_kwargs(IN_PATH="/a/b/")
    with pytest.raises(AtomicCmdBuilderError):
        builder.set_kwargs(IN_PATH="/dst/file")
    assert builder.kwargs == expected


###############################################################################
###############################################################################
# AtomicCmdBuilder: finalize


def test_builder__finalized_call__simple_command():
    builder = AtomicCmdBuilder("echo")
    assert builder.finalized_call == ["echo"]


def test_builder__finalized_call__kwargs_are_instantiated():
    builder = AtomicCmdBuilder(
        ("echo", "%(ARG1)s", "X=%(ARG2)s"), ARG1="/foo/bar", ARG2="zod"
    )
    assert builder.finalized_call == ["echo", "/foo/bar", "X=zod"]


def test_builder__finalized_call__kwargs_are_instantiated__with_temp_dir():
    builder = AtomicCmdBuilder(("echo", "%(ARG)s", "%(TEMP_DIR)s"), ARG="/foo/bar")
    assert builder.finalized_call == ["echo", "/foo/bar", "%(TEMP_DIR)"]


def test_builder__finalized_call__kwargs_are_instantiated__with_non_str_arg():
    builder = AtomicCmdBuilder(("echo", "%(ARG)s", 17), ARG="/foo/bar")
    assert builder.finalized_call == ["echo", "/foo/bar", "17"]


###############################################################################
###############################################################################
# AtomicCmdBuilder: finalize


def test_builder__finalize__returns_singleton():
    builder = AtomicCmdBuilder("echo")
    assert builder.finalize() is builder.finalize()


def test_builder__finalize__calls_atomiccmd():
    builder = AtomicCmdBuilder("echo", set_cwd=True)
    builder.add_option("-out", "%(OUT_FILE)s")
    builder.add_value("%(IN_FILE)s")
    builder.set_kwargs(OUT_FILE="/out/file", IN_FILE="/in/file")
    with patch("paleomix.atomiccmd.builder.AtomicCmd") as mock:
        builder.finalize()

    assert mock.mock_calls == [
        call(
            ["echo", "-out", "%(OUT_FILE)s", "%(IN_FILE)s"],
            IN_FILE="/in/file",
            OUT_FILE="/out/file",
            set_cwd=True,
        )
    ]


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_multiple_options


def test_builder__add_multiple_options():
    values = ("file_a", "file_b")
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", values)

    assert kwargs == expected
    assert builder.kwargs == expected
    assert builder.call == ["ls", "-i", "%(IN_FILE_01)s", "-i", "%(IN_FILE_02)s"]


def test_builder__add_multiple_options_with_sep():
    values = ("file_a", "file_b")
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", values, sep="=")

    assert kwargs == expected
    assert builder.kwargs == expected
    assert builder.call == ["ls", "-i=%(IN_FILE_01)s", "-i=%(IN_FILE_02)s"]


def test_builder__add_multiple_options_with_template():
    values = ("file_a", "file_b")
    expected = {"OUT_BAM_1": "file_a", "OUT_BAM_2": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", values, template="OUT_BAM_%i")

    assert kwargs == expected
    assert builder.kwargs == expected
    assert builder.call == ["ls", "-i", "%(OUT_BAM_1)s", "-i", "%(OUT_BAM_2)s"]


def test_builder__add_multiple_options_multiple_times():
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_options("-i", ("file_a",))
    assert kwargs == {"IN_FILE_01": "file_a"}
    kwargs = builder.add_multiple_options("-i", ("file_b",))
    assert kwargs == {"IN_FILE_02": "file_b"}

    assert builder.kwargs == expected
    assert builder.call == ["ls", "-i", "%(IN_FILE_01)s", "-i", "%(IN_FILE_02)s"]


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_multiple_values


def test_builder__add_multiple_values():
    values = ("file_a", "file_b")
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_values(values)

    assert kwargs == expected
    assert builder.kwargs == expected
    assert builder.call == ["ls", "%(IN_FILE_01)s", "%(IN_FILE_02)s"]


def test_builder__add_multiple_values_with_template():
    values = ("file_a", "file_b")
    expected = {"OUT_BAM_1": "file_a", "OUT_BAM_2": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_values(values, template="OUT_BAM_%i")

    assert kwargs == expected
    assert builder.kwargs == expected
    assert builder.call == ["ls", "%(OUT_BAM_1)s", "%(OUT_BAM_2)s"]


def test_builder__add_multiple_values_multiple_times():
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_values(("file_a",))
    assert kwargs == {"IN_FILE_01": "file_a"}
    kwargs = builder.add_multiple_values(("file_b",))
    assert kwargs == {"IN_FILE_02": "file_b"}

    assert builder.kwargs == expected
    assert builder.call == ["ls", "%(IN_FILE_01)s", "%(IN_FILE_02)s"]


###############################################################################
###############################################################################
# AtomicCmdBuilder: add_multiple_kwargs


def test_builder__add_multiple_kwargs():
    values = ("file_a", "file_b")
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_kwargs(values)

    assert kwargs == expected
    assert builder.kwargs == expected
    assert builder.call == ["ls"]


def test_builder__add_multiple_kwargs_with_template():
    values = ("file_a", "file_b")
    expected = {"OUT_BAM_1": "file_a", "OUT_BAM_2": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_kwargs(values, template="OUT_BAM_%i")

    assert kwargs == expected
    assert builder.kwargs == expected
    assert builder.call == ["ls"]


def test_builder__add_multiple_kwargs_multiple_times():
    expected = {"IN_FILE_01": "file_a", "IN_FILE_02": "file_b"}

    builder = AtomicCmdBuilder("ls")
    kwargs = builder.add_multiple_kwargs(("file_a",))
    assert kwargs == {"IN_FILE_01": "file_a"}
    kwargs = builder.add_multiple_kwargs(("file_b",))
    assert kwargs == {"IN_FILE_02": "file_b"}

    assert builder.kwargs == expected
    assert builder.call == ["ls"]


###############################################################################
###############################################################################
# AtomicJavaCmdBuilder


def test_java_builder__default__no_config():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar")
    assert builder.call == [
        "java",
        "-server",
        "-Djava.io.tmpdir=%(TEMP_DIR)s",
        "-Djava.awt.headless=true",
        "-XX:+UseSerialGC",
        "-Xmx4g",
        "-jar",
        "%(AUX_JAR)s",
    ]


def test_java_builder__defaults__call():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar", temp_root="/disk/tmp")
    assert builder.call == [
        "java",
        "-server",
        "-Djava.io.tmpdir=/disk/tmp",
        "-Djava.awt.headless=true",
        "-XX:+UseSerialGC",
        "-Xmx4g",
        "-jar",
        "%(AUX_JAR)s",
    ]


def test_java_builder__defaults__kwargs():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar")
    assert builder.kwargs == {"AUX_JAR": "/path/Foo.jar"}


def test_java_builder__multithreaded_gc():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar", temp_root="/disk/tmp", gc_threads=3)
    assert builder.call == [
        "java",
        "-server",
        "-Djava.io.tmpdir=/disk/tmp",
        "-Djava.awt.headless=true",
        "-XX:ParallelGCThreads=3",
        "-Xmx4g",
        "-jar",
        "%(AUX_JAR)s",
    ]


def test_java_builder__multithreaded_gc__zero_or_negative_threads():
    with pytest.raises(ValueError):
        AtomicJavaCmdBuilder("/path/Foo.jar", gc_threads=0)
    with pytest.raises(ValueError):
        AtomicJavaCmdBuilder("/path/Foo.jar", gc_threads=-1)


def test_java_builder__multithreaded_gc__non_int_threads():
    with pytest.raises(TypeError):
        AtomicJavaCmdBuilder("/path/Foo.jar", gc_threads="3")


def test_java_builder__kwargs():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar", set_cwd=True)
    assert builder.kwargs == {"AUX_JAR": "/path/Foo.jar", "set_cwd": True}


def test_java_builder__default__jre_options():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar", jre_options=("-XtraOption",))
    assert builder.call == [
        "java",
        "-server",
        "-Djava.io.tmpdir=%(TEMP_DIR)s",
        "-Djava.awt.headless=true",
        "-XX:+UseSerialGC",
        "-XtraOption",
        "-Xmx4g",
        "-jar",
        "%(AUX_JAR)s",
    ]


def test_java_builder__default__jre_options__Xmx():
    builder = AtomicJavaCmdBuilder("/path/Foo.jar", jre_options=("-Xmx8g",))
    assert builder.call == [
        "java",
        "-server",
        "-Djava.io.tmpdir=%(TEMP_DIR)s",
        "-Djava.awt.headless=true",
        "-XX:+UseSerialGC",
        "-Xmx8g",
        "-jar",
        "%(AUX_JAR)s",
    ]


###############################################################################
###############################################################################
# AtomicMPICmdBuilder


def test_mpi_builder__defaults__str():
    builder = AtomicMPICmdBuilder("ls")
    assert builder.call == ["ls"]
    assert builder.kwargs == {"EXEC_MPI": "mpirun"}


def test_mpi_builder__multithreaded__str():
    builder = AtomicMPICmdBuilder("ls", threads=3)
    assert builder.call == ["mpirun", "-n", 3, "ls"]
    assert builder.kwargs == {"EXEC_MAIN": "ls"}


def test_mpi_builder__defaults__complex_cmd():
    builder = AtomicMPICmdBuilder(["python", "/foo/run.py"])
    assert builder.call == ["python", "/foo/run.py"]
    assert builder.kwargs == {"EXEC_MPI": "mpirun"}


def test_mpi_builder__multithreaded__complex_cmd():
    builder = AtomicMPICmdBuilder(["python", "/foo/run.py"], threads=3)
    assert builder.call == ["mpirun", "-n", 3, "python", "/foo/run.py"]
    assert builder.kwargs == {"EXEC_MAIN": "python"}


def test_mpi_builder__kwargs():
    builder = AtomicMPICmdBuilder("ls", set_cwd=True)
    assert builder.kwargs == {"set_cwd": True, "EXEC_MPI": "mpirun"}


def test_mpi_builder__threads__zero_or_negative():
    with pytest.raises(ValueError):
        AtomicMPICmdBuilder("ls", threads=0)
    with pytest.raises(ValueError):
        AtomicMPICmdBuilder("ls", threads=-1)


def test_mpi_builder__threads__non_int():
    with pytest.raises(TypeError):
        AtomicMPICmdBuilder("ls", threads="3")


###############################################################################
###############################################################################
# apply_options


def test_apply_options__single_option__default_pred__set_when_pred_is_true():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"--foo": 17})
    mock.set_option.assert_called_once_with("--foo", 17)


def test_apply_options__single_option__default_pred__ignore_false_pred():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"Other": None})


def _user_pred(key):
    return key.startswith("FOO")


def test_apply_options__single_option__user_pred__set_when_pred_is_true():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"FOO_BAR": 17}, _user_pred)
    mock.set_option.assert_called_once_with("FOO_BAR", 17)


def test_apply_options__single_option__user_pred__ignore_when_pred_is_false():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"BAR_FOO": 17}, _user_pred)


def test_apply_options__single_option__boolean__set_when_value_is_true():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"-v": True})
    mock.set_option.assert_called_once_with("-v")


def test_apply_options__single_option__boolean__set_when_value_is_none():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"-v": None})
    mock.set_option.assert_called_once_with("-v")


def test_apply_options__single_option__boolean__pop_when_value_is_false():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"-v": False})
    mock.pop_option.assert_called_once_with("-v")


def test_apply_options__single_option__boolean__pop_missing_throws():
    mock = Mock(AtomicCmdBuilder("ls"))
    mock.pop_option.side_effect = KeyError("-v")
    with pytest.raises(KeyError):
        apply_options(mock, {"-v": False})
    mock.pop_option.assert_called_once_with("-v")


def test_apply_options__multiple_option():
    mock = Mock(AtomicCmdBuilder("ls"))
    apply_options(mock, {"--foo": [3, 17]})
    assert mock.mock_calls == [
        call.add_option("--foo", 3),
        call.add_option("--foo", 17),
    ]


def test_apply_options__boolean_and_none_is_single_value_only():
    mock = Mock(AtomicCmdBuilder("ls"))
    with pytest.raises(TypeError):
        apply_options(mock, {"--foo": [True]})
    with pytest.raises(TypeError):
        apply_options(mock, {"--foo": [False]})
    with pytest.raises(TypeError):
        apply_options(mock, {"--foo": [None]})


def test_apply_options__unexpected_types_in_values():
    mock = Mock(AtomicCmdBuilder("ls"))
    with pytest.raises(TypeError):
        apply_options(mock, {"--foo": object()})
    with pytest.raises(TypeError):
        apply_options(mock, {"--foo": iter([])})
    with pytest.raises(TypeError):
        apply_options(mock, {"--foo": {}})
    with pytest.raises(TypeError):
        apply_options(mock, {"--foo": set()})


def test_apply_options__non_string_types_in_keys():
    mock = Mock(AtomicCmdBuilder("ls"))
    with pytest.raises(TypeError):
        apply_options(mock, {1: 17})
    with pytest.raises(TypeError):
        apply_options(mock, {("foo",): 17})


def test_apply_options__not_dict_like():
    mock = Mock(AtomicCmdBuilder("ls"))
    with pytest.raises(TypeError):
        apply_options(mock, None)
    with pytest.raises(TypeError):
        apply_options(mock, [1, 2, 3])
