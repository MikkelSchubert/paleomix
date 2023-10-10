#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

import pytest

from paleomix.common.layout import Layout, LayoutError


def test_layout__minimal() -> None:
    layout = Layout({})

    assert layout.kwargs == {}
    assert list(layout) == []


def test_layout__simple_layout() -> None:
    layout = Layout({"{root}": "my_file"}, root="/root")

    assert layout.kwargs == {"root": "/root"}
    assert list(layout) == ["my_file"]
    assert layout["my_file"] == "/root"


def test_layout__nested_layout() -> None:
    layout = Layout(
        {
            "{root}": {
                "a": {"{name}.txt": "deep_file"},
                "b": "other_file",
            }
        },
        root="/root",
        name="foobar",
    )

    assert layout.kwargs == {"root": "/root", "name": "foobar"}
    assert sorted(layout) == ["deep_file", "other_file"]
    assert layout["deep_file"] == "/root/a/foobar.txt"
    assert layout["other_file"] == "/root/b"


def test_layout__update() -> None:
    layout1 = Layout({"{root}": {"{name}": "my_file"}})

    layout2 = layout1.update(root="/tmp")
    assert layout1.kwargs == {}
    assert layout2.kwargs == {"root": "/tmp"}

    layout3 = layout2.update(name="foobar")
    assert layout2.kwargs == {"root": "/tmp"}
    assert layout3.kwargs == {"root": "/tmp", "name": "foobar"}
    assert layout3["my_file"] == "/tmp/foobar"

    layout4 = layout3.update(root="/home/username")
    assert layout3.kwargs == {"root": "/tmp", "name": "foobar"}
    assert layout4.kwargs == {"root": "/home/username", "name": "foobar"}
    assert layout3["my_file"] == "/tmp/foobar"
    assert layout4["my_file"] == "/home/username/foobar"


def test_layout__get_without_updates() -> None:
    layout1 = Layout({"{root}": {"{name}": "my_file"}})

    assert layout1.get("my_file", root="foo", name="bar") == "foo/bar"
    assert layout1.kwargs == {}


def test_layout__get_without_updates__override() -> None:
    layout1 = Layout({"{root}": {"{name}": "my_file"}}, root="bar")

    assert layout1.get("my_file", root="foo", name="bar") == "foo/bar"
    assert layout1.kwargs == {"root": "bar"}


def test_layout__get_without_updates__partial() -> None:
    layout1 = Layout({"{root}": {"{name}": "my_file"}}, root="foo")

    assert layout1.get("my_file", name="bar") == "foo/bar"
    assert layout1.kwargs == {"root": "foo"}


def test_layout__extranous_path_components() -> None:
    layout = Layout({"{root}": {"{name}": "my_file"}})

    assert layout.get("my_file", root="foo/", name="test") == "foo/test"


def test_layout__unnamed_field() -> None:
    with pytest.raises(LayoutError, match="unnamed field are not allowed in"):
        Layout({"{}": "my_file"})


def test_layout__missing_value() -> None:
    layout = Layout({"{root}": "my_file"})

    with pytest.raises(KeyError, match="root"):
        layout["my_file"]


def test_layout__unknown_field__in_init() -> None:
    with pytest.raises(LayoutError, match="unknown key"):
        Layout({"{root}": "my_file"}, foobar="/path/to/somewhere")


def test_layout__unknown_field__in_update() -> None:
    layout = Layout({"{root}": "my_file"})

    with pytest.raises(LayoutError, match="unknown key"):
        layout.update(foobar="/path/to/somewhere")


def test_layout__non_string_name() -> None:
    with pytest.raises(LayoutError, match="invalid key 17"):
        Layout({17: "my_file"})  # pyright: ignore[reportGeneralTypeIssues]


def test_layout__non_string_dict_value() -> None:
    with pytest.raises(LayoutError, match="invalid value 17"):
        Layout({"{root}": 17})  # pyright: ignore[reportGeneralTypeIssues]


def test_layout__duplicate_path_names() -> None:
    with pytest.raises(LayoutError, match="'file_1' used multiple times"):
        Layout({"foo": "file_1", "{root}": {"zod": "file_1"}})


def test_layout__non_unique_key_field_name() -> None:
    with pytest.raises(LayoutError, match="'file_1' used as both key and field"):
        Layout({"foo": "file_1", "{root}": {"{file_1}": "file_2"}})


def test_layout__get_field() -> None:
    layout = Layout({"{root}": {"b": "other_file"}}, root="/root")

    assert layout.get_field("root") == "/root"


def test_layout__get_missing_field() -> None:
    layout = Layout({"{root}": {"b": "other_file"}}, root="/root")

    with pytest.raises(KeyError):
        assert layout.get_field("unknown_key")
