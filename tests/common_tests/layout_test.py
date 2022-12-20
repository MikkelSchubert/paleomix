import pytest

from paleomix.common.layout import Layout, LayoutError


def test_layout__minimal():
    layout = Layout({})

    assert layout.kwargs == {}
    assert list(layout) == []


def test_layout__simple_layout():
    layout = Layout({"{root}": "my_file"}, root="/root")

    assert layout.kwargs == {"root": "/root"}
    assert list(layout) == ["my_file"]
    assert layout["my_file"] == "/root"


def test_layout__nested_layout():
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


def test_layout__update():
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


def test_layout__get_without_updates():
    layout1 = Layout({"{root}": {"{name}": "my_file"}})

    assert layout1.get("my_file", root="foo", name="bar") == "foo/bar"
    assert layout1.kwargs == {}


def test_layout__get_without_updates__override():
    layout1 = Layout({"{root}": {"{name}": "my_file"}}, root="bar")

    assert layout1.get("my_file", root="foo", name="bar") == "foo/bar"
    assert layout1.kwargs == {"root": "bar"}


def test_layout__get_without_updates__partial():
    layout1 = Layout({"{root}": {"{name}": "my_file"}}, root="foo")

    assert layout1.get("my_file", name="bar") == "foo/bar"
    assert layout1.kwargs == {"root": "foo"}


def test_layout__extranous_path_components():
    layout = Layout({"{root}": {"{name}": "my_file"}})

    assert layout.get("my_file", root="foo/", name="test") == "foo/test"


def test_layout__unnamed_field():
    with pytest.raises(LayoutError, match="unnamed field are not allowed in"):
        Layout({"{}": "my_file"})


def test_layout__missing_value():
    layout = Layout({"{root}": "my_file"})

    with pytest.raises(KeyError, match="root"):
        layout["my_file"]


def test_layout__unknown_field__in_init():
    with pytest.raises(LayoutError, match="unknown key"):
        Layout({"{root}": "my_file"}, foobar="/path/to/somewhere")


def test_layout__unknown_field__in_update():
    layout = Layout({"{root}": "my_file"})

    with pytest.raises(LayoutError, match="unknown key"):
        layout.update(foobar="/path/to/somewhere")


def test_layout__non_string_name():
    with pytest.raises(LayoutError, match="invalid key 17"):
        Layout({17: "my_file"})  # type: ignore


def test_layout__non_string_dict_value():
    with pytest.raises(LayoutError, match="invalid value 17"):
        Layout({"{root}": 17})  # type: ignore


def test_layout__duplicate_path_names():
    with pytest.raises(LayoutError, match="'file_1' used multiple times"):
        Layout({"foo": "file_1", "{root}": {"zod": "file_1"}})


def test_layout__non_unique_key_field_name():
    with pytest.raises(LayoutError, match="'file_1' used as both key and field"):
        Layout({"foo": "file_1", "{root}": {"{file_1}": "file_2"}})


def test_layout__get_field():
    layout = Layout({"{root}": {"b": "other_file"}}, root="/root")

    assert layout.get_field("root") == "/root"


def test_layout__get_missing_field():
    layout = Layout({"{root}": {"b": "other_file"}}, root="/root")

    with pytest.raises(KeyError):
        assert layout.get_field("unknown_key")
