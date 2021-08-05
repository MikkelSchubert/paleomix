import copy
import os.path
import string

from typing import Any, Dict, Iterator, Tuple


class LayoutError(Exception):
    pass


class Layout:
    def __init__(self, layout: Dict[str, Any], **kwargs) -> None:
        self.kwargs = kwargs

        self._layout = {}
        for key, value in self._flatten_layout(layout):
            self._layout[key] = value

        fmt = string.Formatter()
        self._fields = set()
        for value in self._layout.values():
            for _, name, _, _ in fmt.parse(value):
                if name == "":
                    raise LayoutError(f"unnamed field are not allowed in {value!r}")

                self._fields.add(name)

        self._validate_kwargs(self.kwargs)

    def update(self, **kwargs) -> "Layout":
        self._validate_kwargs(kwargs)

        layout = copy.deepcopy(self)
        layout.kwargs.update(kwargs)

        return layout

    def __iter__(self):
        return iter(self._layout)

    def __getitem__(self, key):
        return self._layout[key].format(**self.kwargs)

    def _validate_kwargs(self, kwargs):
        for key in kwargs:
            if key not in self._fields:
                raise LayoutError(f"unknown key {key!r}")

    @classmethod
    def _flatten_layout(cls, layout: Dict[str, Any]) -> Iterator[Tuple[str, str]]:
        for key, value in layout.items():
            if not isinstance(key, str):
                raise LayoutError(f"invalid key {key!r}")
            elif isinstance(value, str):
                yield value, key
            elif isinstance(value, dict):
                for label, path in cls._flatten_layout(value):
                    yield label, os.path.join(key, path)
            else:
                raise LayoutError(f"invalid value {value!r}")
