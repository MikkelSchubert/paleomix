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

import copy
import os.path
import string
from typing import Dict, Iterator, List, Set, Tuple, Union

LayoutType = Dict[str, Union[str, "LayoutType"]]


class LayoutError(Exception):
    pass


class Layout:
    kwargs: Dict[str, str]
    _layout: Dict[str, Tuple[str, ...]]
    _fields: Set[str]

    def __init__(self, layout: LayoutType, **kwargs: str) -> None:
        self.kwargs = kwargs
        self._layout = {}
        for key, value in self._flatten_layout(layout):
            if key in self._layout:
                raise LayoutError(f"key {key!r} used multiple times")

            self._layout[key] = value

        self._fields = self._collect_fields(self._layout)

        self._validate_kwargs(self.kwargs)

    def get(self, key: str, **kwargs: str) -> str:
        self._validate_kwargs(kwargs)

        merged_kwargs = dict(self.kwargs)
        merged_kwargs.update(kwargs)

        return self._build_path(key, merged_kwargs)

    def get_field(self, key: str) -> str:
        return self.kwargs[key]

    def update(self, **kwargs: str) -> "Layout":
        self._validate_kwargs(kwargs)

        layout = copy.deepcopy(self)
        layout.kwargs.update(kwargs)

        return layout

    def __iter__(self) -> Iterator[str]:
        return iter(self._layout)

    def __getitem__(self, key: str) -> str:
        return self._build_path(key, self.kwargs)

    def _validate_kwargs(self, kwargs: Dict[str, str]) -> None:
        for key in kwargs:
            if key not in self._fields:
                raise LayoutError(f"unknown key {key!r}")

    def _build_path(self, key: str, kwargs: Dict[str, str]) -> str:
        components: List[str] = []
        for component in self._layout[key]:
            components.append(component.format(**kwargs))

        return os.path.join(*components)

    @classmethod
    def _flatten_layout(
        cls,
        layout: LayoutType,
    ) -> Iterator[Tuple[str, Tuple[str, ...]]]:
        for key, value in layout.items():
            if not isinstance(key, str):
                raise LayoutError(f"invalid key {key!r}")
            elif isinstance(value, str):
                yield value, (key,)
            elif isinstance(value, dict):
                for label, path in cls._flatten_layout(value):
                    yield label, (key,) + path
            else:
                raise LayoutError(f"invalid value {value!r}")

    @staticmethod
    def _collect_fields(layout: Dict[str, Tuple[str, ...]]) -> Set[str]:
        fmt = string.Formatter()
        fields: Set[str] = set()
        for path in layout.values():
            for value in path:
                for _, name, _, _ in fmt.parse(value):
                    if name == "":
                        raise LayoutError(f"unnamed field are not allowed in {value!r}")
                    elif name in layout:
                        raise LayoutError(f"{name!r} used as both key and field name")
                    elif name is not None:
                        fields.add(name)

        return fields
