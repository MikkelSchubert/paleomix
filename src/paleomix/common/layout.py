# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import copy
import os.path
import string
from collections.abc import Iterator
from typing import Union

from typing_extensions import TypeAlias

LayoutType: TypeAlias = dict[str, Union[str, "LayoutType"]]


class LayoutError(Exception):
    pass


class Layout:
    kwargs: dict[str, str | int]
    _layout: dict[str, tuple[str, ...]]
    _fields: set[str]

    def __init__(self, layout: LayoutType, **kwargs: str | int) -> None:
        self.kwargs = kwargs
        self._layout = {}
        for key, value in self._flatten_layout(layout):
            if key in self._layout:
                raise LayoutError(f"key {key!r} used multiple times")

            self._layout[key] = value

        self._fields = self._collect_fields(self._layout)

        self._validate_kwargs(self.kwargs)

    def get(self, key: str, **kwargs: str | int) -> str:
        self._validate_kwargs(kwargs)

        merged_kwargs = dict(self.kwargs)
        merged_kwargs.update(kwargs)

        return self._build_path(key, merged_kwargs)

    def get_field(self, key: str) -> str | int:
        return self.kwargs[key]

    def update(self, **kwargs: str | int) -> Layout:
        self._validate_kwargs(kwargs)

        layout = copy.deepcopy(self)
        layout.kwargs.update(kwargs)

        return layout

    def __iter__(self) -> Iterator[str]:
        return iter(self._layout)

    def __getitem__(self, key: str) -> str:
        return self._build_path(key, self.kwargs)

    def _validate_kwargs(self, kwargs: dict[str, str | int]) -> None:
        for key in kwargs:
            if key not in self._fields:
                raise LayoutError(f"unknown key {key!r}")

    def _build_path(self, key: str, kwargs: dict[str, str | int]) -> str:
        try:
            path = self._layout[key]
        except KeyError:
            raise KeyError(f"unknown layout key {key!r}") from None
        try:
            return os.path.join(*(it.format(**kwargs) for it in path))
        except KeyError as error:
            raise KeyError(f"value {error} missing for layout key {key!r}") from None

    @classmethod
    def _flatten_layout(
        cls,
        layout: LayoutType,
    ) -> Iterator[tuple[str, tuple[str, ...]]]:
        for key, value in layout.items():
            if not isinstance(key, str):
                raise LayoutError(f"invalid key {key!r}")
            elif isinstance(value, str):
                yield value, (key,)
            elif isinstance(value, dict):
                for label, path in cls._flatten_layout(value):
                    yield label, (key, *path)
            else:
                raise LayoutError(f"invalid value {value!r}")

    @staticmethod
    def _collect_fields(layout: dict[str, tuple[str, ...]]) -> set[str]:
        fmt = string.Formatter()
        fields: set[str] = set()
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
