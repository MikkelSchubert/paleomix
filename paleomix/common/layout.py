import copy
import os.path
import string
from typing import Any, Dict, Iterator, Tuple


class LayoutError(Exception):
    pass


class Layout:
    def __init__(self, layout: Dict[str, Any], **kwargs) -> None:
        self.kwargs = kwargs

        self._layout: Dict[str, Tuple[str, ...]] = {}
        for key, value in self._flatten_layout(layout):
            if key in self._layout:
                raise LayoutError(f"key {key!r} used multiple times")

            self._layout[key] = value

        self._fields = self._collect_fields(self._layout)

        self._validate_kwargs(self.kwargs)

    def get(self, key, **kwargs):
        self._validate_kwargs(kwargs)

        merged_kwargs = dict(self.kwargs)
        merged_kwargs.update(kwargs)

        return self._build_path(key, merged_kwargs)

    def get_field(self, key):
        return self.kwargs[key]

    def update(self, **kwargs) -> "Layout":
        self._validate_kwargs(kwargs)

        layout = copy.deepcopy(self)
        layout.kwargs.update(kwargs)

        return layout

    def __iter__(self):
        return iter(self._layout)

    def __getitem__(self, key):
        return self._build_path(key, self.kwargs)

    def _validate_kwargs(self, kwargs):
        for key in kwargs:
            if key not in self._fields:
                raise LayoutError(f"unknown key {key!r}")

    def _build_path(self, key, kwargs):
        components = []
        for component in self._layout[key]:
            components.append(component.format(**kwargs))

        return os.path.join(*components)

    @classmethod
    def _flatten_layout(
        cls, layout: Dict[str, Any]
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
    def _collect_fields(layout):
        fmt = string.Formatter()
        fields = set()
        for path in layout.values():
            for value in path:
                for _, name, _, _ in fmt.parse(value):
                    if name == "":
                        raise LayoutError(f"unnamed field are not allowed in {value!r}")
                    elif name in layout:
                        raise LayoutError(f"{name!r} used as both key and field name")

                    fields.add(name)

        return fields
