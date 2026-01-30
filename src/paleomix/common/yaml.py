# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import warnings
from typing import IO

import ruamel.yaml
import ruamel.yaml.error
from ruamel.yaml.constructor import DuplicateKeyError
from ruamel.yaml.error import YAMLError

__all__ = [
    "DuplicateKeyError",
    "YAMLError",
    "safe_load",
]


def safe_load(stream: IO[str]) -> object:
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)  # pyright: ignore[reportGeneralTypeIssues]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ruamel.yaml.error.MantissaNoDotYAML1_1Warning)

        return yaml.load(stream)
