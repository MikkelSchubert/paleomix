# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2023 Mikkel Schubert <mikkelsch@gmail.com>
from __future__ import annotations

import argparse
from collections.abc import Sequence
from typing import Any

import configargparse

import paleomix

__all__ = [
    "SUPPRESS",
    "ArgumentDefaultsHelpFormatter",
    "ArgumentGroup",
    "ArgumentParser",
    "Namespace",
    "SubParsersAction",
]

Action = configargparse.Action
ArgumentGroup = argparse._ArgumentGroup  # pyright: ignore[reportPrivateUsage]  # noqa: SLF001
ArgumentParserBase = configargparse.ArgumentParser
Namespace = configargparse.Namespace
SubParsersAction = argparse._SubParsersAction  # pyright: ignore[reportPrivateUsage]  # noqa: SLF001
SUPPRESS = configargparse.SUPPRESS


class ArgumentDefaultsHelpFormatter(configargparse.ArgumentDefaultsHelpFormatter):
    """Modified ArgumentDefaultsHelpFormatter that excludes several constants (True,
    False, None) and uses a custom presentation of the default value.
    """

    def __init__(
        self,
        prog: str,
        *,
        indent_increment: int = 2,
        max_help_position: int = 24,
        width: int | None = 79,
    ) -> None:
        super().__init__(
            prog=prog,
            indent_increment=indent_increment,
            max_help_position=max_help_position,
            width=width,
        )

    def _get_help_string(self, action: configargparse.Action) -> str | None:
        # The following values look silly as part of a help string
        if isinstance(action.default, bool) or action.default in [None, [], ()]:
            return action.help

        # The subclass does not allow modification to the defaults string, so instead
        # we access the logic by simply checking if the result was modified.
        if super()._get_help_string(action) == action.help:
            return action.help

        assert action.help is not None
        return action.help + " [%(default)s]"


class ArgumentParser(ArgumentParserBase):
    """Supports keys with underscores instead of dashes, for backwards compatibility
    with old paleomix config files, provided that these do not use per-host sections.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("formatter_class", ArgumentDefaultsHelpFormatter)
        # Workaround for configargparse (1.2.3) not considering abbreviations when
        # applying options from config files, resulting in config file options
        # overriding abbreviated options supplied on the command-line.
        kwargs.setdefault("allow_abbrev", False)

        super().__init__(*args, **kwargs)

        self.add_argument(
            "-v",
            "--version",
            action="version",
            version="%(prog)s v" + paleomix.__version__,
        )

    def get_possible_config_keys(self, action: Action) -> list[str]:
        keys = super().get_possible_config_keys(action)
        for key in keys:
            key = key.strip("-").replace("-", "_")
            if key not in keys:
                keys.append(key)

        return keys

    def convert_item_to_command_line_arg(
        self,
        action: Action,
        key: str,
        value: str | list[str],
    ) -> list[str]:
        # Ignore empty options from old config files
        if action and value in ("", "="):
            return []

        return super().convert_item_to_command_line_arg(action, key, value)

    def add_subparsers(self, **kwargs: Any) -> SubParsersAction[ArgumentParser]:
        subparsers = super().add_subparsers(
            **kwargs,
            parser_class=ArgumentParser,
        )

        # Hack to hide aliases from subcommand help text, since aliases are only used
        # for deprecated commands/command-names
        subparsers._ChoicesPseudoAction = _ChoicesPseudoAction  # noqa: SLF001

        return subparsers


class _ChoicesPseudoAction(configargparse.Action):
    def __init__(self, name: str, _aliases: Sequence[str], help: str | None) -> None:
        super().__init__(option_strings=[], dest=name, help=help, metavar=name)
