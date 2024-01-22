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

import argparse
from typing import Any

import configargparse

import paleomix

__all__ = [
    "ArgumentDefaultsHelpFormatter",
    "ArgumentGroup",
    "ArgumentParser",
    "Namespace",
    "SubParsersAction",
    "SUPPRESS",
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

    def get_possible_config_keys(self, *args: Any, **kwargs: Any) -> list[str]:
        keys = super().get_possible_config_keys(*args, **kwargs)
        for key in keys:
            key = key.strip("-").replace("-", "_")
            if key not in keys:
                keys.append(key)

        return keys

    def convert_item_to_command_line_arg(self, action, key, value) -> list[str]:
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
    def __init__(self, name, aliases, help):  # noqa: A002
        super().__init__(
            option_strings=[],
            dest=name,
            help=help,
            metavar=name,
        )
