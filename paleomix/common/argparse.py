#!/usr/bin/env python3
from typing import Any, List

import configargparse
import paleomix

SUPPRESS = configargparse.SUPPRESS


class ArgumentDefaultsHelpFormatter(configargparse.ArgumentDefaultsHelpFormatter):
    """Modified ArgumentDefaultsHelpFormatter that excludes several constants (True,
    False, None) and uses a custom presentation of the default value.
    """

    def __init__(self, *args: Any, **kwargs: Any):
        # Enable wordwrapping
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)

    def _get_help_string(self, action: configargparse.Action):
        # The following values look silly as part of a help string
        if isinstance(action.default, bool) or action.default in [None, [], ()]:
            return action.help

        # The subclass does not allow modification to the defaults string, so instead
        # we access the logic by simply checking if the result was modified.
        if super()._get_help_string(action) == action.help:
            return action.help

        assert action.help is not None
        return action.help + " [%(default)s]"


class ArgumentParser(configargparse.ArgumentParser):
    """Supports keys with underscores instead of dashes, for backwards compatibility
    with old paleomix config files, provided that these do not use per-host setions.
    """

    def __init__(self, *args: Any, **kwargs: Any):
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

    def get_possible_config_keys(self, *args: Any, **kwargs: Any) -> List[str]:
        keys = super().get_possible_config_keys(*args, **kwargs)
        for key in keys:
            key = key.strip("-").replace("-", "_")
            if key not in keys:
                keys.append(key)

        return keys

    def convert_item_to_command_line_arg(self, action, key, value):
        # Ignore empty options from old config files
        if action and value in ("", "="):
            return []

        return super().convert_item_to_command_line_arg(action, key, value)

    def add_subparsers(self, *args: Any, **kwargs: Any):
        subparsers = super().add_subparsers(*args, **kwargs)
        # Hack to hide aliases from subcommand help text, since aliases are only used
        # for deprecated commands/command-names
        subparsers._ChoicesPseudoAction = _ChoicesPseudoAction

        return subparsers


class _ChoicesPseudoAction(configargparse.Action):
    def __init__(self, name, aliases, help):
        super(_ChoicesPseudoAction, self).__init__(
            option_strings=[],
            dest=name,
            help=help,
            metavar=name,
        )
