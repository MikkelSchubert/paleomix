import argparse
from collections import OrderedDict
from collections.abc import Sequence
from typing import IO, Any, Callable, TypeVar, overload

from _typeshed import Incomplete
from typing_extensions import Protocol, TypeAlias

HelpFormatter = argparse.HelpFormatter
RawDescriptionHelpFormatter = argparse.RawDescriptionHelpFormatter
RawTextHelpFormatter = argparse.RawTextHelpFormatter
ArgumentDefaultsHelpFormatter = argparse.ArgumentDefaultsHelpFormatter
ArgumentError = argparse.ArgumentError
ArgumentTypeError = argparse.ArgumentTypeError
Action = argparse.Action
FileType = argparse.FileType
Namespace = argparse.Namespace

class ArgumentDefaultsRawHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawTextHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
): ...

_ConfigDict: TypeAlias = OrderedDict[str, str | list[str]]

class ConfigFileParser:
    def get_syntax_description(self) -> str: ...
    def parse(self, stream: IO[str]) -> _ConfigDict: ...
    def serialize(self, items: _ConfigDict) -> str: ...

class ConfigFileParserException(Exception): ...  # noqa: N818

class DefaultConfigFileParser(ConfigFileParser):
    def get_syntax_description(self) -> str: ...
    def parse(self, stream: IO[str]) -> _ConfigDict: ...
    def serialize(self, items: _ConfigDict) -> str: ...

class ConfigparserConfigFileParser(ConfigFileParser):
    def get_syntax_description(self) -> str: ...
    def parse(self, stream: IO[str]) -> _ConfigDict: ...
    def serialize(self, items: _ConfigDict) -> str: ...

class YAMLConfigFileParser(ConfigFileParser):
    def get_syntax_description(self) -> str: ...
    def parse(self, stream: IO[str]) -> _ConfigDict: ...
    def serialize(
        self,
        items: _ConfigDict,
        default_flow_style: bool = ...,
    ) -> str: ...

class _FormatterClass(Protocol):
    def __call__(self, *, prog: str) -> HelpFormatter: ...

_ArgumentParserType = TypeVar("_ArgumentParserType", bound=ArgumentParser)

class ArgumentParser(argparse.ArgumentParser):
    def __init__(
        self,
        prog: str | None = None,
        usage: str | None = None,
        description: str | None = None,
        epilog: str | None = None,
        parents: Sequence[ArgumentParser] = [],
        formatter_class: _FormatterClass = ...,
        prefix_chars: str = "-",
        fromfile_prefix_chars: str | None = None,
        argument_default: object = None,
        conflict_handler: str = "error",
        *,
        add_help: bool = True,
        allow_abbrev: bool = True,
        exit_on_error: bool = True,
        # Configargparse arguments
        add_config_file_help: bool = True,
        add_env_var_help: bool = True,
        auto_env_var_prefix: str | None = None,
        default_config_files: Sequence[str] = ...,
        ignore_unknown_config_file_keys: bool = False,
        config_file_open_func: Callable[[str, str], IO[str]] = ...,
        config_file_parser_class: ConfigFileParser = ...,
        args_for_setting_config_path: Sequence[str] = ...,
        config_arg_is_required: bool = False,
        config_arg_help_message: str = ...,
        args_for_writing_out_config_file: Sequence[str] = ...,
        write_out_config_file_arg_help_message: str = ...,
    ) -> None: ...
    def parse_args(
        self,
        args: Sequence[str] | None = ...,
        namespace: Namespace | None = ...,
        config_file_contents: str | None = ...,
        env_vars: dict[Any, Any] = ...,
    ) -> Namespace: ...
    def parse_known_args(
        self,
        *,
        args: Sequence[str] | None = ...,
        namespace: Namespace | None = ...,
        config_file_contents: Incomplete | None = ...,
        env_vars: dict[Any, Any] = ...,
        ignore_help_args: bool = ...,
    ) -> tuple[Namespace, list[str]]: ...
    def write_config_file(
        self,
        parsed_namespace: Namespace,
        output_file_paths: Sequence[str],
        exit_after: bool = ...,
    ) -> None: ...
    def get_possible_config_keys(self, action: Action) -> list[str]: ...
    def convert_item_to_command_line_arg(
        self, action: Action, key: str, value: str | list[str]
    ) -> list[str]: ...
    @overload
    def add_subparsers(
        self: _ArgumentParserType,
        *,
        title: str = ...,
        description: str | None = ...,
        prog: str = ...,
        action: type[Action] = ...,
        option_string: str = ...,
        dest: str | None = ...,
        required: bool = ...,
        help: str | None = ...,
        metavar: str | None = ...,
    ) -> argparse._SubParsersAction[_ArgumentParserType]: ...  # pyright: ignore[reportPrivateUsage]
    @overload
    def add_subparsers(
        self,
        *,
        title: str = ...,
        description: str | None = ...,
        prog: str = ...,
        parser_class: type[_ArgumentParserType],
        action: type[Action] = ...,
        option_string: str = ...,
        dest: str | None = ...,
        required: bool = ...,
        help: str | None = ...,
        metavar: str | None = ...,
    ) -> argparse._SubParsersAction[_ArgumentParserType]: ...  # pyright: ignore[reportPrivateUsage]

ONE_OR_MORE = argparse.ONE_OR_MORE
OPTIONAL = argparse.OPTIONAL
REMAINDER = argparse.REMAINDER
SUPPRESS = argparse.SUPPRESS
ZERO_OR_MORE = argparse.ZERO_OR_MORE

RawFormatter = RawDescriptionHelpFormatter
DefaultsFormatter = ArgumentDefaultsHelpFormatter
DefaultsRawFormatter = ArgumentDefaultsRawHelpFormatter
DefaultsRawFormatter = ArgumentDefaultsRawHelpFormatter

RawFormatter = RawDescriptionHelpFormatter
DefaultsFormatter = ArgumentDefaultsHelpFormatter
DefaultsRawFormatter = ArgumentDefaultsRawHelpFormatter
DefaultsRawFormatter = ArgumentDefaultsRawHelpFormatter
