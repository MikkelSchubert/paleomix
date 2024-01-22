from typing import IO, AnyStr

from typing_extensions import TypeAlias

Color: TypeAlias = str | int | tuple[int, int, int] | None

def ansi_wrap(
    text: str,
    color: Color = None,
    background: Color = None,
    bright: bool = False,
    bold: bool = False,
) -> str: ...
def terminal_supports_colors(stream: IO[AnyStr] | None = None) -> bool: ...
