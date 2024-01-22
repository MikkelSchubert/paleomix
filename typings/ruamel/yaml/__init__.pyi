from pathlib import Path
from typing import IO

from typing_extensions import Literal

class YAML:
    version: tuple[int, int] | None
    def __init__(
        self,
        *,
        typ: Literal["rt", "safe", "unsafe", "base"] | None = None,
        pure: bool = False,
        output: IO[str] | None = None,
        plug_ins: list[str] | None = None,
    ) -> None: ...
    def load(self, stream: Path | IO[str]) -> object: ...
