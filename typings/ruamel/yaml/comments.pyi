class LineCol:
    line: int
    col: int
    def item(self, idx: object) -> tuple[int, int] | None: ...

class CommentedBase:
    @property
    def lc(self) -> LineCol: ...
