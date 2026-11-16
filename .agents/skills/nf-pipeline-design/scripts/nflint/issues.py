"""Issue records emitted by rules."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


@dataclass(frozen=True)
class Issue:
    path: Path
    line: int
    column: int
    code: str
    severity: Severity
    message: str

    def format_text(self, root: Path | None = None) -> str:
        display = self.path
        if root is not None:
            try:
                display = self.path.relative_to(root)
            except ValueError:
                display = self.path
        return (
            f"{display}:{self.line}:{self.column}: "
            f"{self.severity.value.upper():<7} {self.code} {self.message}"
        )

    def to_dict(self, root: Path | None = None) -> dict:
        display = self.path
        if root is not None:
            try:
                display = self.path.relative_to(root)
            except ValueError:
                display = self.path
        return {
            "path": str(display),
            "line": self.line,
            "column": self.column,
            "code": self.code,
            "severity": self.severity.value,
            "message": self.message,
        }
