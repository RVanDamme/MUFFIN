"""Lightweight Nextflow / Groovy scanner and block extractor.

The approach is pragmatic, not exhaustive: we first *mask* string literals and
comments (replacing their interior with spaces while preserving newlines and
byte offsets), then parse structural constructs — ``workflow``, ``process``,
``include``, ``function`` — on the masked view by balancing braces. The
original text is kept alongside so rules that need to inspect shell scripts
(for example, detecting embedded Python) can do so on the raw content.

This is good enough to enforce structural rules reliably across real-world
Nextflow pipelines without requiring a JVM or a Groovy parser.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# Masking: strip strings and comments to space while preserving offsets
# ---------------------------------------------------------------------------


def mask_source(src: str) -> str:
    """Return ``src`` with string/comment bodies replaced by spaces.

    Newlines are preserved so ``line = masked.count('\\n', 0, offset) + 1`` still
    gives the correct source line for any offset.
    """
    out = []
    i = 0
    n = len(src)
    while i < n:
        two = src[i : i + 2]
        three = src[i : i + 3]

        # Line comment --------------------------------------------------------
        if two == "//":
            while i < n and src[i] != "\n":
                out.append(" ")
                i += 1
            continue

        # Block comment -------------------------------------------------------
        if two == "/*":
            out.append("  ")
            i += 2
            while i < n and src[i : i + 2] != "*/":
                out.append("\n" if src[i] == "\n" else " ")
                i += 1
            if i < n:
                out.append("  ")
                i += 2
            continue

        # Triple-double-quoted string ----------------------------------------
        if three == '"""':
            out.append("   ")
            i += 3
            while i < n and src[i : i + 3] != '"""':
                out.append("\n" if src[i] == "\n" else " ")
                i += 1
            if i < n:
                out.append("   ")
                i += 3
            continue

        # Triple-single-quoted string ----------------------------------------
        if three == "'''":
            out.append("   ")
            i += 3
            while i < n and src[i : i + 3] != "'''":
                out.append("\n" if src[i] == "\n" else " ")
                i += 1
            if i < n:
                out.append("   ")
                i += 3
            continue

        # Single-line string --------------------------------------------------
        if src[i] in ('"', "'"):
            quote = src[i]
            out.append(" ")
            i += 1
            while i < n and src[i] != quote and src[i] != "\n":
                if src[i] == "\\" and i + 1 < n:
                    out.append("  ")
                    i += 2
                    continue
                out.append(" ")
                i += 1
            if i < n and src[i] == quote:
                out.append(" ")
                i += 1
            continue

        out.append(src[i])
        i += 1
    return "".join(out)


def mask_comments_only(src: str) -> str:
    """Strip ``//`` and ``/* */`` comments while preserving string bodies.

    Useful for rules that need to see Groovy string interpolations such as
    ``"${params.foo}"`` — those are real references to ``params`` even though
    they live inside a quoted string, so ``mask_source`` (which blanks
    string interiors) would miss them.
    """
    out = []
    i = 0
    n = len(src)
    while i < n:
        two = src[i : i + 2]
        if two == "//":
            while i < n and src[i] != "\n":
                out.append(" ")
                i += 1
            continue
        if two == "/*":
            out.append("  ")
            i += 2
            while i < n and src[i : i + 2] != "*/":
                out.append("\n" if src[i] == "\n" else " ")
                i += 1
            if i < n:
                out.append("  ")
                i += 2
            continue
        out.append(src[i])
        i += 1
    return "".join(out)


def line_of(src: str, offset: int) -> int:
    return src.count("\n", 0, offset) + 1


def column_of(src: str, offset: int) -> int:
    last_nl = src.rfind("\n", 0, offset)
    return offset - last_nl


# ---------------------------------------------------------------------------
# Block extraction
# ---------------------------------------------------------------------------


def find_matching_brace(masked: str, open_idx: int) -> int:
    """Given the index of a ``{``, return the index of its matching ``}``.

    Returns -1 if unbalanced.
    """
    assert masked[open_idx] == "{"
    depth = 0
    i = open_idx
    n = len(masked)
    while i < n:
        c = masked[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


@dataclass
class Block:
    kind: str  # 'workflow' | 'process' | 'function'
    name: str | None
    header_start: int  # start of ``workflow``/``process``/``function`` keyword
    body_start: int  # offset of the opening ``{``
    body_end: int  # offset of the matching ``}``


_BLOCK_KEYWORDS = re.compile(
    r"(?P<kw>\bworkflow\b|\bprocess\b)\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)?\s*\{"
)


def iter_blocks(masked: str) -> Iterator[Block]:
    """Yield all top-level ``workflow`` and ``process`` declarations.

    Block bodies may contain nested braces (closures, maps, scripts); the brace
    matcher handles that. An anonymous ``workflow { ... }`` is the implicit
    entry workflow and is yielded with ``name=None``.
    """
    pos = 0
    while True:
        m = _BLOCK_KEYWORDS.search(masked, pos)
        if not m:
            return
        open_brace = m.end() - 1
        close_brace = find_matching_brace(masked, open_brace)
        if close_brace == -1:
            return
        yield Block(
            kind=m.group("kw"),
            name=m.group("name"),
            header_start=m.start("kw"),
            body_start=open_brace,
            body_end=close_brace,
        )
        pos = close_brace + 1


# ---------------------------------------------------------------------------
# Workflow internals (take / main / emit sections, module calls)
# ---------------------------------------------------------------------------


_SECTION_RE = re.compile(
    r"^[ \t]*(?P<label>take|main|emit|script|shell|exec|output|when|input|workflow|publish)[ \t]*:[ \t]*$",
    re.MULTILINE,
)


@dataclass
class ProcessSections:
    input: tuple[int, int] | None = None
    output: tuple[int, int] | None = None
    when: tuple[int, int] | None = None
    script: tuple[int, int] | None = None
    shell: tuple[int, int] | None = None
    exec: tuple[int, int] | None = None
    publish: tuple[int, int] | None = None


@dataclass
class WorkflowSections:
    take: tuple[int, int] | None = None
    main: tuple[int, int] | None = None
    emit: tuple[int, int] | None = None
    publish: tuple[int, int] | None = None


def _split_sections(masked: str, body_start: int, body_end: int) -> dict[str, tuple[int, int]]:
    """Split the body of a process/workflow into labelled sections.

    Returns a mapping ``label -> (start, end)`` where the start is the offset
    just after the ``label:`` line and the end is the offset of the next label
    or the block's closing brace.
    """
    inner_start = body_start + 1
    inner_end = body_end
    inner = masked[inner_start:inner_end]

    hits: list[tuple[str, int, int]] = []
    for m in _SECTION_RE.finditer(inner):
        label = m.group("label")
        section_head_start = inner_start + m.start()
        section_head_end = inner_start + m.end()
        hits.append((label, section_head_start, section_head_end))

    if not hits:
        return {}

    sections: dict[str, tuple[int, int]] = {}
    for idx, (label, _head_start, head_end) in enumerate(hits):
        next_start = hits[idx + 1][1] if idx + 1 < len(hits) else inner_end
        sections[label] = (head_end, next_start)
    return sections


def workflow_sections(masked: str, block: Block) -> WorkflowSections:
    if block.kind != "workflow":
        raise ValueError("workflow_sections() requires a workflow block")
    s = _split_sections(masked, block.body_start, block.body_end)
    ws = WorkflowSections()
    if "take" in s:
        ws.take = s["take"]
    if "main" in s:
        ws.main = s["main"]
    if "emit" in s:
        ws.emit = s["emit"]
    if "publish" in s:
        ws.publish = s["publish"]
    return ws


def process_sections(masked: str, block: Block) -> ProcessSections:
    if block.kind != "process":
        raise ValueError("process_sections() requires a process block")
    s = _split_sections(masked, block.body_start, block.body_end)
    ps = ProcessSections()
    for attr in ("input", "output", "when", "script", "shell", "exec", "publish"):
        if attr in s:
            setattr(ps, attr, s[attr])
    return ps


# ---------------------------------------------------------------------------
# Include statements
# ---------------------------------------------------------------------------


@dataclass
class Include:
    names: list[str]
    source: str
    addparams: bool
    params_clause: bool
    offset: int


_INCLUDE_RE = re.compile(
    r"\binclude\s*\{\s*(?P<body>[^}]*)\s*\}\s*from\s*['\"](?P<source>[^'\"]+)['\"]"
    r"(?P<rest>[^\n;]*)"
)


def iter_includes(masked: str, raw: str) -> Iterator[Include]:
    """Find all include declarations.

    Runs against the raw source (because masking strips the quoted ``from`` path),
    but verifies via the masked source that the match's ``include`` keyword is
    real code — not part of a comment.
    """
    for m in _INCLUDE_RE.finditer(raw):
        # Reject matches that actually live inside a comment or string: in the
        # masked view, the keyword 'include' would have been blanked out.
        if masked[m.start() : m.start() + 7] != "include":
            continue
        raw_body = raw[m.start("body") : m.end("body")]
        names = [x.strip() for x in re.split(r"\s*;\s*|\s*,\s*", raw_body) if x.strip()]
        # When `include { FOO as BAR }`, the call site uses `BAR`. Record the alias
        # (last part after `as`); fall back to the original name when no alias.
        clean_names = [re.split(r"\s+as\s+", n)[-1].strip() for n in names]
        rest = m.group("rest") or ""
        yield Include(
            names=clean_names,
            source=m.group("source"),
            addparams="addParams" in rest,
            params_clause=bool(re.search(r"(?<!\.)\bparams\s*\(", rest)),
            offset=m.start(),
        )


# ---------------------------------------------------------------------------
# Call detection (module invocation sites inside workflow bodies)
# ---------------------------------------------------------------------------


_CALL_RE = re.compile(r"(?<![A-Za-z0-9_.])([A-Z][A-Z0-9_]{2,})\s*\(")


@dataclass
class Call:
    name: str
    offset: int


def iter_calls(masked: str, start: int, end: int) -> Iterator[Call]:
    for m in _CALL_RE.finditer(masked, start, end):
        yield Call(name=m.group(1), offset=m.start(1))


# ---------------------------------------------------------------------------
# Top-level parse result
# ---------------------------------------------------------------------------


@dataclass
class ParsedFile:
    path: Path
    raw: str
    masked: str
    blocks: list[Block] = field(default_factory=list)
    includes: list[Include] = field(default_factory=list)

    @property
    def is_main_nf(self) -> bool:
        return (
            self.path.name == "main.nf"
            and "modules" not in self.path.parts
            and "subworkflows" not in self.path.parts
        )

    @property
    def is_subworkflow(self) -> bool:
        return "subworkflows" in self.path.parts

    @property
    def is_module(self) -> bool:
        return "modules" in self.path.parts

    @property
    def is_workflow_file(self) -> bool:
        return "workflows" in self.path.parts

    def line_of(self, offset: int) -> int:
        return line_of(self.raw, offset)

    def col_of(self, offset: int) -> int:
        return column_of(self.raw, offset)


def parse_file(path: Path) -> ParsedFile:
    raw = path.read_text(encoding="utf-8")
    masked = mask_source(raw)
    parsed = ParsedFile(path=path, raw=raw, masked=masked)
    parsed.blocks = list(iter_blocks(masked))
    parsed.includes = list(iter_includes(masked, raw))
    return parsed
