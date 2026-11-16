"""Rule implementations.

Each rule is a callable ``(ParsedFile) -> Iterable[Issue]``. Rules are grouped
by applicable file kind (main.nf / subworkflow / module / any) and registered
with an integer code so users can ``--select`` or ``--ignore`` them.

Rule codes:

  NF001 main.nf calls a module (UPPER_SNAKE) directly instead of a subworkflow
  NF002 main.nf contains a deep operator chain (>3 chained ops on a channel)
  NF003 main.nf references params.* outside the workflow setup area
  NF004 main.nf is missing a setup / run section comment structure
  NF005 include statement uses deprecated addParams / params clause
  NF010 subworkflow missing a take: block
  NF011 subworkflow missing an emit: block
  NF012 subworkflow has unnamed emit outputs
  NF013 subworkflow executes >1 module per execution path (umbrella anti-pattern)
  NF014 multi-module subworkflow does not emit a 'versions' channel
  NF020 module references params.* (incl. inside a process script body)
  NF021 process uses a when: block for branching
  NF022 module script contains inline glue in a non-Nextflow language
  NF023 process name is not UPPER_SNAKE_CASE
  NF024 subworkflow references params.* directly (pass as explicit input)
  NF030 deprecated Channel.* call (prefer lowercase channel.*)
  NF031 implicit environment variable (${PWD}) instead of env('PWD')
  NF032 import statement used (disallowed in strict syntax)
  NF033 for / while loop used (prefer higher-order operators)
  NF034 class declaration used (move to lib/)
  NF040 subworkflow/process output consumed by numeric index (.out[N])
  NF041 process script silently swallows failures (set +e, || true, || echo)
  NF042 interpreter invoked with absolute path (call scripts by name on PATH)
  NF043 process script copies a staged input file instead of using it
  NF050 module process missing container directive
  NF051 pipeline missing nextflow_schema.json next to main.nf
  NF052 module process missing resource declaration (label or cpus/memory/time)
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable
from dataclasses import dataclass

from .issues import Issue, Severity
from .parser import (
    Block,
    ParsedFile,
    iter_calls,
    mask_comments_only,
    process_sections,
    workflow_sections,
)


def _in_shell_comment(raw: str, offset: int) -> bool:
    """Return True if ``offset`` falls on a line whose first non-blank char is '#'.

    Used when scanning the raw body of a process ``script:`` section to
    avoid flagging example patterns that live inside a bash-style comment.
    """
    line_start = raw.rfind("\n", 0, offset) + 1
    i = line_start
    n = len(raw)
    while i < n and raw[i] in " \t":
        i += 1
    return i < n and raw[i] == "#"


# ---------------------------------------------------------------------------
# Rule registry
# ---------------------------------------------------------------------------


RuleFn = Callable[[ParsedFile], Iterable[Issue]]


@dataclass(frozen=True)
class Rule:
    code: str
    severity: Severity
    title: str
    fn: RuleFn


_REGISTRY: list[Rule] = []


def rule(code: str, severity: Severity, title: str) -> Callable[[RuleFn], RuleFn]:
    def deco(fn: RuleFn) -> RuleFn:
        _REGISTRY.append(Rule(code=code, severity=severity, title=title, fn=fn))
        return fn

    return deco


def all_rules() -> list[Rule]:
    return list(_REGISTRY)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_issue(parsed: ParsedFile, offset: int, code: str, severity: Severity, msg: str) -> Issue:
    return Issue(
        path=parsed.path,
        line=parsed.line_of(offset),
        column=parsed.col_of(offset),
        code=code,
        severity=severity,
        message=msg,
    )


def _entry_workflow(parsed: ParsedFile) -> Block | None:
    for b in parsed.blocks:
        if b.kind == "workflow" and b.name is None:
            return b
    return None


def _named_workflows(parsed: ParsedFile) -> list[Block]:
    return [b for b in parsed.blocks if b.kind == "workflow" and b.name is not None]


def _processes(parsed: ParsedFile) -> list[Block]:
    return [b for b in parsed.blocks if b.kind == "process"]


# ---------------------------------------------------------------------------
# main.nf rules
# ---------------------------------------------------------------------------


@rule("NF001", Severity.ERROR, "main.nf should call subworkflows, not modules directly")
def nf001_main_direct_module(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_main_nf:
        return
    # Build the set of symbols imported from ``modules/...``. Anything matching
    # an UPPER_SNAKE call inside the entry workflow that was imported from
    # ``modules/`` is a module invocation (which SKILL.md forbids in main.nf).
    module_names: set[str] = set()
    for inc in parsed.includes:
        if "modules" in inc.source.replace("\\", "/").split("/"):
            module_names.update(inc.names)

    entry = _entry_workflow(parsed)
    if entry is None:
        return
    sections = workflow_sections(parsed.masked, entry)
    main_range = sections.main or (entry.body_start + 1, entry.body_end)

    for call in iter_calls(parsed.masked, main_range[0], main_range[1]):
        if call.name in module_names:
            yield _make_issue(
                parsed,
                call.offset,
                "NF001",
                Severity.ERROR,
                f"main.nf calls module '{call.name}' directly; invoke it from a subworkflow instead",
            )


@rule("NF002", Severity.WARNING, "main.nf should avoid deep operator chains")
def nf002_main_operator_chain(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_main_nf:
        return
    entry = _entry_workflow(parsed)
    if entry is None:
        return
    sections = workflow_sections(parsed.masked, entry)
    main_range = sections.main or (entry.body_start + 1, entry.body_end)
    body = parsed.masked[main_range[0] : main_range[1]]
    # Heuristic: ≥4 ``.method`` call sites on one physical line form a "deep"
    # chain. This catches both ``.foo()`` and bare-closure ``.foo { ... }``
    # styles without trying to balance nested braces.
    line_offset = main_range[0]
    for raw_line in body.split("\n"):
        # Count method-call dots: ``.name`` where name is a Groovy identifier.
        dots = re.findall(r"\.\s*[a-zA-Z_][a-zA-Z0-9_]*\s*[\({]", raw_line)
        if len(dots) >= 4:
            yield _make_issue(
                parsed,
                line_offset,
                "NF002",
                Severity.WARNING,
                f"deep operator chain ({len(dots)} calls) in main.nf; "
                "move channel-shaping logic into a subworkflow",
            )
        line_offset += len(raw_line) + 1  # +1 for the stripped newline


@rule("NF003", Severity.ERROR, "main.nf should not access params.* outside the setup area")
def nf003_main_params_outside_setup(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_main_nf:
        return
    entry = _entry_workflow(parsed)
    if entry is None:
        return
    sections = workflow_sections(parsed.masked, entry)
    main_range = sections.main or (entry.body_start + 1, entry.body_end)

    # Split main: the SETUP is everything up to the first named-call site;
    # after that, references to params.* are considered "deep".
    body = parsed.masked[main_range[0] : main_range[1]]
    first_call = None
    for call in iter_calls(parsed.masked, main_range[0], main_range[1]):
        first_call = call.offset
        break
    if first_call is None:
        return

    for m in re.finditer(r"\bparams\.\w+", body):
        abs_off = main_range[0] + m.start()
        if abs_off > first_call:
            yield _make_issue(
                parsed,
                abs_off,
                "NF003",
                Severity.ERROR,
                "params.* used after the setup section; pass explicit inputs to subworkflows instead",
            )


@rule("NF004", Severity.WARNING, "main.nf should have setup / run section comments")
def nf004_main_section_comments(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_main_nf:
        return
    entry = _entry_workflow(parsed)
    if entry is None:
        return
    # Look at the RAW body of the entry workflow (we want to see the comments).
    body = parsed.raw[entry.body_start : entry.body_end]
    has_setup = bool(
        re.search(r"(?i)(parameter\s*setup|setup|parameters?)\s*$", body, re.MULTILINE)
    )
    has_run = bool(
        re.search(r"(?i)(pipeline\s*run|run\s*section|analysis\s*stages?)\s*$", body, re.MULTILINE)
    )
    # Only warn if neither marker is found and the workflow body is non-trivial:
    # at least 5 non-blank, non-brace body lines suggest the author has enough
    # steps to benefit from a 'Parameter setup' / 'Pipeline run' structure.
    meaningful_lines = sum(
        1 for line in body.splitlines() if line.strip() and line.strip() not in ("{", "}")
    )
    if meaningful_lines >= 5 and not (has_setup and has_run):
        yield _make_issue(
            parsed,
            entry.body_start,
            "NF004",
            Severity.WARNING,
            "main.nf workflow should be organised into a 'Parameter setup' and 'Pipeline run' section",
        )


@rule("NF005", Severity.ERROR, "include statements should not use deprecated addParams / params")
def nf005_include_addparams(parsed: ParsedFile) -> Iterable[Issue]:
    for inc in parsed.includes:
        if inc.addparams:
            yield _make_issue(
                parsed,
                inc.offset,
                "NF005",
                Severity.ERROR,
                "include uses deprecated .addParams(...); pass params as explicit inputs",
            )
        if inc.params_clause:
            yield _make_issue(
                parsed,
                inc.offset,
                "NF005",
                Severity.ERROR,
                "include uses deprecated .params(...) clause; pass params as explicit inputs",
            )


# ---------------------------------------------------------------------------
# Subworkflow rules
# ---------------------------------------------------------------------------


@rule("NF010", Severity.ERROR, "subworkflow must declare a take: block")
def nf010_subworkflow_take(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_subworkflow:
        return
    for wf in _named_workflows(parsed):
        s = workflow_sections(parsed.masked, wf)
        if s.take is None:
            yield _make_issue(
                parsed,
                wf.header_start,
                "NF010",
                Severity.ERROR,
                f"subworkflow '{wf.name}' is missing a take: block",
            )


@rule("NF011", Severity.ERROR, "subworkflow must declare an emit: block")
def nf011_subworkflow_emit(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_subworkflow:
        return
    for wf in _named_workflows(parsed):
        s = workflow_sections(parsed.masked, wf)
        if s.emit is None:
            yield _make_issue(
                parsed,
                wf.header_start,
                "NF011",
                Severity.ERROR,
                f"subworkflow '{wf.name}' is missing an emit: block",
            )


@rule("NF012", Severity.WARNING, "subworkflow emit outputs should be named")
def nf012_emit_named(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_subworkflow:
        return
    for wf in _named_workflows(parsed):
        s = workflow_sections(parsed.masked, wf)
        if s.emit is None:
            continue
        emit_body = parsed.masked[s.emit[0] : s.emit[1]]
        for raw_line in emit_body.split("\n"):
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("//"):
                continue
            # A named emit is ``name = expression`` or just ``name`` (bare
            # channel name). Things like ``PROC.out.chan`` or a closing ``}``
            # are not named.
            if stripped in ("}", "{"):
                continue
            if "=" in stripped:
                continue
            if re.fullmatch(r"[a-z_][a-zA-Z0-9_]*", stripped):
                continue
            # Unnamed emit (e.g. ``PROC.out`` or a channel expression).
            offset = s.emit[0] + raw_line_index(emit_body, raw_line)
            yield _make_issue(
                parsed,
                offset,
                "NF012",
                Severity.WARNING,
                "emit output is not named; use 'name = expr' so callers consume by name",
            )


def raw_line_index(body: str, raw_line: str) -> int:
    """Return the offset within ``body`` of the first occurrence of ``raw_line``."""
    idx = body.find(raw_line)
    return max(idx, 0)


def _skip_balanced(masked: str, open_i: int, open_ch: str, close_ch: str, end: int) -> int:
    """Return the offset *past* the ``close_ch`` that balances ``masked[open_i]``.

    Assumes ``masked[open_i] == open_ch``. Works on the masked view, so strings
    and comments have been replaced with whitespace — brace/paren counts reflect
    actual code structure, not text inside literals.
    """
    depth = 0
    i = open_i
    while i < end:
        ch = masked[i]
        if ch == open_ch:
            depth += 1
        elif ch == close_ch:
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return end


_IF_KW_RE = re.compile(r"\bif\s*\(")
_ELSE_KW_RE = re.compile(r"\belse\b")
_WORD_RE = re.compile(r"\w")


def _max_module_calls_per_path(masked: str, start: int, end: int, module_names: set[str]) -> int:
    """Count module calls along a single execution path in ``masked[start:end]``.

    Sequential top-level calls are additive. Calls inside mutually exclusive
    ``if`` / ``else if`` / ``else`` branches are treated as alternatives — only
    one branch runs per execution — so the branch contribution is the MAX across
    branches, not the sum.

    This matches the skill rule: a subworkflow should execute one module per
    execution path; multiple interchangeable modules are fine only when selected
    by workflow-level ``if/else``.
    """
    i = start
    count = 0
    while i < end:
        m = _IF_KW_RE.search(masked, i, end)

        # Count sequential module calls up to the next `if` (or to `end`).
        scan_end = m.start() if m else end
        for call in iter_calls(masked, i, scan_end):
            if call.name in module_names:
                count += 1

        if not m:
            break

        # Found an `if (…) { … }` at m.start(). Parse the condition + body,
        # then any chained `else if` / `else` branches.
        cond_open = masked.find("(", m.start(), end)
        if cond_open == -1:
            i = m.end()
            continue
        cond_close = _skip_balanced(masked, cond_open, "(", ")", end)

        # Locate the body: must be a `{ … }` block (we don't handle single-stmt
        # bodies without braces — not idiomatic in Nextflow workflow bodies).
        j = cond_close
        while j < end and masked[j].isspace():
            j += 1
        if j >= end or masked[j] != "{":
            i = cond_close
            continue
        body_open = j
        body_close = _skip_balanced(masked, body_open, "{", "}", end)

        branch_counts = [
            _max_module_calls_per_path(masked, body_open + 1, body_close - 1, module_names)
        ]

        # Chase any `else` / `else if` branches.
        k = body_close
        while True:
            while k < end and masked[k].isspace():
                k += 1
            em = _ELSE_KW_RE.match(masked, k)
            if not em or (em.end() < end and _WORD_RE.match(masked, em.end())):
                break
            k = em.end()
            while k < end and masked[k].isspace():
                k += 1
            # `else if (...)` — skip the condition before the body brace.
            if masked.startswith("if", k) and (k + 2 == end or not _WORD_RE.match(masked, k + 2)):
                k += 2
                while k < end and masked[k].isspace():
                    k += 1
                if k < end and masked[k] == "(":
                    k = _skip_balanced(masked, k, "(", ")", end)
                while k < end and masked[k].isspace():
                    k += 1
            if k >= end or masked[k] != "{":
                break
            b_open = k
            b_close = _skip_balanced(masked, b_open, "{", "}", end)
            branch_counts.append(
                _max_module_calls_per_path(masked, b_open + 1, b_close - 1, module_names)
            )
            k = b_close

        count += max(branch_counts) if branch_counts else 0
        i = k
    return count


@rule(
    "NF013",
    Severity.WARNING,
    "subworkflow executes >1 module per execution path (umbrella anti-pattern)",
)
def nf013_subworkflow_size(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_subworkflow:
        return
    module_names: set[str] = set()
    for inc in parsed.includes:
        src = inc.source.replace("\\", "/")
        if "modules" in src.split("/"):
            module_names.update(inc.names)
    for wf in _named_workflows(parsed):
        s = workflow_sections(parsed.masked, wf)
        main_range = s.main or (wf.body_start + 1, wf.body_end)
        max_per_path = _max_module_calls_per_path(
            parsed.masked, main_range[0], main_range[1], module_names
        )
        if max_per_path > 1:
            yield _make_issue(
                parsed,
                wf.header_start,
                "NF013",
                Severity.WARNING,
                f"subworkflow '{wf.name}' runs {max_per_path} modules per execution path; "
                "a subworkflow should execute one module per path — put alternatives "
                "inside workflow-level if/else, and split sequential steps into separate "
                "subworkflows",
            )


@rule("NF014", Severity.WARNING, "multi-module subworkflow should emit a versions channel")
def nf014_versions_channel(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_subworkflow:
        return
    module_names: set[str] = set()
    for inc in parsed.includes:
        src = inc.source.replace("\\", "/")
        if "modules" in src.split("/"):
            module_names.update(inc.names)
    for wf in _named_workflows(parsed):
        s = workflow_sections(parsed.masked, wf)
        if s.emit is None:
            continue
        main_range = s.main or (wf.body_start + 1, wf.body_end)
        calls = [c.name for c in iter_calls(parsed.masked, main_range[0], main_range[1])]
        module_calls = set(c for c in calls if c in module_names)
        if len(module_calls) < 2:
            continue
        emit_body = parsed.masked[s.emit[0] : s.emit[1]]
        if not re.search(r"\bversions\b", emit_body):
            yield _make_issue(
                parsed,
                wf.header_start,
                "NF014",
                Severity.WARNING,
                f"subworkflow '{wf.name}' runs ≥2 modules but does not emit a 'versions' channel",
            )


# ---------------------------------------------------------------------------
# Module rules
# ---------------------------------------------------------------------------


@rule("NF020", Severity.ERROR, "module must not reference params.*")
def nf020_module_params(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_module:
        return
    # Groovy-level references (directives, output expressions, etc.) — these
    # are visible in the masked view because they're outside script strings.
    for m in re.finditer(r"\bparams\.\w+", parsed.masked):
        yield _make_issue(
            parsed,
            m.start(),
            "NF020",
            Severity.ERROR,
            "module references params.* — pass values as explicit process inputs",
        )
    # References hidden inside a process script body. The masker blanks
    # triple-quoted script bodies, so a string-interpolated ``${params.foo}``
    # inside a shell script is invisible to the outer scan. Look at the raw
    # script text directly.
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if ps.script is None:
            continue
        raw = parsed.raw[ps.script[0] : ps.script[1]]
        for m in re.finditer(r"\$\{?\s*params\.\w+", raw):
            if _in_shell_comment(raw, m.start()):
                continue
            off = ps.script[0] + m.start()
            yield _make_issue(
                parsed,
                off,
                "NF020",
                Severity.ERROR,
                "module script interpolates params.* — declare it as a process input",
            )


@rule("NF021", Severity.WARNING, "process should not use a when: block")
def nf021_process_when(parsed: ParsedFile) -> Iterable[Issue]:
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if ps.when is not None:
            yield _make_issue(
                parsed,
                ps.when[0],
                "NF021",
                Severity.WARNING,
                f"process '{proc.name}' uses when: — move branching to the calling workflow",
            )


# Languages that SKILL.md explicitly names as inappropriate for channel-glue
# logic embedded in a module script. Python was the original offender, but the
# principle extends to any general-purpose language.
_GLUE_LANGS = r"(?:python\d?|perl|ruby|node|Rscript|julia|lua)"


@rule("NF022", Severity.WARNING, "module script should not embed inline glue in another language")
def nf022_inline_glue(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_module:
        return
    patterns = [
        # `perl -e '...'`, `ruby -e '...'`, `python -c "..."`, `Rscript -e ...`
        rf"\b{_GLUE_LANGS}\s+-[ec]\b",
        # `python3 - <<'PY'`, `Rscript - <<'R'` — stdin heredoc forms
        rf"\b{_GLUE_LANGS}\s+-\s*<<",
    ]
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if ps.script is None:
            continue
        raw_script = parsed.raw[ps.script[0] : ps.script[1]]
        hit = False
        for pat in patterns:
            for m in re.finditer(pat, raw_script):
                if _in_shell_comment(raw_script, m.start()):
                    continue
                yield _make_issue(
                    parsed,
                    ps.script[0] + m.start(),
                    "NF022",
                    Severity.WARNING,
                    f"process '{proc.name}' embeds inline {m.group(0).split()[0]} glue; "
                    "express channel logic in Nextflow operators",
                )
                hit = True
                break
            if hit:
                break


@rule("NF023", Severity.ERROR, "process name must be UPPER_SNAKE_CASE")
def nf023_process_name(parsed: ParsedFile) -> Iterable[Issue]:
    for proc in _processes(parsed):
        if proc.name is None:
            continue
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*", proc.name):
            yield _make_issue(
                parsed,
                proc.header_start,
                "NF023",
                Severity.ERROR,
                f"process name '{proc.name}' should be UPPER_SNAKE_CASE",
            )


@rule("NF024", Severity.WARNING, "subworkflow should not reference params.* directly")
def nf024_subworkflow_params(parsed: ParsedFile) -> Iterable[Issue]:
    if not parsed.is_subworkflow:
        return
    # Scan comment-stripped raw source — that way we catch both plain
    # ``params.foo`` Groovy references and ``"${params.foo}"`` string
    # interpolations, without false-positives from comments.
    scannable = mask_comments_only(parsed.raw)
    for m in re.finditer(r"\bparams\.\w+", scannable):
        yield _make_issue(
            parsed,
            m.start(),
            "NF024",
            Severity.WARNING,
            "subworkflow reads params.* directly — pass the value through take: instead",
        )


# ---------------------------------------------------------------------------
# General strict-syntax hygiene
# ---------------------------------------------------------------------------


@rule("NF030", Severity.WARNING, "prefer lowercase channel.* over Channel.*")
def nf030_channel_case(parsed: ParsedFile) -> Iterable[Issue]:
    for m in re.finditer(r"(?<![A-Za-z0-9_.])Channel\s*\.", parsed.masked):
        yield _make_issue(
            parsed,
            m.start(),
            "NF030",
            Severity.WARNING,
            "use lowercase 'channel.*' factory; 'Channel.*' is deprecated in strict syntax",
        )


@rule("NF031", Severity.WARNING, "avoid implicit environment variables in scripts")
def nf031_implicit_env(parsed: ParsedFile) -> Iterable[Issue]:
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if ps.script is None:
            continue
        # Look at raw script text (interpolations are preserved even if we
        # masked the enclosing triple-quotes).
        raw_script = parsed.raw[ps.script[0] : ps.script[1]]
        for m in re.finditer(r"\$\{(PWD|HOME|USER|HOSTNAME|PATH)\}", raw_script):
            off = ps.script[0] + m.start()
            yield _make_issue(
                parsed,
                off,
                "NF031",
                Severity.WARNING,
                f"implicit env variable '${{{m.group(1)}}}' — prefer env('{m.group(1)}')",
            )


@rule("NF032", Severity.ERROR, "import statements are disallowed in strict syntax")
def nf032_import(parsed: ParsedFile) -> Iterable[Issue]:
    for m in re.finditer(r"^\s*import\s+[a-zA-Z_][\w.]*", parsed.masked, re.MULTILINE):
        yield _make_issue(
            parsed,
            m.start(),
            "NF032",
            Severity.ERROR,
            "import statement is disallowed in strict syntax; use fully-qualified names",
        )


@rule("NF033", Severity.ERROR, "for / while loops are disallowed in strict syntax")
def nf033_loops(parsed: ParsedFile) -> Iterable[Issue]:
    for m in re.finditer(r"\b(for|while)\s*\(", parsed.masked):
        yield _make_issue(
            parsed,
            m.start(),
            "NF033",
            Severity.ERROR,
            f"'{m.group(1)}' loop is disallowed in strict syntax; use .each / .collect / .find",
        )


@rule("NF034", Severity.ERROR, "class declarations are disallowed in strict syntax")
def nf034_class(parsed: ParsedFile) -> Iterable[Issue]:
    for m in re.finditer(
        r"(?m)^\s*(?:public\s+|private\s+|abstract\s+)*class\s+[A-Z]\w*", parsed.masked
    ):
        yield _make_issue(
            parsed,
            m.start(),
            "NF034",
            Severity.ERROR,
            "class declarations are disallowed in strict syntax; move custom types to lib/",
        )


# ---------------------------------------------------------------------------
# Additional structural rules (closing adversarial-review blind spots)
# ---------------------------------------------------------------------------


@rule("NF040", Severity.WARNING, "consume outputs by name, not by numeric index")
def nf040_numeric_out_access(parsed: ParsedFile) -> Iterable[Issue]:
    """Flag ``SUBWF.out[0]`` / ``PROC.out[1]`` access.

    SKILL.md: "Emit named outputs and consume them by name, not by numeric
    position." A named output is safer and survives refactors; an index
    breaks silently if outputs are reordered.
    """
    for m in re.finditer(r"([A-Z][A-Z0-9_]{2,})\s*\.\s*out\s*\[\s*\d+\s*\]", parsed.masked):
        yield _make_issue(
            parsed,
            m.start(),
            "NF040",
            Severity.WARNING,
            f"'{m.group(1)}.out[N]' — consume outputs by name via an emit: label instead",
        )


@rule("NF041", Severity.WARNING, "process script must not silently swallow failures")
def nf041_swallowed_failures(parsed: ParsedFile) -> Iterable[Issue]:
    """Flag common ways a script hides a failing command from Nextflow.

    SKILL.md: "Let it fail loudly... do not hide failures by returning
    placeholder files, swallowing exit codes, emitting generic empty
    channels."
    """
    patterns = [
        (r"\bset\s+\+e\b", "'set +e' disables errexit; the task cannot fail loudly"),
        (r"\|\|\s*true\b", "'|| true' swallows the exit code of the previous command"),
        (r"\|\|\s*:(?!\w)", "'|| :' swallows the exit code of the previous command"),
        (
            r"\|\|\s*echo\b[^|\n]*>\s*\S",
            "'|| echo ... > file' silently substitutes a placeholder for the real output",
        ),
    ]
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if ps.script is None:
            continue
        raw = parsed.raw[ps.script[0] : ps.script[1]]
        for pat, msg in patterns:
            for m in re.finditer(pat, raw):
                # Skip matches that live on a shell comment line — authors
                # often mention these patterns in example prose.
                if _in_shell_comment(raw, m.start()):
                    continue
                yield _make_issue(
                    parsed,
                    ps.script[0] + m.start(),
                    "NF041",
                    Severity.WARNING,
                    f"process '{proc.name}': {msg}",
                )


@rule("NF042", Severity.WARNING, "call scripts by name on PATH, not by absolute path")
def nf042_absolute_script_path(parsed: ParsedFile) -> Iterable[Issue]:
    """Flag ``python3 /opt/.../foo.py`` and similar.

    SKILL.md: "give it a shebang, make it executable, call it directly by
    name." Absolute paths bake deployment layout into the module and
    break portability across executors and Wave containers.
    """
    # Interpreters that SKILL.md's example list applies to.
    langs = r"(?:python\d?|perl|ruby|Rscript|node|bash|sh)"
    pattern = rf"\b{langs}\s+(/\S+)"
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if ps.script is None:
            continue
        raw = parsed.raw[ps.script[0] : ps.script[1]]
        for m in re.finditer(pattern, raw):
            if _in_shell_comment(raw, m.start()):
                continue
            arg = m.group(1)
            # Skip path-like shapes that aren't really "absolute script"
            # invocations (stdin, proc, tmp scratch).
            if arg.startswith(("/dev/", "/proc/", "/tmp/")):
                continue
            yield _make_issue(
                parsed,
                ps.script[0] + m.start(),
                "NF042",
                Severity.WARNING,
                f"process '{proc.name}' invokes an interpreter with absolute path '{arg}'; "
                "place the script in bin/ and call it by name instead",
            )


def _path_input_names(input_body: str) -> set[str]:
    """Extract the variable names of ``path``-typed process inputs.

    Handles the three shapes SKILL.md shows:
      - ``path fasta``
      - ``path(fasta)``
      - ``tuple val(meta), path(fasta)`` (name is captured from the ``path(...)`` form)
    """
    names: set[str] = set()
    for m in re.finditer(
        r"\bpath\s*(?:\(\s*(?P<paren>\w+)(?:\s*,[^)]*)?\s*\)|(?P<bare>\w+))",
        input_body,
    ):
        name = m.group("paren") or m.group("bare")
        if name:
            names.add(name)
    return names


@rule("NF043", Severity.WARNING, "do not copy a staged input file inside the task script")
def nf043_copy_staged_input(parsed: ParsedFile) -> Iterable[Issue]:
    """Flag ``cp ${fasta} local.fa`` / ``cat ${fasta} > local.fa``.

    SKILL.md has an explicit Wrong/Correct example: if a file is declared
    in ``input:`` as a ``path``, Nextflow already stages it into the task
    work directory and the script should consume the staged path directly.
    Copying it defeats caching and wastes disk.
    """
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if ps.input is None or ps.script is None:
            continue
        input_raw = parsed.raw[ps.input[0] : ps.input[1]]
        path_vars = _path_input_names(input_raw)
        if not path_vars:
            continue
        raw = parsed.raw[ps.script[0] : ps.script[1]]
        for var in path_vars:
            pat = rf"\b(?:cp|cat)\s+\$\{{\s*{re.escape(var)}\s*\}}"
            for m in re.finditer(pat, raw):
                if _in_shell_comment(raw, m.start()):
                    continue
                yield _make_issue(
                    parsed,
                    ps.script[0] + m.start(),
                    "NF043",
                    Severity.WARNING,
                    f"process '{proc.name}' copies staged input '{var}' inside the task script; "
                    "use the staged path directly",
                )


# ---------------------------------------------------------------------------
# Environment and operability rules (containers, schema, resources)
# ---------------------------------------------------------------------------


def _directive_exists(parsed: ParsedFile, proc: Block, ps: ProcessSections, directive: str) -> bool:
    """Return True if the named directive appears at the top level of the process body.

    Scans the masked source so that directive-shaped text embedded in a
    quoted string cannot produce a false positive. Excludes matches that
    land inside a script/shell/exec section.
    """
    body_start = proc.body_start + 1
    body_end = proc.body_end
    pattern = rf"(?m)^[ \t]*{re.escape(directive)}\b"
    for m in re.finditer(pattern, parsed.masked[body_start:body_end]):
        abs_off = body_start + m.start()
        in_script = any(
            sec and sec[0] <= abs_off < sec[1] for sec in (ps.script, ps.shell, ps.exec)
        )
        if not in_script:
            return True
    return False


# The ``ProcessSections`` symbol is imported lazily to avoid a circular
# import with the module's header; we re-expose it here.
from .parser import ProcessSections  # noqa: E402


@rule("NF050", Severity.ERROR, "module process must declare a container directive")
def nf050_process_environment(parsed: ParsedFile) -> Iterable[Issue]:
    """Every module process needs an explicit container image.

    SKILL.md: the unit of reproducibility is the container image, not a
    conda solve at run time. A ``container`` directive is mandatory. A
    ``conda`` spec is optional alongside it (Wave can build an image from
    the spec and cache it), but conda alone is not sufficient.
    """
    if not parsed.is_module:
        return
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        if _directive_exists(parsed, proc, ps, "container"):
            continue
        yield _make_issue(
            parsed,
            proc.header_start,
            "NF050",
            Severity.ERROR,
            f"process '{proc.name}' declares no container directive — "
            "every module needs an explicit image (use Wave to build one from conda if needed)",
        )


@rule("NF051", Severity.WARNING, "pipeline should define params via nextflow_schema.json")
def nf051_schema_file(parsed: ParsedFile) -> Iterable[Issue]:
    """Entry-level pipelines should have a JSON schema defining their params.

    SKILL.md: "Use schema-driven validation for user-facing pipelines."
    The schema doubles as documentation for every tunable parameter and
    is what tools like ``nextflow run -params-file`` and the nf-core
    launcher read to validate user input before any task is submitted.
    """
    if not parsed.is_main_nf:
        return
    schema = parsed.path.parent / "nextflow_schema.json"
    if not schema.exists():
        yield _make_issue(
            parsed,
            0,
            "NF051",
            Severity.WARNING,
            "no nextflow_schema.json next to main.nf — add one so params get "
            "schema-driven validation and self-documentation",
        )


_RESOURCE_LABEL_RE = re.compile(
    r"(?m)^[ \t]*label\b",
)
_RESOURCE_DIRECTIVES = ("cpus", "memory", "time", "machineType", "accelerator")


@rule("NF052", Severity.WARNING, "module process must declare compute-cost controls")
def nf052_process_resources(parsed: ParsedFile) -> Iterable[Issue]:
    """Every module needs a handle for operators to tune compute cost.

    The nf-core convention is a ``label 'process_low'`` (or similar) on
    each module; the pipeline's ``conf/base.config`` then maps labels to
    concrete ``cpus`` / ``memory`` / ``time`` / ``machineType`` values.
    Alternatively the module can declare explicit resources directly.
    Either way, a module with no resource signal leaves cloud operators
    unable to scale or constrain the task without editing the module.
    """
    if not parsed.is_module:
        return
    for proc in _processes(parsed):
        ps = process_sections(parsed.masked, proc)
        body_start = proc.body_start + 1
        body_end = proc.body_end
        body_masked = parsed.masked[body_start:body_end]

        # A ``label 'process_...'`` directive anywhere outside a script
        # section satisfies the rule.
        has_label = False
        for m in _RESOURCE_LABEL_RE.finditer(body_masked):
            abs_off = body_start + m.start()
            if not any(
                sec and sec[0] <= abs_off < sec[1] for sec in (ps.script, ps.shell, ps.exec)
            ):
                has_label = True
                break
        if has_label:
            continue

        # Or an explicit resource directive.
        has_resource = any(_directive_exists(parsed, proc, ps, d) for d in _RESOURCE_DIRECTIVES)
        if has_resource:
            continue

        yield _make_issue(
            parsed,
            proc.header_start,
            "NF052",
            Severity.WARNING,
            f"process '{proc.name}' declares no 'label' or resource directives "
            "(cpus/memory/time/machineType) — operators cannot tune compute cost from config",
        )
