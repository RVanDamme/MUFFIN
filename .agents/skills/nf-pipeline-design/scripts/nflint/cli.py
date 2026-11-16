"""Command-line entry point for nflint.

Usage:
    nflint [paths...]                 # lint given files/directories
    nflint --format json <paths>      # machine-readable output
    nflint --select NF001,NF010 ...   # run only specific rules
    nflint --ignore NF030,NF031 ...   # skip specific rules
    nflint --list-rules               # show all registered rules

Exit codes:
    0 — clean
    1 — at least one issue with severity >= min-severity (default: warning)
    2 — internal error (e.g. unreadable file)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .issues import Issue, Severity
from .parser import parse_file
from .rules import all_rules

SEVERITY_ORDER = {Severity.INFO: 0, Severity.WARNING: 1, Severity.ERROR: 2}


def discover_files(targets: list[Path]) -> list[Path]:
    files: list[Path] = []
    for target in targets:
        if target.is_file() and target.suffix == ".nf":
            files.append(target)
        elif target.is_dir():
            files.extend(sorted(target.rglob("*.nf")))
    # Deterministic order
    return sorted(dict.fromkeys(files))


def lint_file(path: Path, active_codes: set[str]) -> list[Issue]:
    parsed = parse_file(path)
    issues: list[Issue] = []
    for rule in all_rules():
        if active_codes and rule.code not in active_codes:
            continue
        for issue in rule.fn(parsed):
            issues.append(issue)
    issues.sort(key=lambda i: (i.path, i.line, i.column, i.code))
    return issues


def filter_by_severity(issues: list[Issue], min_severity: Severity) -> list[Issue]:
    threshold = SEVERITY_ORDER[min_severity]
    return [i for i in issues if SEVERITY_ORDER[i.severity] >= threshold]


def parse_csv(value: str | None) -> set[str]:
    if not value:
        return set()
    return {x.strip() for x in value.split(",") if x.strip()}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="nflint",
        description=(
            "Opinionated linter for the design-nextflow-pipelines skill. Enforces "
            "the structural rules of SKILL.md plus strict-syntax hygiene."
        ),
    )
    ap.add_argument("paths", nargs="*", type=Path, help="Files or directories to lint")
    ap.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Output format (default: text)",
    )
    ap.add_argument(
        "--select",
        default="",
        metavar="CODES",
        help="Comma-separated list of rule codes to run (e.g. 'NF001,NF010')",
    )
    ap.add_argument(
        "--ignore",
        default="",
        metavar="CODES",
        help="Comma-separated list of rule codes to skip",
    )
    ap.add_argument(
        "--min-severity",
        choices=("info", "warning", "error"),
        default="warning",
        help="Only fail on issues at or above this severity (default: warning)",
    )
    ap.add_argument("--list-rules", action="store_true", help="Print all rules and exit")
    ap.add_argument("--root", type=Path, default=None, help="Display paths relative to this root")
    return ap


def cmd_list_rules() -> int:
    for r in all_rules():
        print(f"{r.code:>6}  {r.severity.value:<7}  {r.title}")
    return 0


def run(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.list_rules:
        return cmd_list_rules()

    if not args.paths:
        print("nflint: no paths given; pass a file or directory", file=sys.stderr)
        return 2

    selected = parse_csv(args.select)
    ignored = parse_csv(args.ignore)
    active = {r.code for r in all_rules()}
    if selected:
        active &= selected
    active -= ignored

    files = discover_files(args.paths)
    if not files:
        print("nflint: no .nf files found", file=sys.stderr)
        return 2

    all_issues: list[Issue] = []
    for f in files:
        try:
            all_issues.extend(lint_file(f, active))
        except OSError as exc:
            print(f"nflint: {f}: {exc}", file=sys.stderr)
            return 2

    min_sev = Severity(args.min_severity)
    failing = filter_by_severity(all_issues, min_sev)

    root = args.root
    if args.format == "json":
        payload = {
            "summary": {
                "files": len(files),
                "issues": len(all_issues),
                "failing": len(failing),
            },
            "issues": [i.to_dict(root) for i in all_issues],
        }
        print(json.dumps(payload, indent=2))
    else:
        for i in all_issues:
            print(i.format_text(root))
        n_err = sum(1 for i in all_issues if i.severity == Severity.ERROR)
        n_warn = sum(1 for i in all_issues if i.severity == Severity.WARNING)
        print(
            f"\n{len(all_issues)} issue(s): {n_err} error, {n_warn} warning — "
            f"scanned {len(files)} file(s)",
            file=sys.stderr,
        )

    return 1 if failing else 0


def main() -> None:
    sys.exit(run())


if __name__ == "__main__":
    main()
