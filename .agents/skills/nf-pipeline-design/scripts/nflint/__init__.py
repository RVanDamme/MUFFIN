"""nflint — an opinionated Nextflow linter for the design-nextflow-pipelines skill.

Enforces the structural rules described in SKILL.md: thin main.nf, well-formed
subworkflows with take/emit and version channels, atomic modules that do not
reach into params, and strict-syntax hygiene.
"""

from .issues import Issue, Severity
from .parser import ParsedFile, parse_file

__version__ = "0.1.0"
__all__ = ["Issue", "Severity", "parse_file", "ParsedFile"]
