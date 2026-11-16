<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Template merge guides

One file per `nf-core/tools` release that introduces a template-sync
PR. Each file is a distilled merge-conflict checklist for that
release.

## Lookup

Read the tools version from the open template-sync PR title (e.g.
*"Template update for nf-core/tools v4.0.0"*). Then read the
corresponding file in this directory:

| Tools version | Guide |
|---|---|
| 4.0.x | `tools-4_0.md` |

If no file exists for the tools version in the PR:

1. Fetch the matching blog post at `https://nf-co.re/blog` for the required release.
2. Resolve conflicts following the patterns described there.
