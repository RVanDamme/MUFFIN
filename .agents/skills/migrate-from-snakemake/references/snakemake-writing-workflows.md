<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
<!-- Purpose: Migration-focused notes from official Snakemake writing/workflow grammar docs. -->

# snakemake-writing-workflows

Use this when inventorying a Snakefile before converting to Nextflow.

## Why this source helps

The official grammar gives a reliable checklist of top-level constructs that must be mapped during migration:

- `rule`
- `include`
- `module` / `use rule`
- `configfile`
- `workdir`
- execution blocks (`shell`, `script`, `run`, `notebook`)
- runtime directives (`threads`, `resources`, `conda`, `container`, `log`, `params`)

This reduces missed features during initial translation.

## Migration use

1. Parse top-level statements first.
2. Build a translation inventory:
   - `include`/`module`/`use rule` → DSL2 subworkflows/modules
   - `configfile`/global settings → `nextflow.config` + `params`
   - per-rule directives → process directives and script sections
3. Confirm minimum Snakemake version assumptions (`min_version`) before inferring behavior.

## Sources

- ReadTheDocs: https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html
- Upstream source (authoritative): https://github.com/snakemake/snakemake/blob/main/docs/snakefiles/writing_snakefiles.rst
