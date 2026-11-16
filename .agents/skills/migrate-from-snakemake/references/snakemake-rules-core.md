<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
<!-- Purpose: Migration-focused notes from official Snakemake rules documentation. -->

# snakemake-rules-core

Primary reference for rule semantics to preserve during Snakemake → Nextflow migration.

## High-value concepts for migration

### Rule anatomy
- `input`, `output`, `params`, `log`, `threads`, `resources`, `shell`/`script`/`run`
- Map to Nextflow `process` input/output blocks and process directives.

### Wildcards and ambiguity
- Wildcards drive DAG resolution from output patterns.
- Constrained wildcards (`wildcard_constraints`) and ambiguous-rule behavior are critical when rewriting path logic.
- Migration implication: make tuple/meta channels explicit; avoid hidden filename inference.

### Aggregation/input helpers
- `expand`, `collect`, input functions, `lookup`, `branch` encode fan-in and dynamic dependencies.
- Migration implication: map to channel transforms (`map`, `flatMap`, `collect`, `combine`, etc.) and explicit workflow wiring.

### Compute directives
- `threads` and `resources` define scheduling semantics, including dynamic resources.
- Migration implication: preserve intent with `cpus`, `memory`, `time`, labels, and executor config.

### Execution form
- `script`, `wrapper`, `conda`, `container` tell you whether to reuse nf-core modules vs local modules.

## Migration checklist from this doc

- Capture wildcard patterns and constraints per rule.
- Capture resource/thread semantics per rule.
- Classify each rule body: shell/script/wrapper/run.
- Identify `rule all`/default target behavior.

## Sources

- ReadTheDocs: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
- Upstream source (authoritative): https://github.com/snakemake/snakemake/blob/main/docs/snakefiles/rules.rst
