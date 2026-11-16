---
name: migrate-from-snakemake
description: >
  Convert Snakemake workflows to idiomatic Nextflow DSL2 pipelines with a
  test-first, modular migration loop. Trigger: "snakemake to nextflow",
  "migrate from snakemake", "convert snakefile", "/migrate-from-snakemake".
---
<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Migrate from Snakemake to Nextflow

Convert a Snakemake pipeline to Nextflow **incrementally**. Do not attempt a one-shot rewrite.

## When to Use

Load this skill when the user wants to:
- Convert Snakemake to Nextflow
- Port a `Snakefile` / `*.smk` workflow to DSL2
- Replace Snakemake wrappers/rules with Nextflow modules
- Preserve behavior using tests and golden outputs while migrating

## Migration Strategy (Human-like, test-first)

### 0) Detect nf-core compliance requirements first

If the pipeline contains `.nf-core.yml`, the user requests nf-core compatibility, or `manifest.name = "nf-core/<name>"` is required, preserve the nf-core scaffolding and validate it with `nf-core pipelines lint` alongside the migrated workflow.

Bootstrap these files before any process exists:

- `.nf-core.yml` — pipeline metadata + lint config (`repository_type: pipeline`, `nf_core_version`)
- `nextflow.config` — `manifest.name = "nf-core/<pipeline>"`, `manifest.version`, `manifest.nextflowVersion = "!>=24.04.0"`, plus a `params {}` block
- Resource and retry configuration — `process.errorStrategy`, resource defaults keyed by labels (`process_low`, `process_medium`, `process_high`)
- Per-process publication configuration — `withName:` `publishDir` blocks for every process that emits files
- Subworkflow metadata — each subworkflow declares inputs/outputs and emits a software-versions channel
- Per-process `versions.yml` emission so subworkflows can collect software versions

Check the relevant lint rules, including `nfcore_yml`, `nextflow_config`, `base_config`, `modules_config`, `local_component_structure`, `pipeline_todos`, `template_strings`, and `merge_markers`. Run `nf-core pipelines lint -k <keys> --json <file>` from the pipeline root before declaring the migration done.

For the complete conventions (branches, edit policy, CLI commands), chain `maintain-nf-core-pipeline` and read its `nf-core-conventions.md`.

### 1) Establish source-of-truth before edits

0. **Audit readiness first** — run `audit-conversion-readiness` over
   the Snakefile and its config. Every `config[...]` key with no configfile entry,
   every `script:`/`wrapper:` that is not in the handover, and every reference path
   on an unreachable mount is a blocker that stops a baseline run before it starts.
1. **Visualize the DAG first** — run `snakemake --rulegraph | dot -Tsvg > rulegraph.svg` to see the rule dependency structure at a glance. Use `--dag` for the full job-level DAG or `--filegraph` for file-level dependencies.
2. Find and run existing tests first (if present).
3. Run the Snakemake pipeline end-to-end on a small fixture dataset.
4. Save expected outputs (checksums, row counts, key reports) as migration baseline.
5. If no tests exist, create baseline regression checks from this run.

Treat this as your spec. The rulegraph is your conversion roadmap — each node becomes a Nextflow process, each edge becomes a channel connection.

### 2) Convert tests early

- Convert Snakemake behavior checks into `nf-test` suites:
  - **Spec tests**: expected behavior for core pipeline paths
  - **Regression tests**: known bug reproductions and fixes
- Start with failing tests when reproducing migration bugs.

### 3) Convert top-level workflow first

- Use the rulegraph from step 1 as the blueprint — translate each rule node to a process, each edge to a channel.
- Translate the main `Snakefile` into `main.nf` + `nextflow.config` skeleton.
- Map global config and top-level DAG flow before implementing every rule detail.
- Keep process stubs minimal but executable so the workflow compiles early.
- For non-trivial param surfaces: [parameter triage](../nextflow-schema/references/triage-pipeline-parameters/README.md) before locking `nextflow.config`.

### 4) Decompose by responsibility

- Included `*.smk` files → **subworkflows**
- Snakemake wrappers / rule bodies → **modules**
- Keep one logical operation per module.

### 5) Prefer nf-core modules as starting points

For each wrapper/rule tool, call `search-existing-modules` first — it returns a verdict (reuse nf-core, vendor a community module, or canonical invocation for a project-owned module). Patch only the delta for local behavior.

### 6) Wire dataflow explicitly

- Snakemake wildcards become channel values (often tuple/meta patterns).
- `expand()` fan-in patterns map to `.collect()` / `.toList()` where required.
- Keep inputs/outputs explicit in every process.

### 7) Iterate in small verified chunks

For each converted piece:
1. Make one focused change.
2. Run targeted `nf-test`.
3. Run `nextflow run <pipeline> -preview`.
4. Run `nextflow lint`.
5. Re-run fixture workflow and compare against Snakemake baseline.

## Snakemake → Nextflow Mapping

| Snakemake | Nextflow DSL2 |
|---|---|
| `rule x:` | `process X {}` |
| `input:` / `output:` | `input:` / `output:` blocks |
| Wildcards (`{sample}`) | Channel values / tuple metadata |
| `SAMPLES = [...]` + wildcards | Samplesheet CSV + `splitCsv` + tuple/meta channels |
| `include: "foo.smk"` | Import the selected module or subworkflow using the user's existing project organization |
| `wrapper:` | nf-core module or local module |
| `conda:` / `container:` | `conda` / `container` directives (verify versions for target platform) |
| `threads:` | `cpus` config + `task.cpus` |
| `resources:` | process directives (`memory`, `time`, `disk`) |
| `configfile:` | `params` in `nextflow.config` |
| `rule all:` targets | workflow outputs / reachable DAG |
| `script:` using a user-provided helper | Refactored script with CLI args (replace `snakemake.input` → `sys.argv`) |
| Pre-built indexes (assumed) | Explicit INDEX process before alignment |
| `snakemake --dag/--rulegraph` | `nextflow run <pipeline> -preview` (DAG visualization) |

## Common Pitfalls

### Conda/container version resolution

Snakemake tutorials and wrappers often pin old package versions that may not exist on all platforms.

- **Never copy version pins blindly.** Verify packages exist for the target platform (`osx-arm64`, `linux-64`) before running.
- Run `conda search -c bioconda <pkg>` to check availability.
- Prefer current stable versions over old tutorial pins (e.g., `bioconda::bwa=0.7.19` not `0.7.17`).
- Bioconda ARM64 coverage lags behind x86 — if a version is missing, try the latest available or fall back to a container.
- When using containers, prefer BioContainers mulled images that bundle related tools (e.g., bwa + samtools).

### Implicit prerequisites

Snakemake workflows often assume pre-built artifacts exist without a rule to create them. These **must become explicit Nextflow processes**:

- **Aligner indexes** (BWA, Bowtie2, STAR, HISAT2) — add an INDEX process that runs before alignment. Wire its output to the alignment process input.
- **Genome dicts/faidx** (`samtools faidx`, `samtools dict`, `picard CreateSequenceDictionary`) — add explicit processes if downstream tools require `.fai` or `.dict` files.
- **Database files** (BLAST DBs, Kraken2 DBs) — add download/build processes or accept as params.
- **Clue:** any file referenced in `input:` that no rule produces is an implicit prerequisite.

### Script block conversion

Snakemake `script:` blocks use magic objects (`snakemake.input`, `snakemake.output`, `snakemake.params`, `snakemake.threads`). These don't exist in Nextflow — convert to CLI arguments:

- `snakemake.input[0]` → positional arg (`sys.argv[1]`) or named arg (`argparse`)
- `snakemake.output[0]` → positional arg or hardcoded output name in the process script
- `snakemake.threads` → `task.cpus` (passed via process script template)
- `snakemake.params.foo` → named arg or environment variable

### Sample handling

- Snakemake `SAMPLES = ["A", "B"]` + expand/wildcards → Nextflow samplesheet CSV + `splitCsv` + tuple/meta channels.
- Prefer samplesheet input (`params.samplesheet`) over hardcoded sample lists for reproducibility and flexibility.
- Use `[id: row.sample]` meta maps to carry sample identity through the channel graph.

## Delivery Checklist

- [ ] Baseline Snakemake run captured on fixture data
- [ ] Tests exist (spec + regression)
- [ ] Main workflow compiles (`nextflow run <pipeline> -preview`)
- [ ] Rule groups split into subworkflows/modules
- [ ] nf-core module reuse attempted before custom modules
- [ ] Lint passes (`nextflow lint`)
- [ ] Migrated outputs compared to Snakemake baseline
- [ ] Migration notes document known deltas

## References

- `snakemake-nextflow-playbook.md` — practical migration command patterns and module sourcing flow
- `snakemake-writing-workflows.md` — Snakefile grammar checklist for migration inventory
- `snakemake-rules-core.md` — authoritative rule semantics (wildcards/resources/aggregation)
- `snakemake-testing-unit-tests.md` — baseline test generation and migration-to-nf-test guidance
- `nextflow-conda-troubleshooting.md` — conda resolution failures, platform checks, container fallbacks
