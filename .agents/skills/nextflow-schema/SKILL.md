---
name: nextflow-schema
description: >
  Generate canonical nextflow_schema.json, optional nextflow_schema.minimal.json,
  and sample sheet schema files for Nextflow pipelines. Use when user asks to create,
  update, or validate a pipeline parameter schema or sample sheet schema. Covers type
  inference, parameter grouping, nf-core conventions, Launchpad minimal-form guidance,
  nf-schema v2 (JSON Schema 2020-12) compliance, source-material parameter triage,
  and migration from nf-validation or draft-07 schemas.
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


# Generating Nextflow Schema

Generate canonical and Launchpad-oriented minimal schema files for Nextflow pipelines.

## When to Use

Load this skill when the user wants to:
- Generate a `nextflow_schema.json` from a `nextflow.config`
- Generate a reduced Launchpad form schema (`nextflow_schema.minimal.json`)
- Create or update a sample sheet schema (`schema_input.json`)
- Add parameters to an existing schema
- Validate schema structure
- Migrate schema from draft-07 to 2020-12 (nf-schema v2)
- Identify which parameters to expose when given raw source material (scripts, notebooks, papers, CLI docs)

## Discovering parameters from raw source material

Read [parameter triage](references/triage-pipeline-parameters/README.md) for the bounded source-audit task, including the config outline and representative test-input manifest.

When given source material (scripts, notebooks, paper methods sections, CLI docs) rather than an existing `nextflow.config`, scan exhaustively before writing the schema. Missing a parameter now means a breaking change to surface it later.

### What to scan

- CLI flags and defaults in the tools being wrapped
- Constants and magic numbers hard-coded in scripts or notebook cells
- Configuration files the original authors loaded (YAML, TOML, INI, JSON)
- Per-function keyword arguments with non-trivial defaults
- Paths to reference data (genomes, checkpoints, blast databases)
- Numbers or thresholds discussed in the paper, README, or method section — if the authors mention a value in prose, it matters

### Triage each collected value

| Class | When to use | Example |
|-------|-------------|---------|
| **Surface as parameter** | User would reasonably change it: inputs, outputs, thresholds, tool/algorithm choices, reference paths, any flag whose value the paper discusses | k-mer size, model checkpoint, significance threshold |
| **Pin as internal constant** | Part of the tool's own contract, never touched: fixed random seeds for reproducibility, internal file format versions | `random_seed = 42` |
| **Expose as profile override** | Environment-dependent: executor, container engine, test-input paths | `executor = 'slurm'` |

When in doubt, surface. Promoting a constant to a parameter later is a breaking change; removing an unused parameter is trivial.

### Triaged parameter audit table

Before writing the schema, produce a markdown table with the source location of each parameter so the reviewer can audit the call:

| Name | Source | Default | Type | Decision | Rationale |
|------|--------|---------|------|----------|-----------|
| `msa_kmer_size` | `run_msa.sh:12` | 15 | integer | **surface** | Paper §2.3 varies this per organism |
| `random_seed` | `fold.py:4` | 42 | integer | **pin** | Fixed for reproducibility |
| `executor` | (n/a) | local | string | **profile** | Environment-dependent |

Then use this table as input to the schema generation workflow below.

## References

Use the host's file-reading tool to load detailed docs:
- `schema-specification.md` — Full nf-schema JSON Schema specification
- `minimal-launchpad-schema.md` — When and how to generate `nextflow_schema.minimal.json`
- `sample-sheet-specification.md` — Sample sheet schema specification and examples

## Workflow

1. **Choose target output** — Decide whether the user needs canonical only, or canonical + minimal Launchpad schema.
2. **Extract params** — Parse all `params {}` from `nextflow.config`
3. **Determine types** — Infer JSON Schema types from Nextflow defaults
4. **Write descriptions** — Short `description`, longer `help_text` where needed
5. **Group parameters** — Organize into logical `$defs` groups
6. **Mark required** — Params with `null` default are typically required
7. **Add validation** — Patterns, enums, formats, min/max constraints
8. **Validate** — Ensure valid JSON Schema 2020-12
9. **If minimal requested** — Generate `nextflow_schema.minimal.json` from canonical schema and keep both aligned

## Full vs Minimal Schema Generation

### Canonical schema (`nextflow_schema.json`)

Use this as the source of truth for:
- Full parameter validation and docs
- nf-schema tooling expectations
- Pipeline maintenance over time

### Minimal Launchpad schema (`nextflow_schema.minimal.json`)

Use this only when the user explicitly wants a reduced Seqera Launchpad form, especially for cloud launches.

Rules:
1. Generate and maintain `nextflow_schema.json` first; never replace it with minimal.
2. Create `nextflow_schema.minimal.json` as a curated subset for UX.
3. Keep launch-critical fields explicit in minimal schemas, especially `input` and `outdir`.
4. Keep defaults and descriptions aligned between canonical and minimal schemas.

`nextflow_schema.minimal.json` is an example filename, not an automatically
discovered Launchpad convention. Select its repository path in Platform or set
Launchpad's `schemaName` to that path. Alternatively, upload the reduced schema
through the documented Platform schema flow. Creating the file alone does not
change the form. See [Launchpad configuration](https://docs.seqera.io/platform-cloud/launch/launchpad)
and [pipeline schema selection](https://docs.seqera.io/platform-cloud/pipeline-schema/overview).

If the user asks for "minimal schema", "fewer Launchpad fields", or "cloud launch form cleanup", load `minimal-launchpad-schema.md`.

## Launchpad / Cloud Pitfalls

- A minimal schema does **not** fix stale Launchpad saved defaults. If users have existing launch presets, they may still inject old values (for example a previously selected output destination).
- `outdir` should remain explicit in minimal cloud schemas so users can set cloud-safe destinations.
- `tower.yml` report/timeline/trace settings are not a replacement for runtime output configuration (`params.outdir`, work/output paths).
- Schema updates and Launchpad saved defaults must stay in sync; otherwise UI defaults can differ from schema defaults.

## Schema Structure

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/<org>/<pipeline>/master/nextflow_schema.json",
  "title": "<org>/<pipeline> pipeline parameters",
  "description": "Pipeline description",
  "type": "object",
  "$defs": {
    "input_output_options": {
      "title": "Input/output options",
      "type": "object",
      "fa_icon": "fas fa-terminal",
      "description": "Define where the pipeline should find input data and save output data.",
      "required": ["input", "outdir"],
      "properties": {
        "input": {
          "type": "string",
          "format": "file-path",
          "exists": true,
          "mimetype": "text/csv",
          "pattern": "^\\S+\\.(csv|tsv|yaml|json)$",
          "schema": "schema_input.json",
          "description": "Path to input samplesheet.",
          "fa_icon": "fas fa-file-csv"
        },
        "outdir": {
          "type": "string",
          "format": "directory-path",
          "description": "Output directory for results.",
          "fa_icon": "fas fa-folder-open"
        }
      }
    }
  },
  "allOf": [
    { "$ref": "#/$defs/input_output_options" }
  ]
}
```

## Type Mapping

Infer JSON Schema types from Nextflow config defaults:

| Nextflow Default | JSON Schema Type | Notes |
|-----------------|-----------------|-------|
| `null` | `"string"` | Usually required; add `format` if path-like name |
| `'string'` / `"string"` | `"string"` | Check for enum hints in comments |
| `42` (integer) | `"integer"` | Add `minimum`/`maximum` if sensible |
| `3.14` (float) | `"number"` | Add `minimum`/`maximum` if sensible |
| `true` / `false` | `"boolean"` | |
| `'128.GB'` / `'24.h'` | `"string"` | Memory/time — string type, not numeric |
| `[:]` (Map) | `"object"` | Empty map default |
| `[]` (List) | `"array"` | Empty list default |
| `'25.MB'` | `"string"` | Nextflow memory unit string |

### Special Type Inference

- Param named `*_input`, `*_fasta`, `*_gtf`, `*_bed`, `*_bam` → `"format": "file-path"`
- Param named `*_outdir`, `*_dir` → `"format": "directory-path"`
- Param named `*email*` → `"format": "email"`
- Comments with `// enum:` or `// DEPRECATED:` → extract enum values or add `"deprecated": true`
- Param named `input` with null default → add `"schema": "schema_input.json"` if sample sheet pipeline

## Parameter Grouping

### Standard nf-core Groups

| Group Key | Title | fa_icon | Contents |
|-----------|-------|---------|----------|
| `input_output_options` | Input/output options | `fas fa-terminal` | input, outdir, email, multiqc_title |
| `reference_genome_options` | Reference genome options | `fas fa-dna` | genome, fasta, gtf, gff, index paths, igenomes_* |
| `institutional_config_options` | Institutional config options | `fas fa-university` | custom_config_*, hostnames (all hidden) |
| `max_job_request_options` | Max job request options | `fab fa-acquisitions-incorporated` | max_cpus, max_memory, max_time (all hidden) |
| `generic_options` | Generic options | `fas fa-file-import` | validate_params, monochrome_logs, etc. (most hidden) |

### Pipeline-Specific Groups

Create additional groups for pipeline-specific params:
- **Trimming options** — skip_trimming, min_trimmed_reads, etc.
- **Alignment options** — aligner, pseudo_aligner, etc.
- **Quality control** — skip_fastqc, skip_multiqc, etc.

### Grouping Rules

1. Every param MUST be in a `$defs` group (no ungrouped top-level properties)
2. Every group MUST have a corresponding `allOf` `$ref` entry
3. Group keys use `snake_case`
4. Group titles use Title Case with "options" suffix
5. Internal/developer params → separate group with `hidden` params
6. AWS/cloud-specific params → separate conditional group

## Required Parameters

- Params with `null` default are candidates for `required`
- `input` is almost always required
- `outdir` is almost always required
- List required params in the group's `required` array, NOT per-param

## Common fa_icon Values

| Parameter Type | Icon |
|---------------|------|
| Input files | `fas fa-file-csv`, `fas fa-file`, `far fa-file-code` |
| Output dirs | `fas fa-folder-open` |
| Genome/DNA | `fas fa-dna` |
| Boolean toggles | `fas fa-toggle-on` |
| Numbers/counts | `fas fa-calculator` |
| Email | `fas fa-envelope` |
| Config/settings | `fas fa-cog`, `fas fa-users-cog` |
| Cloud/S3 | `fas fa-cloud-download-alt` |
| Hidden/internal | `fas fa-eye-slash` |
| Tools/processes | `fas fa-tools` |

## Validation Keys

Apply these where appropriate:

```json
// String with regex pattern
{ "type": "string", "pattern": "^\\S+\\.csv$" }

// Enum values
{ "type": "string", "enum": ["star_salmon", "star_rsem", "hisat2"] }

// File path with existence check
{ "type": "string", "format": "file-path", "exists": true }

// Number with range
{ "type": "integer", "minimum": 1, "maximum": 128 }

// Deprecated param
{ "type": "string", "deprecated": true, "errorMessage": "Use --email instead" }

// Hidden param
{ "type": "string", "hidden": true }
```

## Sample Sheet Schema

For `schema_input.json`:

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/<org>/<pipeline>/master/assets/schema_input.json",
  "title": "<pipeline> samplesheet schema",
  "description": "Schema for the samplesheet used by <pipeline>",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "sample": {
        "type": "string",
        "pattern": "^\\S+$",
        "errorMessage": "Sample name must be provided and cannot contain spaces",
        "meta": ["id"]
      },
      "fastq_1": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "errorMessage": "FastQ file for reads 1 must match '*.fq.gz' or '*.fastq.gz'"
      },
      "fastq_2": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$"
      }
    },
    "required": ["sample", "fastq_1"]
  }
}
```

### Sample Sheet Key Differences
- Top-level `type: "array"` with `items` containing `type: "object"`
- Property order in schema defines channel output order
- `meta` key maps fields to Nextflow meta map (e.g., `"meta": ["id"]`)
- Fields in sample sheet but not in schema produce warnings

## Generation Guidelines

1. **Always use 2020-12** — `"$schema": "https://json-schema.org/draft/2020-12/schema"`
2. **Use `$defs`** not `definitions` (nf-schema v2)
3. **Every param in a group** — no ungrouped params
4. **Keep descriptions short** — one sentence max; use `help_text` for details
5. **Match config defaults** — schema `default` must match `nextflow.config`
6. **Omit default key** when value is null — don't set `"default": null`
7. **Add fa_icon** to every param and group for UI navigation
8. **Hidden params** — institutional config, internal params, igenomes_base
9. **Validate enum from comments** — `// enum: a, b, c` → `"enum": ["a", "b", "c"]`
10. **Extract deprecation** — `// DEPRECATED:` → `"deprecated": true` + `errorMessage`
11. **For minimal schema requests** — keep `input` + `outdir` explicit and document that Launchpad saved defaults may need reset
12. **Do not confuse `tower.yml` with runtime output** — reporting config does not set cloud result destinations

## Migration

For migrating existing pipelines from nf-validation to nf-schema, read [schema migration](references/nf-schema-migration/README.md) before editing. Schema syntax validation alone does not prove samplesheet or pipeline equivalence.

## Source-material parameter triage

When deriving parameters from scripts, notebooks, papers or CLI documentation, produce the source audit, schema, config outline and representative test-input manifest. Read [parameter triage](references/triage-pipeline-parameters/README.md) before proceeding.

## Legacy schema migration

When upgrading nf-validation or draft-07 schemas, preserve validation, samplesheet tuple shape, help behavior and representative execution. Read [schema migration](references/nf-schema-migration/README.md) before proceeding.
