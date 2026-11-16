---
name: nextflow-output-patterns
description: >
  Use this skill when writing or reviewing any Nextflow workflow that
  aggregates per-sample outputs into a combined file, joins channels with
  a remainder mode, or groups tuples for per-group processing. Covers the
  operator idioms — `collectFile`, `join`, `groupTuple`, channel-level
  null handling — that determine whether a workflow's output is
  *correct*, not just whether it *runs*. Trigger on phrases like
  `collectFile`, `groupTuple`, `join`, channel aggregation, per-sample
  merging, combined output files, metadata-keyed grouping, optional
  tuple positions, and the broader category of "the workflow finishes
  without error but the merged output is malformed." Pair with
  migrate-nextflow-code (strict v2 syntax) and build-nextflow-pipeline
  (high-level structure and module composition); this skill is specifically about
  output-correctness idioms in the operator layer.
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


# Nextflow output patterns

You are writing or reviewing a Nextflow workflow that aggregates per-sample outputs. Three operator idioms commonly produce subtle output bugs that the agent gets wrong on first attempt. Apply them by default.

## 1. `collectFile(seed:)` for the header

When concatenating per-sample TSV/CSV outputs into one file, **`collectFile` does NOT prepend a header automatically**. A downstream consumer will see the first data row as the header and the file will be off by one line.

### ❌ DO NOT — header missing

```groovy
ch_summaries
    .map { sample, summary -> summary.text }
    .collectFile(name: 'merged.tsv', newLine: true)
```

The output has N rows. The expected output has N+1 rows (header + N data). Off-by-one. Every consumer of this file is wrong.

### ✅ DO — explicit `seed:` header

```groovy
def header = 'sample_id\tgenome\tparsed_fcsgx\ttiara_summary\tadaptor_file'

ch_summaries
    .map { sample, summary -> summary.text }
    .collectFile(name: 'merged.tsv', seed: header, newLine: true)
```

`seed:` writes the given string as the FIRST line of the collected file. Subsequent items are appended after. `newLine: true` ensures the seed line and each item end with a newline.

### Edge case — empty input channel

When the upstream channel is empty (e.g. no samples for an optional output), `collectFile` produces **no output file at all** — not an empty file, not a header-only file. If a downstream consumer requires the file to exist, pair with `.ifEmpty([])` upstream, or assert the channel is non-empty before `collectFile`.

## 2. Null-coalescing after `join(remainder: true)`

`join(remainder: true)` keeps tuples whose join key has no match by filling missing positions with `null`. If you then format the tuple as a tab-separated string, **literal `null` strings leak into the output**.

### ❌ DO NOT — null leakage

```groovy
ch_primary
    .join(ch_optional, remainder: true)
    .map { sample, primary, optional -> "${sample}\t${primary}\t${optional}" }
    // → 'sampleB\tfull-genome-b\tnull'  ← literal 'null' string in the file
```

### ✅ DO — coalesce nulls before formatting

```groovy
ch_primary
    .join(ch_optional, remainder: true)
    .map { values ->
        values.collect { v -> v != null ? v : '' }.join('\t')
    }
```

This pattern works for any-arity tuple. The `collect { ... }` coalesces every position; replace `''` with whatever sentinel your downstream tool expects (`'NA'`, `'.'`, `'-'`, etc.).

### When NOT to use `remainder: true`

If every key MUST be in both channels, drop `remainder: true` and let `join` filter the mismatches. Use `remainder: true` only when you explicitly want partial-match rows in the output, and pair it with null-coalescing every time.

## 3. `groupTuple` aggregation by metadata key

When aggregating per-replicate or per-condition outputs to a group level (e.g. all replicates of one sample), use `groupTuple` keyed on the metadata position, not `collect()` over the full channel.

### ✅ DO — group then aggregate

```groovy
ch_per_replicate    // tuple val(meta), path(file)
    .map { meta, file -> tuple(meta.sample, file) }
    .groupTuple()    // tuple val(sample), list<path>
    .map { sample, files -> tuple(sample, files.sort()) }   // deterministic order
```

`groupTuple` blocks until all upstream items with the same key arrive, then emits a list. Always `.sort()` the list afterward for reproducible output ordering — Nextflow does not guarantee arrival order.

### Edge case — `groupTuple(by: N)` for non-zero key position

If the grouping key is not the first tuple element, specify `by:`:

```groovy
ch_files    // tuple val(sample), val(condition), path(file)
    .groupTuple(by: 1)   // group by condition
```

## Checklist when reviewing aggregation code

- [ ] Every `collectFile` that produces a TSV/CSV with a header has `seed:` set.
- [ ] Every `join(remainder: true)` is followed by a null-coalescing `map { ... }`.
- [ ] Every `groupTuple` aggregation calls `.sort()` on the grouped list.
- [ ] Channels that may be empty either have `.ifEmpty(...)` upstream or the workflow asserts non-emptiness explicitly.

## When this skill is NOT enough

- For `splitFasta`, `splitText`, `splitCsv` (the *input* side of file-channel work), consult the Nextflow docs (the host's web tools restricted to official https://docs.seqera.io/nextflow/ documentation).
- For `collectFile`'s `storeDir:` / `cache:` / `sort:` parameters, consult the current official Nextflow operator reference.
- For metadata-channel shape design (when meta has multiple keys), use `nf-pipeline-design`.
