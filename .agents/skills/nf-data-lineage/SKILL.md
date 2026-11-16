---
name: nf-data-lineage
description: >
  Explore Nextflow data lineage to trace where outputs came from and what inputs
  produced them. Uses the configured lineage store and `nextflow lineage`
  CLI. Requires Nextflow 25.04+ with lineage.enabled=true. Use when the user asks
  about provenance, data origins, reproducibility, or output traceability.
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


# Nextflow Data Lineage

Trace the provenance of pipeline outputs — what inputs, parameters, and processes
produced a given result.

## When to Use

Load this skill when the user wants to:
- Know where a specific output file came from
- Trace back from a result to the original inputs
- Understand the full provenance chain of a dataset
- Check if two outputs were produced by the same run/parameters
- Verify reproducibility of results
- Explore what the lineage store contains

## Prerequisites

Data lineage requires **Nextflow 25.04+** with lineage enabled:

```groovy
// nextflow.config
lineage {
    enabled = true
    store {
        // Set location only if the user chooses a nondefault metadata store.
    }
}
```

Check if lineage is available:

```bash
# Check Nextflow version
nextflow -version

```

Inspect the selected run's effective configuration to confirm lineage is enabled,
then verify that its configured store is accessible in that environment.

If lineage is not enabled, guide the user to enable it and re-run their pipeline.

## Lineage Concepts

### Lineage IDs (LIDs)

Every tracked artifact gets a unique Lineage ID:
- Format: `lid://<hash>` — content-addressable identifier
- LIDs reference: workflow runs, tasks, outputs, inputs, parameters

### What Gets Tracked

- **Workflow runs** — which pipeline ran, when, with what parameters
- **Task executions** — each process invocation with inputs/outputs
- **Data artifacts** — files produced and consumed, with checksums
- **Parameters** — the exact params used for each run

### Lineage Store

The configured lineage store contains:
- Metadata entries indexed by LID
- Relationships between artifacts (input → task → output)

## Analysis Procedure

### Step 1: Check lineage availability

Identify the lineage store configured for the selected run and verify that it is
accessible through the available host tools.

Read `lineage.store.location` from the run's effective configuration rather than
assuming a store location.

### Step 2: Use the `nextflow lineage` CLI

```bash
# List recent lineage entries
nextflow lineage list

# Show details of a specific entry
nextflow lineage log <lid>

# Trace the provenance of an output file
nextflow lineage render <lid>
```

### Step 3: Trace an output back to its source

When the user asks "where did this file come from?":

```bash
# Find the LID for an output file
nextflow lineage list | grep "<filename>"

# Get the full provenance chain
nextflow lineage log <lid>
```

This shows:
- Which task produced the file
- What inputs the task consumed
- Which workflow run it belonged to
- What parameters were used

### Step 4: Compare runs

To check if outputs from different runs used the same inputs/parameters:

```bash
# Get lineage for both outputs
nextflow lineage log <lid_1>
nextflow lineage log <lid_2>

# Compare the parameter sets and input checksums
```

### Step 5: Browse the lineage store directly

If the CLI doesn't provide enough detail:

Inspect the selected store's metadata entries with the host's available tools.
Use only entries associated with the user's requested run or output.

## Response Guidelines

**Do:**
- Explain provenance as a chain: "This VCF was produced by GATK_HAPLOTYPECALLER,
  which consumed a BAM from BWA_MEM, which aligned the original FASTQ files you
  provided as input"
- Mention checksums when relevant for reproducibility
- Note the parameters that were active for a given run
- Help the user understand whether a result can be reproduced
- If lineage isn't enabled, explain how to enable it and what they'll get

**Don't:**
- Dump raw lineage store JSON without interpretation
- Assume lineage is available — always check first
- Confuse lineage (provenance tracking) with the work directory (execution cache)
- Show LID hashes without context — always pair with human-readable descriptions

## Example Response Style

> Your file selected variant file was produced in your run from
> yesterday (run name: `focused_borg`, 2025-03-04 14:32).
>
> The provenance chain:
> - **Input**: `sample_A_R{1,2}.fastq.gz` (checksums match your original files)
> - **TRIM_GALORE** removed adapters → trimmed reads
> - **BWA_MEM** aligned to GRCh38 (the `--genome` param was set to `GRCh38`)
> - **GATK_HAPLOTYPECALLER** called variants with default settings
>
> This was produced with the same parameters as your March 2nd run, but with
> updated input files — the R1 FASTQ checksum differs between the two runs.

## Enabling Lineage

If the user hasn't enabled lineage yet:

```groovy
// Add to nextflow.config
lineage {
    enabled = true
}
```

Then re-run the pipeline. Lineage metadata will be generated for all new runs.
Previous runs without lineage enabled cannot be retroactively tracked.

## Related Skills

- `debug-local-run` — for debugging failures (uses work dir, not lineage store)
- [narrative run-history analysis](../nextflow-history/references/nf-run-history/README.md) — for run history narrative (uses run-history metadata)
- `nextflow-config` — for configuring lineage settings
