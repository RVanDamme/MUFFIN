---
name: nextflow-history
description: >
  Analyze local Nextflow run history and cache. Use when user asks about their
  recent local runs, wants to see what pipelines they've executed, understand
  run lineage, inspect the Nextflow cache, or correlate runs with work directories.
  Uses run history, task cache metadata, and Nextflow logs supplied by the user.
  Also use for narrative recaps, recurring failures and pipeline evolution.
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


# Nextflow Local History & Cache

Use this skill to explain local Nextflow run history, task cache behavior, and
storage use from metadata associated with the user's selected run environment.
Use the available `nextflow log` interface or the files the user provides; do
not assume a working directory or cache location.

## When to use

- Summarize recent runs, their status, duration, pipeline, and command.
- Explain which cached tasks can be reused with `-resume`.
- Compare activity or storage use across selected runs.
- Identify safe cleanup candidates after reviewing the user's retention needs.

For failure diagnosis, use the `debug-local-run` skill. For Platform history,
use MCP discovery/execution to retrieve workflows; read [MCP connection](../launch-workflow/references/seqera-mcp/README.md) if connection or discovery guidance is needed.

## Run-history fields

Nextflow's run history describes each launch with these fields:

| Field | Meaning |
|---|---|
| Timestamp | When the run started |
| Duration | Wall-clock duration |
| Run name | Generated or user-selected run identifier |
| Status | Whether the run completed successfully |
| Hash | Cache identifier recorded for the run |
| Session ID | Identifier linking the run to its task cache |
| Command | Pipeline and launch arguments |

Read the metadata for the requested environment. Summaries can group by
pipeline, status, date, or duration without prescribing where the metadata lives.

## Built-in history interface

If Nextflow is available in the selected environment, use its history interface:

```bash
nextflow log
nextflow log <run-name>
nextflow log <run-name> -f 'name,status,hash,duration,realtime,%cpu,rss'
nextflow log <run-name> -f name -F 'status == "COMPLETED"'
nextflow log <run-name> -f hash
```

Replace `<run-name>` with an identifier obtained from the actual history.
Run these operations in the environment that owns the selected history rather
than assuming local access also provides access to remote runs.

## Log entries

For versions supporting `nextflow logfile`, use the selected run identifier
and requested severity or entry count:

```bash
nextflow logfile <run-name> -level ERROR
nextflow logfile <run-name> -level WARN
nextflow logfile <run-name> -n 10
```

Otherwise inspect the log supplied by the user. Preserve multiline exceptions
and correlate entries with the requested run. Report unavailable history or
logs explicitly rather than inferring successful execution.

## Task cache and resume

A session's task cache associates task hashes with results. Run indexes link
the task records to a run identifier. With `-resume`, Nextflow reuses eligible
results when the task's inputs, script, and other hashed settings still match
and the required output files remain accessible.

Inspect the actual cache and work files associated with the selected run.
Compare task identifiers, generated commands, and storage usage before deciding
whether to resume or clean up. Preserve data needed for current runs and future
resume operations.

Use the available `nextflow clean` interface to review cleanup candidates.
Check its dry-run behavior in the installed version before deletion, and act
only within the user's authorized retention and cleanup scope.

## Reporting

For a narrative recap or recurring patterns, read [narrative run-history analysis](references/nf-run-history/README.md). First establish the selected environment, time range and completeness of the records. Use logs/trace evidence to distinguish confirmed causes from correlations.

Lead with the observed run count, success/failure status, and requested time
range. Include the selected pipeline versions and identifiers when useful.
Explain cache reuse or cleanup implications from the actual evidence. Distinguish
historical failures from current failures and avoid presenting stale logs as
the latest run.

Related skills: `nextflow-config`, `debug-local-run`,
`debug-seqera-failed-run`, and [MCP connection](../launch-workflow/references/seqera-mcp/README.md).

## Narrative history analysis

When the user wants a recap, recurring failure patterns or how a pipeline evolved over time, analyze the selected records and explain gaps without inferring causality from timing alone. Read [narrative run-history analysis](references/nf-run-history/README.md) before proceeding.
