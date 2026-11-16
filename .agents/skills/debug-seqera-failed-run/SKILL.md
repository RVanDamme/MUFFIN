---
name: debug-seqera-failed-run
description: >
  Debug failed Nextflow Platform pipeline runs. Fetches workflow details, failed
  tasks, and logs to identify root causes and suggest fixes. Use when the user
  asks why a Platform run failed, to debug the last run, to diagnose a
  Seqera workflow error, or to triage recurring failures across runs and
  propose a fix for review.
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


# Debug Seqera Failed Run

Systematic analysis of Nextflow Platform pipeline failures using MCP tools.

## When to Use

Load this skill when the user wants to:
- Debug why their last Platform run failed
- Analyze failed tasks on Nextflow Platform
- Understand resource, storage, credential, configuration, or input-data issues
- Get actionable fixes for a workflow error

## Core Approach

1. Resolve the workflow and workspace.
2. Fetch workflow details/progress before task logs.
3. Inspect failed tasks and their logs.
4. Use workflow-level logs for launch-time, staging, credential, or cascade
   failures.
5. Extract the failing process and first real error line to form a failure
   fingerprint, especially when many tasks failed.
6. If multiple attempts exist for this run/session, compare them before
   recommending a fix.
7. Identify the root cause, not just the first visible symptom.
8. Classify the failure type and provide an actionable fix.

## Progressive References

- Read `debug-workflow.md` when you need exact Platform MCP calls or
  a step-by-step investigation checklist.
- Read `failure-patterns.md` when logs are noisy, multiple tasks
  failed, or you need to distinguish infrastructure failures from pipeline,
  configuration, or bioinformatics/input-data errors.
- Read `resume-and-attempt-history.md` when the run was resumed, the
  user mentions retrying, or multiple attempts exist for the same session/run.
- Read `fingerprints-and-cascades.md` when many tasks failed, several
  errors mention the same path/artifact, or downstream tasks may be cascade
  symptoms of an upstream failure.
- Read `recurring-failures-and-fix-handoff.md` when several runs
  failed, the user asks which failures recur, the pipeline is a customized
  fork of nf-core, or the user wants a fix proposed as a branch, pull
  request, or ticket.
- Read `fix-confidence.md` when logs are incomplete, the category is
  ambiguous, or the answer needs to distinguish proven fixes from hypotheses.

## Diagnosis Principles

- Prefer explicit log evidence over exit-code guesses.
- Use category-first triage: launch/config, storage/staging/auth, transient
  compute, resource exhaustion, container/registry, then bioinformatics/tool.
- Distinguish infrastructure/resource failures from bioinformatics/tool/input
  errors.
- For cascades, report the original missing artifact, permission issue, or
  launch/config validation failure rather than blaming downstream tasks.
- Quote short evidence excerpts only; avoid dumping verbose logs back to the
  user.
- State confidence when evidence is incomplete; ask for `.command.err`, task
  stderr/stdout, or workflow logs instead of inventing a fix.
- Treat successful resumed-attempt deltas as strong clues, but do not overclaim
  causality unless the changed field directly matches the failed log evidence.

## Response Format

1. **Run summary** — pipeline name, run name, status, duration, tasks succeeded/failed
2. **Root cause** — the specific error with context and failure classification
3. **Failed tasks** — which process(es) failed, with exit codes and process tags
4. **Evidence** — relevant log excerpts, kept brief
5. **Attempt history** — if a later resumed attempt succeeded, summarize what
   changed between the failed and successful attempts
6. **Fix** — concrete next steps such as retry/resume, more memory/disk, IAM or
   data-link changes, sample sheet correction, or pipeline/config changes

## External References

- [Nextflow Platform troubleshooting](https://docs.seqera.io/platform/latest/troubleshooting)
- [Nextflow process error strategies](https://www.nextflow.io/docs/latest/process.html#error-strategy)
