---
name: debug-local-run
description: >
  Debug a local Nextflow pipeline run. Analyzes .nextflow.log, work directories,
  and task error logs to identify failures and suggest fixes. Use when user asks
  to debug their last run, diagnose pipeline errors, or understand why a run failed.
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


# Debug a Local Nextflow Run

Diagnose the user's run from actual logs, task evidence, active configuration,
and history. Use files and task locations identified by the user or returned by
the run metadata. Do not assume a fixed working directory or search layout.

## Procedure

1. Identify the intended run and its log. Ask for the relevant evidence when the
   host cannot access it. Common log filenames include `.nextflow.log`, but
   the user may have configured another location.
2. Read the run summary: pipeline, revision, status, timing, and completed or
   failed task counts.
3. Find the first meaningful failure and identify the affected task. Use the
   task location returned by the logs or history rather than inventing it.
4. Inspect the actual command, stdout, stderr, exit status, and task trace when
   available. Nextflow evidence filenames may include `.command.sh`,
   `.command.log`, `.command.err`, `.command.out`, `.exitcode`, and
   `.command.trace`.
5. Compare requested and observed CPU, memory, disk, and runtime. Check for OOM,
   scheduler kills, quota limits, timeouts, or staging failures.
6. Inspect the user's active configuration, including `nextflow.config` and
   any configuration it includes. Check profile, executor, container, process
   settings, and credentials appropriate to the selected environment.

## Log commands

When supported by the installed Nextflow version, `nextflow logfile` filters
multi-line entries without splitting stack traces:

```bash
nextflow logfile -level ERROR
nextflow logfile -level WARN
nextflow logfile -n 20
```

These commands apply only when the host has the intended run context. Otherwise,
read the user-provided log with the available file tools. For older Nextflow
versions, inspect the final error entries and the lines containing ERROR, WARN,
Exception, or Caused by.

## Compile errors

For `Script compilation error`, `MultipleCompilationErrorsException`,
`No such variable`, `No such property`, `MissingMethodException`,
`Unknown process directive`, or `Unexpected input`:

1. Inspect the named file and line first.
2. Trace the failing channel or value back through its immediate construction.
   Check tuple shapes and destructuring around branch, combine, mix, and join.
3. Make the smallest change supported by the error.
4. Validate the correction with `nextflow lint` or the user's selected test
   execution. Report the new error precisely if it moves.

Do not edit a downstream process that has not run unless the compile evidence
proves it is involved.

## Common failures

| Evidence | Likely category | Next step |
|---|---|---|
| Exit status 137 with OOM evidence | Memory exhaustion | Adjust process memory and Java heap consistently |
| Exit status 139 | Segmentation fault | Check tool input, version, and container |
| Exit status 1 | Tool or script failure | Read the meaningful stderr before diagnosing |
| Missing file | Input or staging failure | Check the user's input references and channel logic |
| Permission denied | Access or container user mismatch | Verify permissions on the actual task storage |
| Command not found | Environment dependency | Check the container, executable availability, and module setup |
| Disk quota exceeded | Storage exhaustion | Check actual storage usage and required capacity |
| Connection refused | Network or executor issue | Check endpoints, configuration, and credentials |

Distinguish a root cause from downstream cascades. Do not treat the first warning
or a generic wrapper message as the explanation.

## Report

Provide the run summary, failing process and exit status, relevant evidence,
diagnosis with confidence, and concrete correction. State which checks were
actually run and which remain unverified.

References:
[Nextflow error handling](https://www.nextflow.io/docs/latest/process.html#error-strategy)
and [Nextflow tracing](https://www.nextflow.io/docs/latest/tracing.html).
