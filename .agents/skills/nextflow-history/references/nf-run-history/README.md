<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../../../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Nextflow Run History Analysis

Analyze the user's available Nextflow execution history and explain patterns, progress, recurring issues, and changes over time.

## Gather the Available Evidence

Ask which local or remote execution environment the user wants analyzed. Use history records, Nextflow logs, and trace reports that the user provides or that the host can locate in that environment. Do not assume a particular directory or that local records describe remote runs.

A Nextflow history record can include a timestamp, duration, run name, status, run hash, session UUID, and complete invocation. Logs explain failures; trace reports add task-level resource and timing evidence when tracing was enabled.

If no history is available, explain what evidence is needed and continue with any supplied logs or reports.

## Analysis Workflow

1. Establish the time range, pipelines, number of runs, and completeness of the records.
2. Compare success and failure patterns, command changes, profiles, parameters, resume usage, and durations.
3. Inspect logs for recent or recurring failures. Tie suspected causes to the actual errors and invocations.
4. Distinguish correlations from confirmed causes. A configuration change near an improvement does not prove causation.
5. Summarize the trajectory and current state in natural language.

For Nextflow releases with logfile support, use the command's run-name and level filters:

The bundled baseline introduces this command in Nextflow 26.05.0-edge. Check
the installed version's help before using it.

~~~bash
nextflow logfile -level WARN
~~~

A particular run name can be supplied when known. For older releases, inspect the log records the user provided using available host tools.

## Report Style

Lead with the overall activity, then explain the arc: early failures, what changed, turning points, and current stability. Cite a few relevant run names or timestamps rather than listing every record.

Use task and resource evidence to make specific suggestions. For example, repeated memory failures in one process justify examining that process's resource request before raising every task's allocation.

Do not dump the entire history, show uninterpreted logs, or claim causes that the evidence does not support. Explain gaps when records are incomplete.
