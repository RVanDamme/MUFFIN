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


# Nextflow Pipeline Debugger

Run supported diagnostic passes and report evidence. Apply fixes only when the user's request authorizes repair, then re-run the affected checks. Check the installed Nextflow version/help before using version-dependent commands.

## Diagnostic Steps

### 1. Lint

```bash
nextflow lint <pipeline>
```

Catches syntax errors, deprecated patterns, undeclared variables.

### 2. Config Validation

```bash
nextflow config -show-profiles
```

Lists available profiles. Inspect the merged configuration for the selected pipeline/profile separately to diagnose unresolved parameters, profile conflicts and duplicate settings. Profile listing alone does not prove the intended configuration is correct.

### 3. Preview (DAG compilation)

```bash
nextflow run <pipeline> -preview
```

Compiles the pipeline and generates the DAG without executing. Surfaces channel wiring errors, missing process inputs, and type mismatches.

## Workflow

1. Run all three commands, collect output.
2. Analyze errors — group by file and severity.
3. Report the diagnosis; if repair was requested, make the smallest evidence-based source edit.
4. Re-run the failing diagnostic to confirm the fix.
5. Report which supported checks passed, failed or could not run. For a repair, run representative tests/output checks as well; static success is not runtime verification.

## Common Fixes

| Symptom | Likely cause |
|---------|-------------|
| `Missing process` | Typo in `include` or workflow block |
| `Channel has no origin` | Process output not connected |
| `Invalid config` | Unquoted string or missing closure brace |
| `Deprecated operator` | Replace with current DSL2 equivalent |

## Notes

- Assumes `nextflow` is on PATH.
- If the repo has a custom entry script (not `main.nf`), adjust paths in `nextflow run`.
- When additional context is provided (e.g. specific errors), focus diagnostics on those areas first.
