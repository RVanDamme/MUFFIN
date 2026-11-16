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


# Triage Pipeline Parameters and Draft the Schema

Given scripts, a repository, notebook, methods, or tool documentation supplied by the user, identify the parameter surface and representative test inputs. Do not assume where the source or generated artifacts are stored.

Produce four artifacts:

1. A triaged parameter list.
2. A draft nextflow_schema.json.
3. A companion nextflow.config configuration outline.
4. A representative test-input manifest.

Use the user's chosen output destination. Do not write workflow source or design channels during this parameter-triage task.

## What to Inspect

Collect values that affect observable results: tool flags and defaults, constants, configuration values, reference resources, model settings, thresholds, and choices discussed in methods.

For each value, record its actual source reference, default, type, effect, and uncertainty. Do not invent a source filename or location.

## Triage

- Surface input resources, output choices, result-affecting thresholds, algorithms, and settings the user would reasonably vary as top-level parameters.
- Pin tool-internal conventions and justified reproducibility constants when they should not be exposed.
- Put executor, container-engine, and environment-dependent defaults in profiles rather than biological parameters.

Explain every decision and keep ambiguous cases visible for review.

## Draft the Schema

Follow the current nf-schema specification:
https://nextflow-io.github.io/nf-schema/latest/nextflow_schema/nextflow_schema_specification/

For nf-schema v2, use JSON Schema 2020-12, $defs, and corresponding $ref values. Group parameters by user intent, such as inputs/outputs, analysis settings, tool arguments, references, and general options.

Each property needs a type and clear description. Add an appropriate default or required membership. Add enums, numeric bounds, and supported file/directory formats when they describe the actual contract.

~~~json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "object",
  "$defs": {
    "analysis_options": {
      "type": "object",
      "properties": {}
    }
  },
  "allOf": [
    {"$ref": "#/$defs/analysis_options"}
  ]
}
~~~

The JSON Pointer is part of the schema contract, not a local filesystem location.

## Configuration Outline

Mirror applicable schema defaults in the user's Nextflow configuration and identify required plugin and profile settings. Use the existing project's configuration conventions without inventing include locations, local output folders, or test data filenames.

Separate script parameters from configuration-only settings when strict syntax or typed params requires it. Use meaningful parameter names and established nf-core conventions when they fit the analysis.

## Representative Test Inputs

For each heavy tool, identify an existing user-provided fixture or explain how to derive a small representative input. Record the input's verified reference, format, size, biological context, and expected output behavior.

A test input must exercise the real algorithm while remaining practical for validation. If no suitable input is available, state the gap rather than generating a fixed file location.

## Return Format

Present the four artifacts with a short uncertainty note. Use a table for parameter decisions and test-input provenance. Explain missing defaults, ambiguous effects, naming choices, and unavailable representative inputs so the caller can review them.
