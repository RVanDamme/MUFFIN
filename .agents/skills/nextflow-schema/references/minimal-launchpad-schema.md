<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Minimal Launchpad Schema Guidance

Use this guide when the user wants a reduced Seqera Launchpad form, especially for cloud launches.

## What "minimal" means

- **Canonical schema** (`nextflow_schema.json`) remains the full source of truth.
- **Minimal schema** (`nextflow_schema.minimal.json`) is an optional reduced UI projection selected for Launchpad; the filename is an example.
- The minimal file should not replace, rename, or weaken the canonical schema.

Creating `nextflow_schema.minimal.json` alone does not change Launchpad. Select
the file's repository path in Platform or set Launchpad's `schemaName` to that
path. Alternatively, upload it through the documented Platform schema flow.
See [Launchpad configuration](https://docs.seqera.io/platform-cloud/launch/launchpad)
and [pipeline schema selection](https://docs.seqera.io/platform-cloud/pipeline-schema/overview).

## When to create `nextflow_schema.minimal.json`

Create a minimal schema when users ask for one of these:
- "Keep only a few launch parameters in Launchpad"
- "Hide advanced options for cloud users"
- "Simplify Launchpad form for templates"

Do **not** create a minimal schema when users only ask for correctness/validation updates to the pipeline schema.

## Authoring workflow

1. Generate/update `nextflow_schema.json` first.
2. Select a minimal subset of parameters for Launchpad UX.
3. Save that subset as `nextflow_schema.minimal.json`.
4. Keep labels/descriptions/defaults aligned with canonical schema.
5. Select the reduced schema's repository path/`schemaName`, or upload it in Platform.
6. Revisit minimal schema whenever canonical launch-critical defaults change.

## Cloud Launchpad pitfalls

### Keep `input` and `outdir` explicit

For cloud launches, minimal schemas usually still need both:
- `input` (samplesheet or input manifest)
- `outdir` (cloud-safe destination for outputs)

Hiding `outdir` often causes users to keep stale or local defaults accidentally.

### Minimal schema does not clear saved Launchpad defaults

Reducing visible fields does not reset values stored in saved Launchpad launches/templates.
If users still see an unexpected output destination, they usually need to edit/reset saved defaults.

### Keep schema defaults and Launchpad saved defaults aligned

If schema defaults change but saved Launchpad defaults do not, users can launch with unexpected values.
When updating schema defaults, explicitly call out that saved defaults may need a refresh.

### `tower.yml` is not runtime output configuration

`tower.yml` report/timeline/trace settings affect monitoring artifacts.
They do not replace runtime output configuration (`params.outdir`, cloud storage destinations, work directory choices).

## Minimal schema template (example)

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/<org>/<pipeline>/master/nextflow_schema.minimal.json",
  "title": "<org>/<pipeline> launchpad minimal parameters",
  "type": "object",
  "$defs": {
    "launchpad_required": {
      "title": "Launchpad required options",
      "type": "object",
      "required": ["input", "outdir"],
      "properties": {
        "input": {
          "type": "string",
          "format": "file-path",
          "schema": "schema_input.json",
          "description": "Samplesheet input for the run."
        },
        "outdir": {
          "type": "string",
          "format": "directory-path",
          "description": "Cloud output directory for results."
        }
      }
    }
  },
  "allOf": [{ "$ref": "#/$defs/launchpad_required" }]
}
```

## Review checklist

- Canonical schema still complete and up to date.
- Minimal schema exists only when user requested form reduction.
- Minimal schema includes cloud launch-critical fields (`input`, `outdir`).
- Platform/Launchpad explicitly selects the reduced schema; the filename alone is not enough.
- Team understands that Launchpad saved defaults may need manual reset.
- Notes distinguish `tower.yml` reporting from runtime output path settings.
