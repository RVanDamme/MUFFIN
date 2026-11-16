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


# nf-validation to nf-schema Migration

Migrate a user's pipeline from nf-validation and JSON Schema draft-07 to nf-schema v2 and JSON Schema 2020-12. Work with the pipeline and schema files the user supplies; do not assume where they are stored.

## Identify the Legacy Format

Inspect nextflow_schema.json and any additional parameter or samplesheet schemas supplied by the user. Legacy signals include draft-07, definitions, references to definitions, property-level unique, and the old dependentRequired representation.

Inspect the user's Nextflow source for imports from plugin/nf-validation, Channel.fromSamplesheet, manual paramsHelp or paramsSummaryLog calls, and validation settings stored as flat params.validation values. Inspect nextflow.config for the old plugin declaration.

## Migration Workflow

1. Inventory the supplied schemas, source imports, samplesheet conversion, help behavior, and validation configuration.
2. Replace the nf-validation plugin declaration and imports with an appropriate nf-schema v2 release.
3. Change each schema to JSON Schema 2020-12, rename definitions to $defs, and update the corresponding $ref values.
4. Migrate special keywords deliberately; simple text replacement is insufficient.
5. Move validation settings from pipeline parameters into the validation configuration scope.
6. Replace fromSamplesheet with samplesheetToList while preserving the channel's contents and expected tuple shape.
7. Update automatic help configuration and remove obsolete manual help calls.
8. Validate the schema and run the user's existing representative pipeline test.

Do not apply global replacement commands across an assumed checkout. Edit only the files identified in the inventory.

## Plugin Declaration

~~~groovy
plugins {
    id 'nf-schema@2.x.x'
}
~~~

Replace the version placeholder with the release appropriate for the user's pipeline. Imports identify the public plugin:

~~~groovy
include { samplesheetToList } from 'plugin/nf-schema'
~~~

## Schema Changes

Use the public schema URI https://json-schema.org/draft/2020-12/schema.

For complete-object uniqueness, use array-level uniqueItems. For uniqueness of selected columns or fields, use nf-schema's uniqueEntries after confirming the current specification. Do not convert field uniqueness into whole-row uniqueness accidentally.

Move dependentRequired to the object containing the dependent properties. For example:

~~~json
{
  "dependentRequired": {
    "fastq_2": ["fastq_1"]
  }
}
~~~

## Validation Configuration

| Legacy parameter | Validation setting |
|---|---|
| validationMonochromeLogs | monochromeLogs |
| validationLenientMode | lenientMode |
| validationFailUnrecognisedParams | failUnrecognisedParams |
| validationShowHiddenParams | showHiddenParams |
| validationIgnoreParams | defaultIgnoreParams or ignoreParams, according to intent |

Place these settings inside the validation scope of the user's Nextflow configuration. defaultIgnoreParams expresses pipeline-defined exclusions; ignoreParams expresses user-selected exclusions.

For help features supported by the selected nf-schema release, configure validation.help and remove obsolete manual help code.

## Samplesheet Conversion

samplesheetToList returns a list. Construct the channel from that list using the user's actual samplesheet and schema references. Confirm required fields, types, field order, and metadata construction remain correct.

## Completion Checklist

- The selected nf-schema release supports the user's Nextflow version.
- Plugin imports and configuration are updated.
- All supplied schemas use the intended draft and $defs references.
- Special uniqueness and dependency semantics are preserved.
- Samplesheet conversion produces the expected channel contract.
- Help behavior is correct.
- The user's representative validation and pipeline test pass.

Consult [migration details](references/migration-guide.md) for legacy keyword examples. This is a historical reference: verify current nf-core tooling compatibility rather than treating its dated support warning as current.
