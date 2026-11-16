<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Sample Sheet Schema Specification

Sample sheet schema files are used by the nf-schema plugin for validation of sample sheet contents and type conversion / channel generation.

Based on JSON Schema 2020-12. See [JSON Schema docs](https://json-schema.org/understanding-json-schema).

## Schema Structure

Sample sheets are typically CSV/TSV files parsed as arrays of objects. The schema must match this structure:

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "field_1": { "type": "string" },
      "field_2": { "type": "string" }
    }
  }
}
```

- Top-level type is `array`
- `items` contains a single `object`
- `properties` keys must match CSV headers
- Fields in sample sheet but not in schema produce warnings
- **Property order in schema defines output channel order**

## Common Keys

Most keys are identical to the Nextflow schema specification:
`type`, `pattern`, `format`, `errorMessage`, `exists`, `enum`, `default`, `deprecated`

See the Nextflow schema specification reference for details.

## Sample Sheet-Specific Keys

### `meta`

Type: `List` or `String`

Marks a field as a meta value for channel generation. Should contain the meta field name(s).

```json
{
  "sample": {
    "type": "string",
    "meta": ["id"]
  }
}
```

This converts the field value to a meta map entry: `[[id:value]...]`

## Example: RNA-seq Samplesheet

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/nf-core/rnaseq/master/schema_input.json",
  "title": "nf-core/rnaseq pipeline - params.input schema",
  "description": "Schema for the file provided with params.input",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "sample": {
        "type": "string",
        "pattern": "^\\S+$",
        "errorMessage": "Sample name must be provided and cannot contain spaces",
        "meta": ["my_sample"]
      },
      "fastq_1": {
        "type": "string",
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "format": "file-path",
        "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
      },
      "fastq_2": {
        "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'",
        "type": "string",
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "format": "file-path"
      },
      "strandedness": {
        "type": "string",
        "errorMessage": "Strandedness must be provided and be one of 'forward', 'reverse' or 'unstranded'",
        "enum": ["forward", "reverse", "unstranded"],
        "meta": ["my_strandedness"]
      }
    },
    "required": ["sample", "fastq_1", "strandedness"]
  }
}
```

## Example: Feature-Complete Test Schema

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/nextflow-io/nf-schema/master/src/testResources/schema_input.json",
  "description": "Schema for sample sheet",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "sample": {
        "type": "string",
        "pattern": "^\\S+$",
        "errorMessage": "Sample name cannot contain spaces",
        "meta": ["id"]
      },
      "fastq_1": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$"
      },
      "fastq_2": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$"
      },
      "strandedness": {
        "type": "string",
        "enum": ["forward", "reverse", "unstranded"]
      }
    },
    "required": ["sample", "fastq_1"]
  }
}
```

## Design Tips

1. **Order matters** — property order defines channel output order
2. **Use `meta`** — map sample identifiers to Nextflow meta maps
3. **Validate paths** — use `format: "file-path"` + `exists: true` for input files
4. **Error messages** — add `errorMessage` to every field for user-friendly validation
5. **Required fields** — list in `items.required`, not per-property
6. **Patterns** — validate file extensions and disallow whitespace
