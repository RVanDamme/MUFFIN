---
description: Schema specification for Nextflow pipeline configuration parameter validation
---
<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->

# Nextflow schema specification

The Nextflow schema file contains information about pipeline configuration parameters.
The file is typically saved in the workflow root directory and called `nextflow_schema.json`.

The Nextflow schema syntax is based on the JSON schema standard, with some key differences.
You can find more information about JSON Schema here:

- Official docs: <https://json-schema.org>
- Excellent "Understanding JSON Schema" docs: <https://json-schema.org/understanding-json-schema>

!!! warning

    This file is a reference specification, _not_ documentation about how to write a schema manually.

    Please see the `nextflow-schema` skill for instructions on how to create these files
    (and don't be tempted to do it manually in a code editor!)

!!! note

    The nf-schema plugin, as well as several other interfaces using Nextflow schema, uses a stock JSON schema library for parameter validation.
    As such, any valid JSON schema _should_ work for validation.

    However, please note that graphical UIs (docs, launch interfaces) are largely hand-written and may not expect JSON schema usage that is
    not described here. As such, it's safest to stick to the specification described here and _not_ the core JSON schema spec.

## Definitions

A slightly strange use of a JSON schema standard that we use for Nextflow schema is `$defs`.

JSON schema can group variables together in an `object`, but then the validation expects this structure to exist in the data that it is validating.
In reality, we have a very long "flat" list of parameters, all at the top level of `params.foo`.

In order to give some structure to log outputs, documentation and so on, we group parameters into `$defs`.
Each `def` is an object with a title, description and so on.
However, as they are under `$defs` scope they are effectively ignored by the validation and so their nested nature is not a problem.
We then bring the contents of each definition object back to the "flat" top level for validation using a series of `allOf` statements at the end of the schema,
which reference the specific definition keys.

<!-- prettier-ignore-start -->
```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "object",
  // Definition groups
  "$defs": { // (1)!
    "my_group_of_params": { // (2)!
      "title": "A virtual grouping used for docs and pretty-printing",
      "type": "object",
      "required": ["foo", "bar"], // (3)!
      "properties": { // (4)!
        "foo": { // (5)!
          "type": "string"
        },
        "bar": {
          "type": "string"
        }
      }
    }
  },
  // Contents of each definition group brought into main schema for validation
  "allOf": [
    { "$ref": "#/$defs/my_group_of_params" } // (6)!
  ]
}
```
<!-- prettier-ignore-end -->

1. An arbitrary number of definition groups can go in here - these are ignored by main schema validation.
2. This ID is used later in the `allOf` block to reference the definition.
3. Note that any required properties need to be listed within this object scope.
4. Actual parameter specifications go in here.
5. Shortened here for the example, see below for full parameter specification.
6. A `$ref` line like this needs to be added for every definition group

Parameters can be described outside of the `$defs` scope, in the regular JSON Schema top-level `properties` scope.
However, they will be displayed as ungrouped in tools working off the schema.

## Nested parameters

!!! example "New feature in v2.1.0"

Nextflow config allows parameters to be nested as objects, for example:

```groovy
params {
    foo {
        bar = "baz"
    }
}
```

or on the CLI:

```bash
nextflow run <pipeline> --foo.bar "baz"
```

Nested parameters can be specified in the schema by adding a `properties` keyword to the root parameters:

```json
{
  "type": "object",
  "properties": {
    "thisIsNested": {
      // Annotation for the --thisIsNested parameter
      "type": "object", // Parameters that contain subparameters need to have the "object" type
      "properties": {
        // Add other parameters in here
        "deep": {
          // Annotation for the --thisIsNested.deep parameter
          "type": "string"
        }
      }
    }
  }
}
```

There is no limit to how deeply nested parameters can be. Mind however that deeply nested parameters are not that user friendly and will create some very ugly help messages. It's advised to not go deeper than two levels of nesting.

## Required parameters

Any parameters that _must_ be specified should be set as `required` in the schema.

!!! tip

    Make sure you do set `null` as a default value for the parameter, otherwise it will have a value even if not supplied by the pipeline user and the required property will have no effect.

This is not done with a property key like other things described below, but rather by naming
the parameter in the `required` array in the definition object / top-level object.

For more information, see the [JSON schema documentation](https://json-schema.org/understanding-json-schema/reference/object.html#required-properties).

```json
{
  "type": "object",
  "properties": {
    "name": { "type": "string" },
    "email": { "type": "string" },
    "address": { "type": "string" },
    "telephone": { "type": "string" }
  },
  "required": ["name", "email"]
}
```

## Parameter name

The `properties` object _key_ must correspond to the parameter variable name in the Nextflow config.

For example, for `params.foo`, the schema should look like this:

```json
// ..
"type": "object",
"properties": {
    "foo": {
        "type": "string",
        // ..
    }
}
// ..
```

## Keys for all parameters

### `type`

Variable type, taken from the [JSON schema keyword vocabulary](https://json-schema.org/understanding-json-schema/reference/type.html):

- `string`
- `number` (float)
- `integer`
- `boolean` (true / false)
- `object` (currently only supported for file validation, see [Nested parameters](#nested-parameters))
- `array` (currently only supported for file validation, see [Nested parameters](#nested-parameters))

Validation checks that the supplied parameter matches the expected type, and will fail with an error if not.

This JSON schema type is _not_ supported:

- `null`

### `default`

Default value for the parameter.

Should match the `type` and validation patterns set for the parameter in other fields.

!!! tip

    If no default should be set, completely omit this key from the schema.
    Do not set it as an empty string, or `null`.

    However, parameters with no defaults _should_ be set to `null` within your Nextflow config file.

!!! note

    When creating a schema using `nf-core schema build`, this field will be automatically created based
    on the default value defined in the pipeline config files.

    Generally speaking, the two should always be kept in sync to avoid unexpected problems and usage errors.
    In some rare cases, this may not be possible (for example, a dynamic groovy expression cannot be encoded in JSON),
    in which case try to specify as "sensible" a default within the schema as possible.

### `description`

A short description of what the parameter does, written in markdown.
Printed in docs and terminal help text.
Should be maximum one short sentence.

### `help_text`

!!! example "Non-standard key"

A longer text with usage help for the parameter, written in markdown.
Can include newlines with multiple paragraphs and more complex markdown structures.

Typically hidden by default in documentation and interfaces, unless explicitly clicked / requested.

### `errorMessage`

!!! example "Non-standard key"

If validation fails, an error message is printed to the terminal, so that the end user knows what to fix.
However, these messages are not always very clear - especially to newcomers.

To improve this experience, pipeline developers can set a custom `errorMessage` for a given parameter in a the schema.
If validation fails, this `errorMessage` is printed after the original error message to guide the pipeline users to an easier solution.

For example, instead of printing:

```
* --input (samples.yml): "samples.yml" does not match regular expression [^\S+\.csv$]
```

We can set

```json
"input": {
  "type": "string",
  "pattern": "^\S+\.csv$",
  "errorMessage": "File name must end in '.csv' cannot contain spaces"
}
```

and get:

```
* --input (samples.yml): "samples.yml" does not match regular expression [^\S+\.csv$] (File name must end in '.csv' cannot contain spaces)
```

### `deprecated`

!!! example "Extended key"

A boolean JSON flag that instructs anything using the schema that this parameter/field is deprecated and should not be used. This can be useful to generate messages telling the user that a parameter has changed between versions.

JSON schema states that this is an informative key only, but in `nf-schema` this will cause a validation error if the parameter/field is used.

!!! tip

    Using the [`errorMessage`](#errormessage) keyword can be useful to provide more information about the deprecation and what to use instead.

### `enum`

An array of enumerated values: the parameter must match one of these values exactly to pass validation.

- See the [JSON schema docs](https://json-schema.org/understanding-json-schema/reference/generic.html#enumerated-values)
  for details.
- Available for strings, numbers and integers.

```json
{
  "enum": ["red", "amber", "green"]
}
```

### `fa_icon`

!!! example "Non-standard key"

A text identifier corresponding to an icon from [Font Awesome](https://fontawesome.com/).
Used for easier visual navigation of documentation and pipeline interfaces.

Should be the font-awesome class names, for example:

```json
"fa_icon": "fas fa-file-csv"
```

### `hidden`

!!! example "Non-standard key"

A boolean JSON flag that instructs anything using the schema that this is an unimportant parameter.

Typically used to keep the pipeline docs / UIs uncluttered with common parameters which are not used by the majority of users.
For example, `--plaintext_email` and `--monochrome_logs`.

```json
"hidden": true
```

## String-specific keys

### `pattern`

Regular expression which the string must match in order to pass validation.

- See the [JSON schema docs](https://json-schema.org/understanding-json-schema/reference/string.html#regular-expressions)
  for details.
- Use <https://regex101.com/> for help with writing regular expressions.

For example, this pattern only validates if the supplied string ends in `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz`:

```json
{
  "type": "string",
  "pattern": ".*.f(ast)?q(.gz)?$"
}
```

### `minLength`, `maxLength`

Specify a minimum / maximum string length with `minLength` and `maxLength`.

- See the [JSON schema docs](https://json-schema.org/understanding-json-schema/reference/string.html#length)
  for details.

```json
{
  "type": "string",
  "minLength": 2,
  "maxLength": 3
}
```

### `format`

Formats can be used to give additional validation checks against `string` values for certain properties.

!!! example "Non-standard key (values)"

    The `format` key is a [standard JSON schema key](https://json-schema.org/understanding-json-schema/reference/string.html#format),
    however we primarily use it for validating file / directory path operations with non-standard schema values.

Example usage is as follows:

```json
{
  "type": "string",
  "format": "file-path"
}
```

The available `format` types are below:

`file-path`
: States that the provided value is a file. Does not check its existence, but it does check if the path is not a directory.

`directory-path`
: States that the provided value is a directory. Does not check its existence, but if it exists, it does check that the path is not a file.

`path`
: States that the provided value is a path (file or directory). Does not check its existence.

`file-path-pattern`
: States that the provided value is a glob pattern that will be used to fetch files. Checks that the pattern is valid and that at least one file is found.

### `exists`

When a format is specified for a value, you can provide the key `exists` set to true in order to validate that the provided path exists. Set this to `false` to validate that the path does not exist.

Example usage is as follows:

```json
{
  "type": "string",
  "format": "file-path",
  "exists": true
}
```

!!! note

    If the parameter is an S3, Azure or Google Cloud URI path, this validation is ignored.

!!! warning

    Make sure to only use the `exists` keyword in combination with any file path format. Using `exists` on a normal string will assume that it's a file and will probably fail unexpectedly.

### `mimetype`

MIME type for a file path. Setting this value informs downstream tools about what _kind_ of file is expected.

Should only be set when `format` is `file-path`.

- See a [list of common MIME types](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/MIME_types/Common_types)

```json
{
  "type": "string",
  "format": "file-path",
  "mimetype": "text/csv"
}
```

### `schema`

Path to a JSON schema file used to validate _the supplied file_.

Should only be set when `format` is `file-path`.

!!! tip

    Setting this field is key to working with sample sheet validation and channel generation, as described in the next section of the nf-schema docs.

These schema files are typically stored in the pipeline `assets` directory, but can be anywhere.

```json
{
  "type": "string",
  "format": "file-path",
  "schema": "foo_schema.json"
}
```

!!! note

    If the parameter is set to `null`, `false` or an empty string, this validation is ignored. The file won't be validated.

## Numeric-specific keys

### `minimum`, `maximum`

Specify a minimum / maximum value for an integer or float number length with `minimum` and `maximum`.

- See the [JSON schema docs](https://json-schema.org/understanding-json-schema/reference/numeric.html#range)
  for details.

> If x is the value being validated, the following must hold true:
>
> - `x` ≥ `minimum`
> - `x` ≤ `maximum`

```json
{
  "type": "number",
  "minimum": 0,
  "maximum": 100
}
```

!!! note

    The JSON schema doc also mention `exclusiveMinimum`, `exclusiveMaximum` and `multipleOf` keys.
    Because nf-schema uses stock JSON schema validation libraries, these _should_ work for validating keys.
    However, they are not officially supported within the Nextflow schema ecosystem and so some interfaces may not recognise them.

## Array-specific keys

### `uniqueItems`

All items in the array should be unique.

- See the [JSON schema docs](https://json-schema.org/understanding-json-schema/reference/array#uniqueItems)
  for details.

```json
{
  "type": "array",
  "uniqueItems": true
}
```

### `uniqueEntries`

!!! example "Non-standard key"

The combination of all values in the given keys should be unique. For this key to work you need to make sure the array items are of type `object` and contains the keys in the `uniqueEntries` list.

```json
{
  "type": "array",
  "uniqueEntries": ["foo", "bar"],
  "items": {
    "type": "object",
    "properties": {
      "foo": { "type": "string" },
      "bar": { "type": "string" }
    }
  }
}
```

This schema tells `nf-schema` that the combination of `foo` and `bar` should be unique across all objects in the array.
