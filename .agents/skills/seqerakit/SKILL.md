---
name: seqerakit
description: >
  Write seqerakit YAML configuration files for automating Nextflow Platform setup.
  Covers pipelines, launch, compute-envs, datasets, credentials, and other entity types.
  Trigger: "seqerakit", "seqerakit YAML", "write seqerakit config", "automate platform setup",
  "tw YAML", "infrastructure as code for Seqera".
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


# Seqerakit YAML Configuration

Draft seqerakit configuration using the user's intended Seqera workspace, resources, inputs, and output choices. Do not assume a local configuration location or prescribe storage folders.

## Critical Rules

Confirm fields using the installed tw CLI's entity-specific help or current seqerakit documentation. Unknown fields can be silently ignored, so do not guess field names.

Use hyphenated field names, including compute-env, work-dir, params-file, and file-path. Use revision for a Git branch, tag, or commit. Verify the workspace identifier from the user's existing context.

The entity-field lists below describe the bundled reference baseline and may differ from the installed release.

## Entity Fields

### Pipelines

Use name, url, workspace, description, compute-env, revision, profile, params-file, work-dir, pre-run, and labels only when supported by current help.

The pipeline URL must identify the requested public or private pipeline. File and storage fields must come from user-provided resources or verified Platform context.

### Launch

Use the selected pipeline, workspace, compute environment, revision, profiles, and parameters. Supply params-file, config, pre-run, or work-dir only when the user has provided or selected those resources.

### Compute Environments

Common fields include name, workspace, credentials, type, file-path, and wait. Seqera Compute may not require user-supplied cloud credentials. Treat provider configuration as specific to the user's cloud environment.

### Datasets

Common fields include name, workspace, description, header, and file-path. Establish the user's actual dataset and format instead of inventing a dataset filename or folder.

### Credentials

Common fields include name, workspace, and type, with provider-specific fields. Reference existing credential identifiers whenever possible; do not expose secrets in generated examples or responses.

## Generic Configuration Example

~~~yaml
pipelines:
  - name: 'rnaseq'
    url: 'https://github.com/nf-core/rnaseq'
    revision: 'master'
    on_exists: 'fail'

launch:
  - name: 'rnaseq-run'
    pipeline: 'https://github.com/nf-core/rnaseq'
    params: {}
~~~

This illustrates the entity structure, not a ready-to-run launch. Add the workspace, compute environment, required parameters, and any file references from the user's verified context.

## Existing Resources

Check the installed release's on_exists behavior before choosing fail, ignore, or overwrite. Do not overwrite existing resources merely to make a configuration succeed.

The older overwrite setting is deprecated in the bundled baseline; use the
installed release's supported on_exists behavior instead.

Use only the requested changes, validate the YAML, and explain any missing values. Running the configuration requires the user's intended scope and an authenticated CLI available in the execution environment.
