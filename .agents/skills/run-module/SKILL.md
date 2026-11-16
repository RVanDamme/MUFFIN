---
name: run-module
description: Run Nextflow Registry modules natively using `nextflow module` commands. Use when running, listing, or getting info about Nextflow modules, including process modules (a single tool) and workflow modules (subworkflows that chain several modules).
---
<!-- Modified for the Nextflow plugin: generic host tools and companion guidance. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Run Nextflow Registry Modules

Use native Nextflow module commands without assuming an input location, output folder, or project layout. Check that the user's Nextflow release supports the module command family.

The bundled command baseline requires Nextflow 26.04 or later. Verify command
availability in the installed release before recommending execution. Workflow
modules require Nextflow 26.10 or later.

### Process modules vs workflow modules

A module has a **kind**, reported as `Kind:` by `nextflow module search` and `nextflow module view`:

| Kind | What it is | Example |
|------|------------|---------|
| `Process` | A single process wrapping one tool | `nf-core/star/align` |
| `Workflow` | A named workflow (subworkflow) that chains several modules; its dependencies are installed with it | `nf-core/fastq_align_star` (STAR + samtools sort/index/stats) |

## Discover and Inspect First

~~~bash
nextflow module search "quality control"
nextflow module view nf-core/fastqc
~~~

Registry identifiers such as nf-core/fastqc name public modules; they are not local filesystem locations.

`nextflow module search` ranks workflow modules below process modules, so they rarely show up in its results. When the task spans several steps (e.g. "align with STAR then sort and index"), also query the registry API with the `kind=Workflow` filter:

~~~bash
curl -s "https://registry.nextflow.io/api/v1/modules?query=align%20reads%20STAR&kind=Workflow&limit=10" \
  | python3 -c 'import json,sys; [print(r["name"], "-", r["description"]) for r in json.load(sys.stdin)["results"]]'
~~~

URL-encode spaces in `query` as `%20`. For a workflow module, `nextflow module view` lists the workflow's `take:` inputs (in call order) and `emit:` outputs.

Read the module's current input contract and run template. Ask for the user's actual input files, parameters, execution environment, and desired result location when they are not established.

## Run the Selected Module

Use nextflow module run with the selected registry identifier and arguments derived from its current template. Supply the user's input references and output choices only after verifying them in the execution environment.

Do not prescribe sample filenames, storage folders, or a wrapper workflow. Modules are fetched on demand when supported by the installed Nextflow release.

If inputs contain several files, follow the template's accepted representation and quoting rules. Resolve any wildcard or file selection from the user's supplied data before running.

### Running Workflow Modules

A workflow module runs with the same `nextflow module run` command, with these differences:

- **It must be typed.** Only workflow modules with typed `take:` inputs can be run directly. Most nf-core workflow modules are untyped and fail with:
  ```
  Workflow `FASTQ_ALIGN_STAR` cannot be executed directly because it is not typed
  ```
  This is not a missing-argument problem. The failed run has already installed the module and its dependencies under `./modules/`. When invoked from registry composition (`build-nextflow-pipeline`) module validation, defer it to the end-to-end run. Otherwise, tell the user the workflow module cannot be run on its own, then offer:
  - running its constituent process modules one by one. They are listed under `requires.modules` in `./modules/<namespace>/<name>/meta.yml`, and any listed workflow module has its own nested `meta.yml` to expand. Or
  - using it inside a pipeline via `build-nextflow-pipeline`.
- **Channel inputs take a samplesheet.** A `Channel<...>` input accepts a CSV, JSON, or YAML samplesheet path; each row becomes one channel item. A `Value<...>` input accepts a single value (e.g. a file path).
- **Outputs are not published.** Each emitted output is reported with its work directory path; nothing is copied to `--outdir`.

## Container Provisioning

Use the module's declared container or dependencies and the user's available runtime. Wave and Conda can provision containers when supported and configured appropriately. A configuration may enable wave.enabled, wave.strategy, and the selected container engine, but do not assume Docker is installed or remote execution uses the local machine.

## Command Reference

| Command | Purpose |
|---|---|
| nextflow module search | Find modules by capability |
| nextflow module view | Inspect a module's inputs and invocation |
| nextflow module run | Execute the selected module with verified arguments |
| nextflow module list | List installed modules with their version and kind |

Use current command help for exact arguments and release support.

## When a Run Fails

1. Inspect the actual error.
2. Reread the selected module's current template and input contract.
3. Compare supplied arguments, file availability, runtime, and parameters.
4. Correct the direct invocation and retry the smallest representative input.

Do not write a wrapper workflow merely to recover from missing or incorrect module arguments.

## Report the Result

Read stdout and summarize status, useful metrics, warnings, and errors. Inspect produced files when the module contract or user's task requires it, using the locations reported by the run.

Suggest a logical next analysis step based on the result. Distinguish an unexecuted recommendation from a completed module run.

## Public Reference

https://registry.nextflow.io
