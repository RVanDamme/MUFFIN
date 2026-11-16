---
name: migrate-nextflow-code
description: Migrate a Nextflow pipeline to new language patterns.
disable-model-invocation: true
---

# Migrate Nextflow Code

Migrate a Nextflow pipeline to new language patterns. Follow best practices from the Nextflow docs. Use `nextflow lint` to validate code as you go. Verify that the pipeline produces the same results before/after migration.

## Requirements

**Requires Nextflow 26.10 or later** for the `nextflow lint` command with type checking.

**Requires nf-core template 3.0.0 or later, if applicable.** If the pipeline has a `.nf-core.yml`, check which version of nf-core/tools last generated its template. If `nf_core_version` is unspecified or less than 3.0.0, stop and tell the user to upgrade their template first. This will resolve syntax errors in the template code and provide a cleaner baseline.

## Migrations

Apply each migration requested by the user, one at a time, in the following order. Each migration requires the ones before it.

Each migration links the relevant docs. Consult these docs when you aren't sure how to do something. For example, which type to use, which method to use for a given type, or how to use a specific channel operator. Don't just guess or try to work from past knowledge. The linter can guide you, but the docs provide examples and best practices.

### Strict syntax

Update a pipeline to comply with `nextflow lint`. Pipelines written before the introduction of the strict syntax parser may contain patterns that are no longer valid.

1. Run the linter over the entire project: `nextflow lint -o concise .`
2. Fix each issue without changing behavior. Keep changes minimal.
3. Re-run the linter as needed and iterate until there are no errors.

Reference:
- [Syntax reference](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/reference/syntax.mdx)
- [Preparing for strict syntax](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/strict-syntax.mdx)

### Topic channels

Replace the `ch_versions` plumbing in a pipeline with a `versions` topic channel.

1. Send process version outputs to the `versions` topic. Replace `emit: versions` with `topic: versions`. Replace the YAML versions file with a tuple of process name, tool name, and eval command.
2. Consume the `versions` topic with `channel.topic('versions')` when writing the versions file. Group elements by process.
3. Delete the `ch_versions` plumbing.

```nextflow
// before
output:
path "versions.yml", emit: versions

// after
output:
tuple val("${task.process}"), val('bash'), eval('bash --version'), topic: versions

// after (typed)
topic:
tuple(task.process, 'bash', eval('bash --version')) >> 'versions'
```

Any process downstream of the topic channel must not emit to that topic, or else the pipeline deadlocks.

Reference:
- [Collecting values with topic channels](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/tutorials/topic-channels.mdx)

### Static typing

Migrate a pipeline to use typed processes, typed workflows, typed parameters, and records. Migrate components from the bottom up. Use `nextflow lint` to type-check as you go.

1. Migrate processes to [typed processes](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/process-typed.mdx)
2. Migrate workflows to [typed workflows](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/workflow-typed.mdx)
3. Migrate the pipeline parameters to [typed parameters](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/typed-parameters.mdx)
4. Re-run the linter as needed until there are no errors.

Rules:
- Use destructured record inputs/outputs in processes.
- Check carefully whether a process `path` input/output should be a single file (`Path`) or a file collection (`List<Path>`, `Set<Path>`, etc).
- Every process should output a single wide record. Every workflow should join related outputs into a wide record channel. This makes it easier to migrate to workflow outputs.
- Use record types to describe channel inputs/outputs for workflows. Define record types with their workflow, not a shared types module.
- If the pipeline uses the `nf-schema` plugin, update it to version 2.7.2 or later to avoid runtime issues with static typing.
- Params used by the script must be declared in the script. Params that only affect config settings should be declared in the config file. Config can override script param defaults, but ideally in a profile (e.g. `test`).
- The global `params` record should only be used in the entry workflow. Pass params to processes and workflows as explicit inputs. Bundle related params into a single workflow input with a record type (e.g. `workflow ALIGNER` -> `record AlignerParams`) so that you can pass `params` directly from the entry workflow.
- Use type coercion (`x as Type`) only when the type is genuinely unknown, such as `channel.topic()`.

Reference:
- [Static typing](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/static-typing.mdx)
- [Standard types](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/reference/stdlib-types.mdx)
- [Migrating to static typing](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/tutorials/static-types.mdx)
- [Using operators with static typing](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/tutorials/static-types-operators.mdx)

### Workflow outputs

Migrate a pipeline to use workflow outputs. Replace all `publishDir` directives with a top-level `output` block that produces the same output directory tree.

If the pipeline hasn't been migrated to static typing and records, stop and warn the user about this. The `output` block requires process outputs to be propagated to the entry workflow, which leads to an excessive number of workflow emits if the pipeline still uses narrow tuples (e.g. `(meta, file)`). Combining narrow tuple channels into wide record channels makes it feasible to publish channels through the `output` block.

1. Identify all `publishDir` directives in processes and config files.
2. For each publisher, determine: which process(es) are targeted, which files are published, and where they are published.
3. For each published output, propagate the corresponding process outputs to the entry workflow (via `emit:`) and publish them (via `publish:`).
4. Add a top-level `output` block with a declaration and `path` directive for each published channel.
5. Define global publishing behavior in config (`outputDir`, `workflow.output.mode`, etc). Override publish settings in the `output` block as needed.
6. Delete the `publishDir` directives.

Rules:
- Consult the docs to learn the `publishDir` options, output `path` directive, and `workflow.output` options. Figure out how to map the former to the latter.
- Prefer one combined per-sample output over many per-tool outputs. Create an index file only if the pipeline already does it (see `collectFile` below).
- You don't need to manually extract or flatten files when publishing a channel. Nextflow automatically extracts files from published values (lists, maps, records, tuples).
- The `collectFile` operator with `storeDir:` is another way to publish a file. Remove `storeDir:` and publish it to the `output` block instead. You might even be able to replace it with an output `index` directive, e.g. if it is a samplesheet of the per-sample outputs.
- Some pipelines have a catch-all `publishDir` in the main config. This can likely be deleted.

Reference:
- [Workflow outputs](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/workflow.mdx#outputs)
- [Migrating to workflow outputs](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/tutorials/workflow-outputs.mdx)
- [Publish configuration](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/reference/config/workflow.mdx)
- [publishDir](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/reference/process/directives/publish-dir.mdx)

### Pipeline composition

Make a pipeline composable, so that it can be executed directly or included by another pipeline.

If the pipeline hasn't been migrated to static typing, or doesn't have the `params` and `output` blocks, stop and warn the user about this. Pipeline composition depends heavily on these language features.

1. Review the `params` block. Determine which params can be refactored as dataflow types. A samplesheet input can be refactored as a channel of records (`Channel<Sample>`). An index file input can be wrapped as a dataflow value (`Value<Path>`). Unless the user says otherwise, refactor samplesheet inputs and leave the rest.

2. Check for global params usage. The pipeline should not use global params outside the entry workflow and `output` block. It should already be publishing all outputs through the `output` block instead of `publishDir` directives.

3. Check for project-level assets (`projectDir`, `bin`, `lib`). The pipeline should avoid these in favor of composable alternatives (params, module `resources/` and `moduleDir`, helper functions).

4. Check for process configuration (`conf/modules`, `ext.args`, `withName`). Refactor it by moving it into pipeline code (see below).

The biggest issue you will likely encounter is `ext` config. nf-core modules use `ext` for certain inputs such as tool args and filename prefix/suffix. Users can then override these settings easily from config. A pipeline may specify its own defaults in config with process selectors. This makes the pipeline harder to compose.

The solution is to refactor the pipeline's `ext` config as pipeline code:
- Modules take `args` / `prefix` / `suffix` as explicit inputs and resolve them as `task.ext.<name> ?: <input> ?: <default>`. Users can still override any tool at runtime with config.
- Each `ext` setting moves into pipeline code. Add the setting to each sample record with `map` or `combine` depending on whether it varies per-sample.
- Use helper functions to separate complex `ext` settings from workflow logic. Build an options map first, then render it to a CLI string with a separate `cli()` helper.

Reference:
- [Pipeline composition](https://raw.githubusercontent.com/nextflow-io/nextflow/master/docs/workflow-typed.mdx#pipeline-composition)
- [fetchngs -> rnaseq example](https://github.com/nextflow-io/nextflow/tree/master/examples/pipeline-composition)
- [Modules](https://raw.githubusercontent.com/nextflow-io/nextflow/refs/heads/master/docs/modules/modules.mdx)

## Validation

After each migration, verify behavior against baseline using the project's test suite. For example:

- Run before/after with test profile: `nextflow run . -profile test,docker`
- nf-test: `nf-test test`
- Stub run: `nextflow run . -stub`

## Critical Rules

1. **One migration at a time.** Finish one migration before starting another one.
2. **Read the docs.** Consult the linked docs when you aren't sure how to do something.
3. **Preserve behavior.** The pipeline should produce the same outputs before/after the migration.
