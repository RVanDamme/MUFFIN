---
name: nf-pipeline-design
description: >
  Make sure to use this skill whenever you are designing, writing, reviewing,
  or refactoring a Nextflow DSL2 pipeline with a focus on STRUCTURE — the
  shape of `main.nf`, where subworkflows belong, when something is a module
  vs. an operator, how channels and metadata flow. Trigger phrases include
  "design a pipeline", "structure my Nextflow code", "review main.nf",
  "refactor this workflow", "is this a module or a subworkflow", "where does
  this logic belong", "clean up AI-generated Nextflow", "make this more
  cloud-efficient", "module boundary", "subworkflow boundary", "channel
  shape", "tuple shape", "operator vs module". Also use it to analyze how an
  existing pipeline is organized — processes, modules, subworkflows, channels
  and data flow — when the user asks how a pipeline works or before changing
  it. Use this skill before writing any new `.nf` file in an unfamiliar layout.
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


# Rules for a clean Nextflow pipeline

## Companion skills

When a structural decision depends on a non-structural one, chain to:

- `find-alternative-tools` — tool choice at a branch
- `search-existing-modules` — reuse vs author a module
- [parameter triage](../nextflow-schema/references/triage-pipeline-parameters/README.md) — param schema before `main.nf`

## Canonical Nextflow baseline

For new pipelines, target **Nextflow 26.04+** unless the user explicitly asks for an older runtime. Treat the current Seqera / Nextflow docs as source of truth when this skill, older nf-core examples, or model memory disagree.

Required defaults for code generation from scratch:

- write code that passes the **strict syntax parser**; it is the default in Nextflow 26.04+
- set `nextflow.enable.types = true` in every `.nf` script that defines typed processes or typed workflows
- keep typed and legacy process/workflow definitions in separate scripts; `nextflow.enable.types` does not allow mixing them in one script
- use typed `params {}` blocks for pipeline parameters when writing typed scripts
- use lowercase `channel.*` factories
- use typed `input:` / `output:` / `topic:` sections in typed processes
- avoid `.out` access in typed workflows; assign call returns to variables instead
- do not add deprecated `nextflow.enable.dsl` or `nextflow.enable.strict` flags

> Nextflow 26.04 enables the strict syntax parser by default. Do **not** instruct agents to set `NXF_SYNTAX_PARSER=v2` for 26.04+ runs — it is redundant. For transitional 25.04 / 25.10 validation, set `NXF_SYNTAX_PARSER=v2` explicitly so `nextflow lint` and test runs exercise the strict parser.

Canonical sources to check when updating this section:

- Strict parser defaults and migration rules: https://docs.seqera.io/nextflow/strict-syntax
- Feature flags, especially `nextflow.enable.types`: https://docs.seqera.io/nextflow/reference/feature-flags
- Workflow outputs and typed workflow call behavior: https://docs.seqera.io/nextflow/workflow
- Modules and registry modules: https://docs.seqera.io/nextflow/module

## Verification — always run checks before finishing

At minimum, run **`nextflow lint` against the user-selected pipeline sources** with the target Nextflow version. It catches parse-level and strict-syntax issues: deprecated DSL1 forms, disallowed Groovy constructs, `Channel.*` vs `channel.*`, `import` statements, `${PWD}` implicit env refs, and other strict-syntax violations. Use the lint formatter only after reviewing what it will change.

If available, use the bundled `nflint` structural linter with Python 3.10+ to
check architectural rules (NF001–NF052) against the user's selected pipeline.
Select relevant rules or request its rule list through the available execution tools.

If the repository provides additional validators, run them too. Examples include nf-core lint and nf-test. A clean parse/lint run does not prove the design is good — semantic items such as umbrella subworkflow names, vague variable names, hardcoded `meta.*` fields, module heft, hardcoded optional args, complex arg strings in `main.nf`, and unnecessary Python wrappers still require human review against the rules below.

## Pipeline layout

Separate entrypoint orchestration, reusable workflow composition, and individual
compute processes. Keep configuration, helper commands, tests, and supporting
data distinct where that helps the user's project, without imposing a layout. For file structure expectations, read https://nf-co.re/docs/contributing/pipelines/pipeline_file_structure when you need to check the canonical layout.

## Reference-pipeline patterns to copy

For canonical reference, mine **strict-syntax-clean nf-core pipelines** — pinned versions that pass strict-parser linting of the selected sources with zero errors. The [`nf-core/strict-syntax-health`](https://github.com/nf-core/strict-syntax-health) dashboard tracks which release tags (and dev-branch commits) are currently clean. Representative pipelines: `nf-core/rnaseq`, `nf-core/ampliseq`, `nf-core/demo`, `nf-core/raredisease`, `nf-core/pixelator`. Browse the full catalog at https://nf-co.re/pipelines. Use these pipelines as examples of structure, not as a license to copy legacy syntax blindly. The strongest recurring patterns are:

- `main.nf` stays thin: it imports the main workflow and initialization/completion helpers, and then calls the pipeline workflow from a short entry workflow.
- Initialization owns parameter validation, run summaries, samplesheet parsing, reference setup, and profile sanity checks. Completion owns MultiQC/report assembly, versions collation, and final run metadata.
- Samplesheets are converted to typed channel shapes at the boundary, commonly via `samplesheetToList` from `plugin/nf-schema` or pipeline utility subworkflows. Do not let ad-hoc CSV parsing leak through modules.
- Cross-cutting streams use topic channels where available: `channel.topic('versions')` for software versions and `channel.topic('multiqc_files')` for report fragments. Avoid plumbing seeded `ch_versions = channel.empty()` through every subworkflow in new code.
- Newer reference pipelines use top-level `output {}` blocks for workflow outputs. Existing nf-core pipelines may still use config-driven `publishDir`; follow the existing pipeline convention when contributing, but prefer workflow outputs for new standalone strict-syntax examples.
- Utility plugins such as `nf-schema` and `nf-core-utils` should be used deliberately: import only the helpers you call, and call validation helpers near the entrypoint instead of leaving plugins decorative.

When mining strict-syntax-clean nf-core pipelines, copy the shape and separation of concerns: entrypoint → workflow → subworkflow → module. Do not copy stale config flags, deprecated strict-mode flags, or older topic/version patterns if current docs say otherwise.

## `main.nf`: keep the entrypoint thin

`main.nf` is the entrypoint of the pipeline. It should be readable top-to-bottom and answer these questions immediately:

- what are the pipeline inputs
- what are all analysis stages
- what are the high-level outputs

`main.nf` should usually outline the pipeline structure:

- initialize and validate parameters
- create top-level channels and values
- call named workflows
- perform only shallow orchestration logic

`main.nf` should usually avoid:

- deep chains of operators
- large blocks of general-purpose Groovy
- tool-specific argument construction
- process-level branching hidden in `when:`
- direct use of `params.*` deep inside the pipeline

The entry workflow can contain channels and operators, but it should stay shallow enough that reading `main.nf` is enough to understand the pipeline story.

When touching workflow structure, `take:`, `emit:`, workflow outputs, or entry workflow behavior, read https://docs.seqera.io/nextflow/workflow.

### `main.nf` structure

`main.nf` should be organized in two conceptual subsections:

- a setup section
- a run section

The setup section should expose all parameters that are actually used by the pipeline and assign them to explicit channels or values with descriptive, jargon-light, pipeline-adapted names.

For example, in a protein design pipeline, `fasta` is often too vague, while `protein_sequence_to_be_folded` is much more informative.

The run section should use names and comments that explain the flow at a glance.

Bad:

```nextflow
MODEL_INFERENCE(
    ch_fold_ready_sequences,
    ch_template_hits,
    folding_model_name
)
```

Better:

```nextflow
/*
Run folding inference.
Requires the prepared inputs, template hits, and model selection.
This stage performs final structure prediction.
*/
FOLDING_MODEL_INFERENCE(
    ch_fold_ready_sequences,
    ch_template_hits,
    folding_model_name
)
```

### Good `main.nf` shape (classic DSL2)

```nextflow

workflow {
    // ----------------------------
    // Parameter setup
    // ----------------------------
    ch_input_sequences = channel.fromPath(params.input, checkIfExists: true)

    // ----------------------------
    // Pipeline run
    // ----------------------------
    /*
    Normalize and validate the raw sequence inputs so downstream folding steps
    consume one predictable channel shape.
    */
    PREPARE_INPUTS(ch_input_sequences)

    /*
    Run the folding stage on prepared sequences using the selected folding model.
    This stage consumes fold-ready inputs and produces structure predictions.
    */
    RUN_FOLDING(
        PREPARE_INPUTS.out.fold_inputs,
        params.model_name
    )

    /*
    Gather the prediction outputs into the final reporting structure and publish
    the pipeline-facing result artifacts.
    */
    COLLECT_REPORTS(RUN_FOLDING.out.predictions)
}
```

### Good `main.nf` shape (typed processes/workflows — Nextflow 26.04+)

In typed workflows, `.out` is not available. Assign the return value of each
workflow/process call to a variable and access named outputs on it. In legacy workflow bodies calling typed processes, `.out` may still work, but new typed examples should avoid it.

**Important:** When a subworkflow has only **one** `emit:`, the return value IS
that channel directly — do not use `.name` on it (it returns `null`). Use `.name`
only for multi-emit subworkflows. See "Typed output access" below for the full rules.

```nextflow
nextflow.enable.types = true


workflow {
    // ----------------------------
    // Parameter setup
    // ----------------------------
    ch_input_sequences = channel.fromPath(params.input, checkIfExists: true)

    // ----------------------------
    // Pipeline run
    // ----------------------------
    /*
    Normalize and validate the raw sequence inputs so downstream folding steps
    consume one predictable channel shape.
    PREPARE_INPUTS has a single emit, so return value is the channel directly.
    */
    ch_fold_ready = PREPARE_INPUTS(ch_input_sequences)

    /*
    Run the folding stage on prepared sequences using the selected folding model.
    RUN_FOLDING has multiple emits (predictions, metrics), so use named access.
    */
    folding = RUN_FOLDING(ch_fold_ready, params.model_name)

    /*
    Gather the prediction outputs into the final reporting structure and publish
    the pipeline-facing result artifacts.
    */
    COLLECT_REPORTS(folding.predictions)
}
```

### Rules for `main.nf`

- Use explicit names that reflect pipeline meaning, not raw file format names, for example `protein_sequence_to_be_folded` instead of `fasta`, `ch_template_hits` instead of `hits`, and `folding_model_name` instead of `model`.
- Use comments that explain the analysis stage, not comments that narrate syntax.
- Pass explicit inputs to workflows and processes instead of reading `params.*` everywhere.
- If the repository supports strict syntax, use the modern `params` block and type annotations.
- Keep pipeline-specific validation near the entrypoint or schema, not scattered through modules.
- `main.nf` should ideally use only subworkflow calls.
- `main.nf` should avoid channel manipulation logic because it obscures the top-level story.
- `main.nf` should avoid generic Groovy helper code for the same reason.
- Reading `main.nf` should usually be enough to understand what the pipeline does without descending into implementation details.

For strict-syntax-era workflow best practices, especially around `params` usage and where to place conditional logic, read https://docs.seqera.io/nextflow/strict-syntax.

## Subworkflows: own the glue logic

Subworkflows are where most pipeline-specific dataflow logic should live.

That includes:

- channel reshaping
- joining metadata and files
- grouping or branching data
- choosing between alternative execution paths
- adapting one module interface to another

This is the right level for heavy use of native operators.

Before writing operator-heavy logic, read the actual Nextflow operator documentation online rather than relying on memory:

- typed operators: `https://docs.seqera.io/nextflow/reference/operator-typed`
- legacy operators: `https://docs.seqera.io/nextflow/reference/operator`

Use the documented operators that actually match the problem, such as `map`, `filter`, `join`, `mix`, `combine`, `branch`, `collect`, `flatMap`, `reduce`, and the appropriate grouping operators.

### Subworkflow shape

A good subworkflow is often built in three parts:

- channel manipulation logic
- module execution
- optional output reshaping

Use small, specific subworkflows over vague umbrella subworkflows. A subworkflow should execute a single module and use a name that mirrors the operation it represents. If multiple modules are interchangeable, select between them with workflow-level `if/else`, but only one module should run in a given execution path.

For example, this is usually less informative:

```nextflow
FOLDING_PIPELINE(...)
```

And this is more informative:

```nextflow
GET_TEMPLATES(...)
RUN_MSA(...)
FOLDING_MODEL_INFERENCE(...)
```

This naming preference exists mainly so that `main.nf` remains informative.

### Branch in workflows, not in processes

Put conditional logic in the calling workflow with `if`, `filter`, or `branch` rather than in process `when:` blocks.

Read https://docs.seqera.io/nextflow/strict-syntax before implementing branching rules.

```nextflow
workflow RUN_FOLDING {
    take:
    ch_fold_inputs
    model_name

    main:
    if( model_name == 'boltz' ) {
        BOLTZ_PREDICT(ch_fold_inputs)
        ch_predictions = BOLTZ_PREDICT.out.predictions
    }
    else if( model_name == 'openfold3' ) {
        OPENFOLD3_PREDICT(ch_fold_inputs)
        ch_predictions = OPENFOLD3_PREDICT.out.predictions
    }
    else {
        error "Unsupported model_name: ${model_name}"
    }

    emit:
    predictions = ch_predictions
}
```

Version aggregation is intentionally not shown here. New code should emit tool
versions to the `versions` topic from each module and subscribe once in the entry
workflow.

### Rules for subworkflows

- Give subworkflows clear inputs and outputs with `take:` and `emit:`.
- Emit named outputs and consume them by name, not by numeric position.
- Put non-trivial channel logic here rather than bloating `main.nf`.
- Use native operators, and when needed light workflow-level Groovy, for metadata reshaping, grouping, merging, and filtering.
- Use subworkflows to express analysis steps, not vague buckets like `RUN_PIPELINE_STAGE_2`.
- Subworkflows should own the branching logic.
- Select tools with workflow-level `if/else` on an explicit input such as `model_name`, instead of burying that choice in scripts or process directives.

### Avoid `channel.empty()` except for genuinely optional branches

Do not teach agents to seed routine dataflow with `channel.empty()`. Prefer explicit branches, typed workflow returns, and topic channels. In particular, do **not** thread routine version aggregation as `ch_versions = channel.empty(); ch_versions = ch_versions.mix(...)` in new code — use the versions topic channel below.

Use `channel.empty()` only when a branch is genuinely optional and downstream consumers are designed to handle zero emissions:

```nextflow
workflow OPTIONAL_QC {
    take:
    ch_reads
    run_qc

    main:
    if (run_qc) {
        FASTQC(ch_reads)
        ch_reports = FASTQC.out.reports
    } else {
        ch_reports = channel.empty()
    }

    emit:
    reports = ch_reports
}
```

Key points:
- Do not use `channel.empty()` as a silent fallback to hide missing data. If a process should always produce output, let it fail.
- If you must use `channel.empty()`, document why zero emissions are valid for that branch.
- Prefer topic channels for cross-cutting streams such as versions and MultiQC fragments.

### Cloud execution cost

When running on the cloud, extra modules are expensive.

Each module execution becomes a scheduled task:

- it has to be submitted to the batch queue
- it may wait for resources before it starts
- it adds task overhead even if the job is only doing simple glue logic

By contrast, channel operations and light workflow-level Groovy run in the Nextflow runtime on the head job and do not require another batch submission.

Because of that:

- use channel operators instead of creating modules for simple metadata reshaping
- use workflow-level `if/else` and branching instead of wrapper modules whose only job is routing
- use light Groovy in subworkflows when it keeps the dataflow clear and avoids spawning a pointless task

Do not create a module just to rename fields, regroup tuples, merge metadata, select one branch, or perform other cheap glue operations that Nextflow can handle directly.

### Versions topic channel

Use a `versions` topic channel to aggregate tool version information across every module. This is the Nextflow 25.04+ / nf-core 2026 idiom — it replaces the older pattern of threading a `ch_versions` channel through every subworkflow with `mix()`.

#### Classic DSL2 topic syntax

Each module publishes a **tuple** containing process name, tool name, and version to the `versions` topic directly in its `output:` block. Use the `eval` output qualifier to capture the version string at runtime — no `versions.yml` file or heredoc needed:

```nextflow
process SAMTOOLS_SORT {
    // ...

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -1 | sed 's/samtools //'"), topic: versions

    // ...
}
```

For tools where the version is known statically (e.g. a Python script bundled in the container), use `val` instead of `eval`:

```nextflow
output:
tuple val("${task.process}"), val('custom_tool'), val("1.2.0"), topic: versions
```

The entry workflow subscribes to the topic once and formats the tuples into a YAML report:

```nextflow
workflow {
    // ... module / subworkflow calls ...

    channel.topic('versions')
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            def dedup = tool_versions.unique().sort()
            "${process}:\n${dedup.join('\n')}"
        }
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
            newLine: true,
        )
}
```

#### Typed process topic syntax (Nextflow 26.04+)

With typed processes (`nextflow.enable.types = true`), the topic is a **separate `topic:` section** — not part of `output:`. Tuples are emitted to a named topic with the `>>` operator:

```nextflow
nextflow.enable.types = true

process SAMTOOLS_SORT {
    // ...

    output:
    sorted_bam = tuple(meta, file("*.bam"))

    topic:
    tuple('samtools', eval('samtools --version | head -1')) >> 'versions'

    // ...
}
```

Key differences from classic DSL2:
- `topic:` is its own section, separate from `output:`
- Tuples use the `>> 'topic_name'` operator to specify the target topic
- The process name is **not** included in the tuple — the runtime attaches it automatically
- Use `eval(...)` for runtime version capture

The entry workflow subscription is simpler because each tuple has only `(tool, version)`:

```nextflow
channel.topic('versions')
    .unique()
    .map { tool, version ->
        "${tool}: ${version}"
    }
    .collectFile(
        name: 'software_versions.yml',
        storeDir: "${params.outdir}/pipeline_info",
        newLine: true,
        sort: true,
    )
```

#### Common rules for both styles

A subworkflow no longer needs to declare a `versions` output or thread `ch_versions = ch_versions.mix(MOD.out.versions)` through every step. The topic collects emissions implicitly across the whole pipeline.

> **Deadlock rule:** any process that consumes a topic must not also emit to it. The canonical versions pattern is safe because the topic is consumed only by `collectFile()` in the entry workflow, not by any module.

> **Backwards compatibility:** nf-core supports both the old `path "versions.yml", emit: versions` file-based pattern and the new tuple+topic pattern simultaneously during the transition. Pipelines can mix both styles — the entry workflow can branch on `instanceof Path` to handle legacy modules that still emit files. For **new code**, always use the typed process `topic:` section above.

For subworkflow conventions, version aggregation, and minimum subworkflow size, read https://nf-co.re/docs/guidelines/components/subworkflows.

## Modules: keep them atomic and predictable

Modules are the unit of heavy computation and the main reuse boundary.

Before writing or editing modules, read:

- https://docs.seqera.io/nextflow/module
- https://nf-co.re/docs/guidelines/components/modules

Good modules are:

- small
- explicit
- easy to test
- easy to containerize
- boring to read

As a default, one module should wrap one tool invocation or one tightly coupled subcommand. Do not write giant "do everything" modules unless there is a real runtime or I/O reason.

### Module registry and single-module validation

Nextflow 26.04 introduced registry modules managed by the native `nextflow module` command. For nf-core-backed work, prefer this workflow before composing a pipeline:

1. Search by capability, not by guessed name:
   ```bash
   nextflow module search "FASTQ quality control"
   ```
2. Inspect the selected module's contract and run template:
   ```bash
   nextflow module view nf-core/fastqc
   ```
3. Run the module directly with `nextflow module run`, using the selected
   module identifier, a user-provided test input, and the chosen output destination.
   Follow the arguments shown by `nextflow module view`.
4. Only after every candidate module has been validated should you compose them into subworkflows and the final pipeline.

Never write a wrapper workflow just to test one registry module or to recover from missing arguments. If `nextflow module run` fails, go back to `nextflow module view`, compare your command to the generated template, fix the arguments, and rerun the module directly. Wrapper workflows are for composing multiple validated modules, not for discovering a single module's CLI contract.

When integrating a registry module into pipeline code, use registry include syntax unless you intentionally vendor or patch the module locally:

```nextflow
include { FASTQC } from 'nf-core/fastqc'
```

For project-owned or intentionally modified vendored modules, resolve the include
from the source files supplied by the user rather than imposing a component location.


### Fan-out modules: raw outputs in modules, shaping in subworkflows

When a module splits one input into many outputs, keep the module reusable by
emitting raw files. Attach caller-specific metadata in the subworkflow that owns
the pipeline context.

```nextflow
// main.nf — module emits raw records
process SPLIT_FASTA {
    input:
    path fasta: Path

    output:
    records: Set<Path> = files("*.fa")

    script:
    """
    split_fasta.py --input ${fasta}
    """
}
```

```nextflow
// split_input_fasta.nf — subworkflow shapes metadata
workflow SPLIT_INPUT_FASTA {
    take:
    ch_fasta

    main:
    ch_records = SPLIT_FASTA(ch_fasta)
        .flatten()
        .map { record -> tuple([id: record.baseName], record) }

    emit:
    records = ch_records
}
```

Agents often try to attach parent metadata inside the module. That couples the
module to one caller and makes fan-out awkward. Let the module do the heavy tool
work; let the subworkflow shape output tuples for downstream consumers.

### Rules for modules

- Required files belong in `input:` definitions.
- Required non-file arguments should come through explicit value inputs.
- Optional tool arguments should usually come from `task.ext.args` configured in `modules.config`, not hardcoded in `main.nf`.
- Keep module naming and channel naming consistent and predictable.
- Use tuples and meta maps to move structured state through the pipeline instead of dumping JSON sidecars for everything.
- Do not hardcode custom `meta.*` fields inside reusable modules.
- Modules are the place for heavy computation, not for lightweight data reshaping.

### `stageAs` to avoid filename collisions

When a process receives multiple files that may share the same name (e.g., an index and a secondary index from different upstream tasks), Nextflow will error on staging. Use the `stageAs` option (or wildcard patterns) to rename files on the way into the task directory:

Configure distinct staging patterns for the two input sets, then pass the staged
file references into the merge command. Preserve original names or introduce
stable prefixes as needed to avoid collisions.

Choose staging patterns from the user's actual input cardinality and collision
requirements. Nextflow supports original names, indexed filename patterns, and
separate staged groups. For typed processes, express the same decision in the
`stage:` section. Consult the installed version's `stageAs` documentation before
changing these patterns.

### Fail loudly instead of emitting silent emptiness

If a module is supposed to produce an output, let it fail loudly when that output is missing.

Do not hide failures by:

- returning placeholder files
- swallowing exit codes
- emitting generic empty channels as a substitute for proper logic

Use optional outputs only when absence is a valid, expected state of the computation.

### Tag every process for useful logs

The `tag` directive labels each task in the execution log and the Nextflow Platform run view. Without it, a failed task shows up as `FORMAT_FASTA (1)` and it's unclear which sample broke. With it, you get `FORMAT_FASTA (sample_42)` — immediately actionable.

```nextflow
tag "${meta.id}"
```

Always add a `tag` to any process that operates per-sample or per-record. Use the most human-readable identifier available — `meta.id`, a file base name, or a combination.

### `errorStrategy`: prefer config over modules

Keep `errorStrategy` in config, not in module files. A module-level directive can still be overridden from config with `process { withName: 'FOO' { errorStrategy = ... } }`, but spreading retry policy across many modules makes per-site tuning harder to reason about — there's no single place to read or change it.

A common nf-core convention is to set one global policy in the shared process configuration:

```groovy
// Shared process configuration
process {
    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'terminate' }
    maxRetries    = 1
}
```

These are Linux signal-derived exit codes (128 + signal number), not cloud-specific: 130 (SIGINT), 137 (SIGKILL — typically OOM-killed), 143 (SIGTERM), and so on. Exit code 104 is ECONNRESET, common on preempted spot / interruptible instances. The retry covers transient infrastructure failures regardless of executor.

**Caveat:** many bioinformatics tools do not respect Linux exit codes — they return `1` for every failure mode, transient or not. Retrying on `1` would mask deterministic tool failures (bad input, missing reference, malformed parameters), so this pattern intentionally only retries on signal-derived codes. If a tool you rely on uses non-standard exits for transient errors, handle it with a targeted `withName:` override rather than widening the global rule.

Only declare `errorStrategy` in a specific module when that process has a well-understood, non-transient failure mode that is genuinely different from the global policy (e.g. a download process that should `ignore` 404s). When you do, document why in a comment beside the directive.

### `stub:` blocks for dry-run testing

Modules should include a `stub:` block alongside `script:`. The stub creates empty output files so the pipeline can be exercised end-to-end without running the real tool:

```nextflow
stub:
def prefix = task.ext.prefix ?: meta.id
"""
touch ${prefix}.fa
touch versions.yml
"""
```

**Never use `-stub-run` as your verification step.** A passing stub run does not mean the pipeline works — it only checks that outputs are declared, not that anything actually works. When verifying your work, always run the pipeline with real (possibly minimal) data.

AI assistants reach for `-stub-run` to get a green exit code without solving the actual problem: code doesn't work? Stubs. Docker doesn't run? Stubs. Workflow doesn't make sense? Stubs. This is not debugging — it is hiding the bug. If the real run is failing, the correct response is to diagnose the script, container, or input problem, not to switch to a mode that bypasses all of them.

Invest effort in slicing test data down to a minimal synthetic set that runs the full pipeline end-to-end in under 5 minutes. That is harder than adding stubs, but it actually proves the pipeline works.

### Example module pattern (classic DSL2)

```nextflow
process FORMAT_FASTA {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/format-fasta:1.0.0--py310h1234567_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val("${task.process}"), val('format_fasta'), eval("format_fasta.py --version"), topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    format_fasta.py \
        ${args} \
        --input ${fasta} \
        --output ${prefix}.fa
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.fa
    """
}
```

### Example module pattern (typed process — Nextflow 26.04+)

With typed processes, `input:` uses `name: Type` syntax, `output:` uses named assignments,
`topic:` is a separate section, and `stage:` replaces `stageAs:`. Enable with
`nextflow.enable.types = true` in each `.nf` script that uses typed processes/workflows. Strict syntax is the default in 26.04+.

```nextflow
nextflow.enable.types = true

process FORMAT_FASTA {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/format-fasta:1.0.0--py310h1234567_0'

    input:
    tuple(meta: Map, fasta: Path)

    output:
    formatted = tuple(meta, file("*.fa"))

    topic:
    tuple('format_fasta', eval('format_fasta.py --version')) >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    format_fasta.py \
        ${args} \
        --input ${fasta} \
        --output ${prefix}.fa
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.fa
    """
}
```

Key typed-process rules:

- **Input types**: `name: Type` — e.g. `fasta: Path`, `kmer_size: Integer`, `meta: Map`. Tuples: `tuple(meta: Map, fasta: Path)`.
- **Output assignments**: `name = expression` — e.g. `result: Path = file("*.txt")`, `bam = tuple(meta, file("*.bam"))`. The name becomes the output channel accessor.
- **`topic:` section**: separate from `output:`. Emit with `tuple(...) >> 'topic_name'`.
- **`stage:` section**: replaces `stageAs:` qualifier — e.g. `stageAs files, user_staging_pattern`.
- **Collections**: use `List<Path>` for collected outputs (from `.collect()`), `Set<Path>` only for `files("glob")` which returns a set at the process level.
- **No `.out` property**: assign the process call return to a variable and access named outputs on it: `result = MY_PROCESS(ch_input); result.bam`.
- **`nextflow.enable.types = true`**: required in each `.nf` script that uses typed processes/workflows. Strict syntax is the default in 26.04+. If you must support older transitional builds, verify the exact feature flag expected by that installed Nextflow version before running lint/tests.

### Typed output access: when to use `.name` vs direct pass

The return value of a process or subworkflow call behaves differently depending on how many outputs/emits it declares. Getting this wrong is one of the most common typed-process bugs — the code looks correct but channels are empty or mistyped at runtime.

**Rules:**

| Call type | Outputs declared | Access pattern |
|---|---|---|
| Process with **multiple** named outputs | 2+ | `result = MY_PROCESS(ch); result.bam; result.stats` |
| Process with **one** named output | 1 | `result = MY_PROCESS(ch)` — the return IS that output channel directly. `.name` also works but is redundant. |
| Subworkflow with **multiple** emits | 2+ | `result = MY_SUBWORKFLOW(ch); result.aligned; result.versions` |
| Subworkflow with **one** emit | 1 | `result = MY_SUBWORKFLOW(ch)` — the return IS the emitted channel. Do NOT use `.name` — it returns `null`. |

**Example — single-emit subworkflow (common case):**

```nextflow
// validate_sequences.nf
workflow VALIDATE_SEQUENCES {
    take:
    ch_records

    main:
    result = VALIDATE_SEQUENCE(ch_records)

    emit:
    sequences = result.validated   // process has named output 'validated'
}

// main.nf — WRONG: subworkflow has one emit, .sequences returns null
ch_validated = VALIDATE_SEQUENCES(ch_records).sequences  // ❌

// main.nf — CORRECT: single-emit subworkflow returns the channel directly
ch_validated = VALIDATE_SEQUENCES(ch_records)            // ✅
```

**Example — multi-emit subworkflow:**

```nextflow
// align_reads.nf
workflow ALIGN_READS {
    take:
    ch_reads

    main:
    result = BWA_MEM(ch_reads)

    emit:
    aligned  = result.bam
    unmapped = result.unmapped
}

// main.nf — multi-emit requires named access
align_result = ALIGN_READS(ch_reads)
ch_bam       = align_result.aligned     // ✅ — must use .name
ch_unmapped  = align_result.unmapped    // ✅
```

**Collection types — `List<Path>` vs `Set<Path>`:**

The `.collect()` operator produces an `ArrayBag` which is `List`-compatible, NOT `Set`-compatible. Use the correct type annotation in process inputs:

```nextflow
// CORRECT — .collect() feeds this input
input:
all_stats: List<Path>

// WRONG — .collect() does not produce a Set
input:
all_stats: Set<Path>
```

Reserve `Set<Path>` only for process outputs using `files("glob")`, which genuinely returns an unordered set of matched files.

**Static topic versions:**

For tools where the version is known at pipeline-authoring time (e.g., a pinned container), use a static string in the `topic:` section instead of `eval(...)`:

```nextflow
topic:
tuple('awk', '5.1') >> 'versions'
```

Use `eval(...)` only when the version is genuinely dynamic (e.g., captured from a tool's `--version` flag inside a container that may be updated independently).

## Helper commands and containers

Many bad Nextflow pipelines come from treating the workflow like a Python application runner instead of a dataflow system.

### Use scripts correctly

If a script is needed, provide it as a project-owned command or a module-owned
binary using the packaging mechanism supported by the user's execution environment.

Then:

- give it a shebang
- make it executable
- call it directly by name

Do this:

```bash
format_fasta.py --input input.fa --output output.fa
```

Avoid hardcoded host locations or dependence on a runtime absent from the task's
execution environment.

### Cloud nuance

- Ensure project helper commands are actually packaged for the selected executor;
  local availability does not prove they will be available on remote compute.
- Module-owned binaries require `nextflow.enable.moduleBinaries = true`, and on cloud executors they require Wave.

This matters for cache behavior, portability, and container fingerprints.

When working with project helper commands, module binaries, or cloud packaging behavior, read https://docs.seqera.io/wave/nextflow/bundle-scripts.

### Container rule

Every module process must declare an explicit container image. The unit of reproducibility is the image — a hashed, cacheable artefact — not a run-time package solve.

```nextflow
// split_fasta.nf
process SPLIT_FASTA {
    container 'quay.io/nf-core/ubuntu:20.04'
    ...
}
```

```nextflow
// nextflow.config — enable Docker (or Singularity, Podman, etc.) plus Wave.
docker { enabled = true }
wave   { enabled = true }
```

Wave (https://docs.seqera.io/wave) is the recommended path for new pipelines. It mirrors and caches whatever image the `container` directive names, handles private-registry auth, and gives every task a stable, content-addressed fingerprint. On HPC, swap `docker` for `singularity` in the config without editing any module.

**`wave.strategy` when processes declare both `conda` and `container`:** nf-core modules typically specify both a `conda` directive and a `container` directive. By default Wave resolves containers in priority order `container,dockerfile,conda` — so it uses the `container` image and ignores `conda`. If you want Wave to build containers from the `conda` directive instead (e.g. to get a custom Conda environment), set:

```groovy
wave.strategy = ['conda']
```

This instructs Wave to use the `conda` directive and ignore `container` and Dockerfile. See https://docs.seqera.io/nextflow/latest/wave#build-conda-based-containers for details. You can also choose a build template for smaller images:

```groovy
wave.strategy = ['conda']
wave.build.template = 'conda/pixi:v1'   // multi-stage, ~30-50% smaller (Nextflow 26.04+)
```

The linter enforces the `container` directive (`NF050`). Leaving process environments implicit, or relying only on a conda spec that has to be re-solved per task, is a hard error.

## Compute-cost tuning

Every module process must expose a handle for operators to tune compute cost without editing the module. That means **at least one of**:

- a `label` directive (the nf-core convention: `process_single`, `process_low`, `process_medium`, `process_high`, `process_high_memory`, `process_long`) that a config file maps to concrete resources, or
- explicit resource directives (`cpus`, `memory`, `time`, `machineType`, `accelerator`) on the process itself.

The label approach is strongly preferred because it centralises cost control in one place (the shared process configuration), so a site-specific override can scale every `process_high` tier at once without touching any module.

```nextflow
// compute_kmer_profile.nf
process COMPUTE_KMER_PROFILE {
    label 'process_low'
    container 'quay.io/nf-core/ubuntu:20.04'
    ...
}
```

```groovy
// Shared process configuration
process {
    withLabel: 'process_single' { cpus = 1;  memory = 1.GB;  time = 30.min }
    withLabel: 'process_low'    { cpus = 2;  memory = 6.GB;  time = 2.h    }
    withLabel: 'process_medium' { cpus = 6;  memory = 24.GB; time = 6.h    }
    withLabel: 'process_high'   { cpus = 12; memory = 64.GB; time = 12.h   }
}
```

On a cloud executor (AWS Batch, Google Batch, Azure), add `machineType` per label in a site-specific config so each tier lands on an appropriately-sized instance. The linter enforces the declaration (`NF052`); choosing the right tier per module is a human call.

### `resourceLimits`: cap requests to available infrastructure

The `resourceLimits` directive (Nextflow 24.04+) caps resource requests at specified maximums. If a process (or a dynamic retry expression) requests more than the limit, Nextflow silently caps it — the task still runs, just at the limit. Without this, a `memory { 8.GB * task.attempt }` retry on a machine with 64 GB max would fail on attempt 9 with an impossible request.

```groovy
// nextflow.config — global limits matching your infrastructure
process {
    resourceLimits = [
        cpus: 72,
        memory: 512.GB,
        time: 168.h
    ]
}
```

Set these to match the largest instance type (or node class) available in your compute environment:

- **AWS Batch Forge** — match the `maxCpus` and the memory of the largest allowed instance family
- **Google Batch** — match the largest machine type in the allowed list
- **HPC / Slurm** — match the max allocatable resources per job on your partition
- **Local laptop** — match your physical hardware (e.g. `cpus: 8, memory: 16.GB`)

Use profiles to vary limits per environment:

```groovy
profiles {
    local {
        process.resourceLimits = [cpus: 8, memory: 16.GB, time: 24.h]
    }
    cloud {
        process.resourceLimits = [cpus: 72, memory: 512.GB, time: 168.h]
    }
}
```

**Interaction with retries:** The standard nf-core pattern multiplies resources by `task.attempt`:

```groovy
process {
    withLabel: 'process_high' {
        cpus   = { 12    * task.attempt }
        memory = { 64.GB * task.attempt }
    }
}
```

Without `resourceLimits`, attempt 3 requests 192 GB — which may exceed your largest node and cause an immediate scheduling failure. With `resourceLimits = [memory: 128.GB]`, attempt 3 is silently capped to 128 GB and can still run.

**Where to declare it:** Always in config (not in module files). A single global declaration in `nextflow.config` or the shared process configuration is usually enough. If needed, per-label overrides are valid too.

## Parameters and configuration

Parameter handling should be boring and standardized.

When naming parameters or deciding where they should be declared and validated, read:

- https://docs.seqera.io/nextflow/strict-syntax
- https://nf-co.re/docs/guidelines/pipelines/requirements/parameters

### Rules for parameters

- Define defaults once.
- Use standard CLI names such as `--input`, `--outdir`, and `--genome` where they fit.
- Pass parameters from the entry workflow into downstream workflows and processes as explicit inputs.
- Use `-params-file` for large launch surfaces.
- Every user-facing pipeline ships a `nextflow_schema.json` next to `main.nf` that defines every `params.*` key, its type, its default, and a description. The schema is the source of truth: it drives `--help`, validates user input before the first task is submitted, and documents the parameter surface. The linter enforces its presence (`NF051`).

What to avoid:

- reading `params.*` from modules
- inventing new parameter names for common concepts
- building complex tool argument strings inside `main.nf`

### Schema-driven validation with `nf-schema`

Use the `nf-schema` plugin for parameter validation. It reads `nextflow_schema.json` at launch, fails loudly if required parameters are missing or types are wrong, and auto-generates `--help` output.

Wire it in `nextflow.config`:

```groovy
plugins {
    id 'nf-schema@2.7.2'
}
```

Then call `validateParameters()` at the top of the entry workflow:

```nextflow
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {
    validateParameters()
    log.info paramsSummaryLog(workflow)
    // ...
}
```

Schema format: JSON Schema 2020-12 (`"$schema": "https://json-schema.org/draft/2020-12/schema"`), with parameter groups in `$defs` referenced from a top-level `allOf`. Do **not** use the older draft-07 format (`definitions`, `http://json-schema.org/draft-07/schema`) — nf-schema 2.x requires 2020-12.

### `publishDir` vs `storeDir`

These are different directives with different semantics — do not confuse them:

- **`publishDir`** — copies or symlinks outputs to a user-visible results directory after the task completes. Use this for final outputs you want the user to find. The working directory is still the source of truth for caching; `publishDir` is a side-effect.
- **`storeDir`** — uses an external directory as the task's backing store. If the expected outputs already exist there, the task is **skipped entirely** (like a persistent cache). Use this for expensive downloads or reference-data prep steps that should never re-run.

```nextflow
// publishDir — always runs; copies results afterward
process ALIGN_READS {
    publishDir params.outdir, mode: 'copy'
    ...
}

// storeDir — skips the task if outputs already exist at the path
process DOWNLOAD_REFERENCE {
    storeDir "${params.genome_cache_dir}"
    ...
}
```

Misusing `storeDir` for regular outputs means Nextflow can never detect when an input changed and must re-run the task.

**Avoid `cleanup = true` in development.** The top-level `cleanup` config option deletes all work-directory files for completed tasks when a run finishes successfully. This breaks `-resume` on subsequent runs and is not supported for remote work directories (S3, GCS, Azure). Use `nextflow clean` for selective cleanup instead.

### Workflow outputs: the replacement for `publishDir`

For new pipelines targeting Nextflow 26.04+, prefer **workflow-level outputs** over per-process `publishDir`. This centralizes all output publishing in one declarative block and works with channels rather than glob patterns. The entry workflow declares what to publish; a top-level `output {}` block declares where and how:

```nextflow
workflow {
    main:
    ALIGN(ch_reads)
    QUANT(ALIGN.out.bam)

    publish:
    alignments = ALIGN.out.bam
    counts     = QUANT.out.counts
}

output {
    directory params.outdir

    alignments {
        path params.alignments_destination
        mode 'copy'
    }
    counts {
        path params.quantification_destination
        index {
            path 'counts.csv'
            header true
        }
    }
}
```

Benefits over `publishDir`:
- **Single source of truth** — all publishing decisions live in one place, not scattered across config `withName:` blocks.
- **Index files** — automatic CSV/JSON/YAML manifests of published outputs with metadata, usable as samplesheets for downstream pipelines.
- **Dynamic paths** — a `path` closure using user-selected destinations and sample identifiers closures route files per-sample without relying on filename parsing.
- **Conditional publishing** — `sample.bam >> (params.save_bams ? params.alignments_destination : null)` to skip outputs based on params.
- **Data lineage integration** — workflow outputs feed directly into the Nextflow lineage store and can be labeled for cross-run queries.

The config-driven `publishDir` in `withName:` blocks (as used by existing nf-core pipelines) remains valid and is the current nf-core convention. Use workflow outputs for new standalone pipelines; follow whatever convention the pipeline already uses when contributing to existing code.

See https://docs.seqera.io/nextflow/workflow for the full syntax.

### Make published reports visible in Nextflow Platform

For a new pipeline, provide a `tower.yml` that exposes its primary
human-readable outputs in the run's Reports tab. Platform can only list files
that the pipeline actually publishes, so each report pattern must match the
published layout exactly:

```yaml
reports:
  "multiqc_report.html":
    display: "MultiQC report"
  "*_summary.tsv":
    display: "Sample summary"
    mimeType: "text/tab-separated-values"
```

- Include at least the pipeline's main HTML, PDF, CSV, TSV, or TXT reports by
  default; do not wait for the user to ask after the first run.
- Match published paths or filenames, not task work-directory paths and not
  `${params.outdir}` expressions. Quote glob patterns containing `*`.
- Keep `display` user-facing and specific. Omit `mimeType` when the file
  extension is sufficient.
- Before declaring the pipeline ready, run it with representative data, confirm
  the report files exist in the published output tree, and check every
  `tower.yml` pattern against those paths. After a Platform launch, verify the
  workflow reports endpoint returns the expected entries.

Do not list arbitrary large data files merely to populate the tab. Report
previews support HTML, PDF, CSV, TSV, and TXT; other outputs remain available
from their published storage paths.

## Use Nextflow primitives before reaching for Python

A recurring anti-pattern is writing Python for tasks that Nextflow already does well.

Use Nextflow primitives for:

- joining tables or metadata
- grouping samples
- choosing branches
- mixing channels
- carrying tuples and meta maps
- passing directories as paths

Do not:

- serialize everything to JSON
- unpack directories into many individual files unless the tool truly requires that
- copy staged files manually into the task directory
- write long inline Python blocks for simple channel logic

If two lines of bash or a native operator solve the problem, use that.

### Do this

Use workflow logic and channel operators for glue code:

```nextflow
ch_fold_inputs = ch_raw_sequences
    .map { meta, fasta -> tuple(meta, fasta) }

if( params.model_name == 'boltz' ) {
    BOLTZ_PREDICT(ch_fold_inputs)
}
else if( params.model_name == 'openfold3' ) {
    OPENFOLD3_PREDICT(ch_fold_inputs)
}
else {
    error "Unsupported model: ${params.model_name}"
}
```

Pass directories as actual path inputs when the tool consumes a directory:

```nextflow
process SCORE_MSA_DIRECTORY {
    input:
    tuple val(meta), path(msa_directory)

    output:
    tuple val(meta), path("scores.json"), emit: scores

    script:
    """
    score_msa_directory \
        --msa-dir ${msa_directory} \
        --output scores.json
    """
}
```

Use an available project helper command directly when a small helper script is genuinely needed:

```nextflow
process NORMALIZE_MANIFEST {
    input:
    path(manifest_csv)

    output:
    path("normalized_manifest.csv"), emit: manifest

    script:
    """
    normalize_manifest.py \
        --input ${manifest_csv} \
        --output normalized_manifest.csv
    """
}
```

### Do not do this

Do not write Python, or any other non-Nextflow language, for channel reshaping, branch selection, or other glue logic that Nextflow can express directly:

```nextflow
process PREPARE_FOLD_INPUTS {
    input:
    tuple val(meta), path(fasta)

    output:
    path("prepared.json"), emit: prepared

    script:
    """
    python3 - <<'PY'
    import json

    payload = {
        "sample_id": "${meta.id}",
        "fasta": "${fasta}",
        "model": "${params.model_name}",
    }

    with open("prepared.json", "w") as handle:
        json.dump(payload, handle)
    PY
    """
}
```

Do not unpack directories into many files if the tool accepts a directory:

```nextflow
process SCORE_MSA_FILES {
    input:
    tuple val(meta), path(a3m_file), path(hhm_file), path(cs219_file), path(index_file)

    script:
    """
    score_msa_directory \
        --a3m ${a3m_file} \
        --hhm ${hhm_file} \
        --cs219 ${cs219_file} \
        --index ${index_file}
    """
}
```

Do not copy, rewrite, or recreate an input file inside the task script when that file was already declared in `input:` as a `path`.

If a file appears in `input:`, Nextflow stages it into the task working directory and provides the correct path through the input variable. The script should consume that staged path directly.

Wrong:

```nextflow
script:
"""
python3 -c "from pathlib import Path; Path('input.fa').write_text(Path('${fasta}').read_text())"
run_tool --input input.fa
"""
```

Correct:

```nextflow
input:
path(fasta)

script:
"""
run_tool --input ${fasta}
"""
```
