---
name: build-nextflow-pipeline
description: >
  Plan and scaffold a new Nextflow DSL2 pipeline from source material —
  scripts, repositories, Jupyter/R notebooks, papers, or informal
  instructions. Drives the planning phases before any .nf code is written:
  mapping the data flow, critiquing it, shaping channels and metadata,
  scoping subworkflows, hunting for existing modules, and only then building
  containers and writing code. Use this skill whenever the user asks to
  "build a pipeline", "compose registry modules", "port this to Nextflow", "turn this notebook/repo/script
  into a pipeline", or otherwise wants to go from an existing body of work to
  a clean Nextflow pipeline. Pair with `nf-pipeline-design`, which owns the
  code-level rules (layout, main.nf, subworkflow/module shape); this skill
  owns the planning process that precedes them.
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


# Build a Nextflow Pipeline

## Choose the construction path

For composing existing Registry modules, read [registry composition](references/create-workflow/README.md) and follow its plan → approval → module validation → composition → end-to-end verification loop. For a single module, use `run-module`; a wrapper workflow is unnecessary.

For scripts, notebooks, repositories or analysis descriptions, use the source-material phases below. Load `nf-pipeline-design` for code-level rules in either path.

This skill guides the **planning** of a Nextflow pipeline that is being built from something that already exists: a set of scripts, a GitHub repository, a Jupyter or R notebook, a methods section in a paper, or a half-formed description of an analysis.

**Read `nf-pipeline-design` first.** That skill defines what clean Nextflow code looks like (`main.nf` shape, subworkflow/module boundaries, operators vs. modules, containers). This skill defines the *process of getting there* — everything that must happen **before and during** the decision to write a single `.nf` file. The two are meant to be used together: plan with this skill, write with the other.

The phases below are ordered. Do not skip ahead. The single most common failure mode when translating existing work into Nextflow is writing modules before the data flow is clear — the resulting pipeline looks like the original script transcribed into processes, not a pipeline designed for parallel execution and reuse.

## Delegate context-free steps to subagents

Several steps in this flow do not need the full planning conversation to do their job — they just need a tight input (source material, a tool name, an analysis step) and return a structured artifact. Hand those steps off to subagents in parallel whenever you can. It is faster, keeps the main conversation focused on decisions that actually require context, and lets you fan out over many tools or analysis steps at once.

Delegate only when the host supports subagents; otherwise perform the same steps sequentially using the companion skills below. Include the relevant skill instructions and tight input in each delegated task.

| Phase step | Delegate? | Companion skill |
|---|---|---|
| Auditing what the source material actually points at | **Yes** — hand the source material to one subagent | `audit-conversion-readiness` |
| Mapping the data flow | No — this is the core act of planning, keep it in the main conversation | — |
| Triaging parameters and drafting `nextflow_schema.json` | **Yes** — hand the source material to one subagent | `nextflow-schema` |
| Designing channels and metadata | No — needs the data flow map as live context | — |
| Comparing alternatives, only when requested or the selected tool is blocked | **Optional** — delegate the bounded comparison when useful | `find-alternative-tools` |
| Searching for existing Nextflow modules for each chosen tool | **Yes — one subagent per tool, in parallel** | `search-existing-modules` |
| Building and verifying the container for each tool branch | **Yes — one subagent per tool, in parallel** | `create-container` |
| Writing `main.nf`, subworkflows, and modules | No — this is a synthesis step that uses every prior artifact | — |

The main conversation's job is to orchestrate: produce the data flow map with the user, spawn subagents with the right skills and inputs, collect their artifacts, integrate them, and keep the planning document coherent. Subagents do the bounded, parallelizable research and scaffolding.

When delegating, each subagent prompt should contain:

- A clear statement of the step's bounded task (paraphrased from the skill).
- The tight input (paths to source material, the tool or analysis step name).
- Any constraints the main conversation has already committed to (target executor, license restrictions, must-use-this-container).

Do not include the full conversation history. The subagent skills are designed to work from the tight input alone — that is the point.

## Planning artifacts

Keep planning artifacts in the user's chosen project locations, or provide them
in the conversation when file access is unavailable. Respect existing project
organization. The useful artifacts are:

- **the planning document** — the living planning document. One markdown file with sections for: data flow map, samplesheet design, channel shapes, subworkflow scoping with chosen tools per step, alternative-tool rationale, test-inputs table, container verification status. Every phase adds to it; subagent outputs get pasted or linked into the relevant section.
- **`nextflow_schema.json`** — the parameter contract. First-pass draft from Phase 1 (produced by `nextflow-schema`), augmented in Phase 3 with approved implemented tool choices (if any) and relevant tool arguments revealed by module search. Always the source of truth for what the pipeline accepts.
- **`nextflow.config`** — mirrors schema defaults, declares profiles, loads `modules.config`. First-drafted in Phase 1 alongside the schema; stays in sync with it.
- **`modules.config`** — plumbs surfaced parameters into `ext.args` per module. Drafted in Phase 3 once tools are picked.
- **`tower.yml`**, when Platform report integration is requested, maps the pipeline's primary published human-readable outputs to Nextflow Platform report entries. Draft it from the planned output layout and verify every pattern against a real run.
- **the sample-sheet schema** — samplesheet schema, validates each row of the samplesheet at launch. Separate from `nextflow_schema.json` (which validates parameters, not rows). Drafted in Phase 1 once "what is a sample" is settled.
- **the user's selected test data** — small, representative test inputs per tool. Identified during Phase 1 (the triage subagent returns a test-inputs manifest as a secondary artifact) and consumed by `create-container` subagents in Phase 4.

Treat `nextflow_schema.json` as a **living artifact** across phases. The Phase 1 draft is not the final version — each later phase may add new parameters (`enum` values for approved implemented alternatives, new `ext.args` knobs after module search), and the main conversation is responsible for keeping the schema, `nextflow.config`, and the planning document in sync as those additions come in.

## Conversion contract

For Python, R or Jupyter source, read the matching conversion playbook before scaffolding. Audit external files, helper code, dependencies and notebook state. Run the source on a small representative fixture when possible and record expected scientific outputs before changing execution structure. If the source cannot run, state that equivalence is unverified.

Preserve tested algorithms by wrapping them. Split only at meaningful independently executable boundaries; file loading or a notebook cell alone does not justify a process. Preserve within-task threading and global/statistical operations. Validate the converted pipeline against the baseline outputs using exact or format-aware/tolerance-based checks appropriate to the algorithm. Compilation and lint alone do not establish equivalence.

Use output syntax supported by the selected Nextflow version and the project's publication contract; legacy examples using `publishDir` are not a reason to undo workflow outputs.

## Phase 0 — Audit what the source material points at

Source material routinely references data and helper code that is not actually
available: absolute paths into a lab NAS, reference bundles nobody can reach from
here, a helper script that never made it into the handover. Finding that out after
the modules exist is the most expensive way to learn it.

Run `audit-conversion-readiness` on the source material first. It
returns a blockers-first report — unreachable paths, missing tool code, tools with
no distributable container — plus one consolidated list of what to ask the user.
Fold the surviving references into Phase 1's data flow map and parameter triage.

## Phase 1 — Map the data flow

Before anything else, produce an explicit map of how data moves through the intended pipeline. This map is the single most important artifact of the planning phase: every later decision depends on it, and errors here propagate everywhere.

The map should answer, concretely:

- **What is the input data?** Not "sequences" — *which* sequences, in which format, from where, at what scale.
- **What are the pipeline parameters?** Every knob in the source material that has an observable effect on output is a candidate parameter. Triage them into user-facing vs. internal constants, then draft the initial `nextflow.config` with default values (see "Design the parameter surface" below).
- **What is the first subworkflow?** What does it accomplish? Why is it there (what invariant does it establish for the rest of the pipeline)?
- **Same for every subsequent subworkflow**, in order. Each one needs a one-sentence rationale. If you can't articulate why a subworkflow exists, it probably shouldn't.
- **What are the outputs?** Files, reports, derived data — and who consumes them.

### The samplesheet is the input contract

Pipelines should be driven by a **samplesheet** — a CSV or TSV whose rows correspond to samples, where a *sample* is defined as the unit the pipeline processes in parallel. The samplesheet is how Nextflow parallelizes work; without it, you don't have a pipeline, you have a script with extra steps.

When mapping the flow, decide early:

- What is a "sample" in this pipeline? (one FASTQ pair, one patient, one protein, one image stack…)
- What per-sample metadata is needed? (sample ID, condition, reference genome, model name…)
- Which inputs are per-sample vs. global?

Global inputs (a reference genome, a model checkpoint, a blastdb) are parameters, not samplesheet columns.

Once "what is a sample" is settled, draft **the sample-sheet schema** — a JSON schema that validates each row of the samplesheet. nf-schema uses this file to validate the CSV at launch (e.g. required columns, allowed values, file-path patterns). Keep it distinct from `nextflow_schema.json`, which validates pipeline parameters, not samplesheet rows.

### Design the parameter surface

A Nextflow pipeline's parameters are its public API. Every time someone runs the pipeline, these are the knobs they see. Getting them right at the planning stage is much cheaper than renaming them after downstream configs, tests, and docs exist.

Go through the source material (scripts, notebook cells, CLI flags in the original tool, hard-coded constants, tool arguments buried in helper functions) and list **every** value that has an observable effect on the output. Then triage:

- **Surface as a parameter** — anything the user would reasonably want to change: inputs, output directory, thresholds that affect results, tool/algorithm choices, reference genome paths, **and any tool flag whose value matters enough to be discussed in the original paper, README, or notebook** (e.g. k-mer size, model checkpoint, significance threshold). These get declared at the top level so they are discoverable, validated, and overridable from the CLI.
- **Pin as a true internal constant** — defaults that are genuinely never touched and are part of the tool's own contract (e.g. a fixed random seed for reproducibility, an internal file format version). These can stay hardcoded in `ext.args`.
- **Expose as a profile, not a param** — environment-dependent choices (executor, container engine, test inputs) belong in `profiles {}`.

When in doubt, surface. A parameter that turns out to be unused is trivial to remove; a constant that turns out to matter requires a breaking change to promote.

**Delegate this step when available.** Parameter triage is exhaustive, source-driven, and does not need the planning conversation's context — only the source material. Use the `nextflow-schema` skill with that source material. Produce four artifacts: the triaged parameter list, the drafted `nextflow_schema.json` (first pass), the `nextflow.config` skeleton, **and a test-inputs manifest** listing a minimal representative test input per heavy tool. Copy those artifacts into the repo and into the planning document's relevant sections, then review the uncertainty flags at the top before accepting.

The test-inputs manifest is the spec that `create-container` subagents will consume in Phase 4 — every heavy tool needs a row. If the triage subagent could not locate a test input for some tool, flag it now; do not leave it for Phase 4 to discover mid-run.

### `nextflow_schema.json` is the top-level source of truth

All surfaced parameters — including the ones that feed into `ext.args` — go in **`nextflow_schema.json`** in the user's pipeline. This is the canonical declaration of the pipeline's parameter contract: types, defaults, descriptions, allowed values, required fields, grouping for docs and launch UIs. `nextflow.config` just mirrors the defaults for runtime; the schema owns the contract.

Draft `nextflow_schema.json` *before* writing any `.nf` code — but understand this as a **first pass**. Later phases will add `enum` values (once alternative tools are enumerated in Phase 3) and occasionally new tool-argument parameters (revealed during module search). The schema is a living artifact across phases; the Phase 1 draft is the skeleton, not the finished version. The act of writing types, enums, and descriptions forces decisions you'd otherwise punt on, and produces documentation for free.

When unsure about schema syntax, params patterns, or how the schema drives validation and launch UIs, read:

- nf-schema plugin and schema specification: https://nextflow-io.github.io/nf-schema/latest/nextflow_schema/nextflow_schema_specification/
- Main config reference: https://docs.seqera.io/nextflow/config
- Strict-syntax rules for `params`: https://docs.seqera.io/nextflow/strict-syntax
- nf-core parameter conventions (names, defaults, `--input`/`--outdir`/`--genome`): https://nf-co.re/docs/guidelines/pipelines/requirements/parameters
- nf-core pipeline template (reference layout including a real `nextflow_schema.json`): https://github.com/nf-core/tools/tree/master/nf_core/pipeline-template

#### Minimal `nextflow_schema.json` shape at the planning stage

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "title": "my-pipeline parameters",
  "type": "object",
  "$defs": {
    "input_output_options": {
      "title": "Input/output options",
      "type": "object",
      "required": ["input", "outdir"],
      "properties": {
        "input": {
          "type": "string",
          "format": "file-path",
          "pattern": "^\\S+\\.csv$",
          "description": "Path to the samplesheet CSV describing the samples to process."
        },
        "outdir": {
          "type": "string",
          "format": "directory-path",
          "description": "The output directory where results will be published."
        }
      }
    },
    "pipeline_knobs": {
      "title": "Pipeline-specific options",
      "type": "object",
      "properties": {
        "model_name": {
          "type": "string",
          "default": "boltz",
          "enum": ["boltz", "openfold3"],
          "description": "Folding model to use for structure prediction."
        },
        "min_sequence_length": {
          "type": "integer",
          "default": 30,
          "minimum": 1,
          "description": "Sequences shorter than this are filtered out before folding."
        }
      }
    },
    "tool_arguments": {
      "title": "Important tool arguments (surfaced from defaults)",
      "type": "object",
      "properties": {
        "msa_kmer_size": {
          "type": "integer",
          "default": 15,
          "description": "K-mer size passed to the MSA search tool. Exposed because the original method varies this per organism."
        },
        "folding_recycles": {
          "type": "integer",
          "default": 3,
          "minimum": 1,
          "description": "Number of recycling iterations for the folding model."
        }
      }
    },
    "reference_data": {
      "title": "Reference data",
      "type": "object",
      "properties": {
        "reference_genome": { "type": "string", "format": "file-path" },
        "blastdb":          { "type": "string", "format": "directory-path" }
      }
    }
  },
  "allOf": [
    { "$ref": "#/$defs/input_output_options" },
    { "$ref": "#/$defs/pipeline_knobs" },
    { "$ref": "#/$defs/tool_arguments" },
    { "$ref": "#/$defs/reference_data" }
  ]
}
```

#### The companion `nextflow.config`

`nextflow.config` mirrors the schema defaults in its `params {}` block, declares profiles, and wires the surfaced tool parameters into `ext.args` so modules stay parameter-free:

```groovy
// nextflow.config
plugins {
    id 'nf-schema'
}

params {
    // mirrors nextflow_schema.json defaults — keep in sync
    input                = null
    outdir               = null
    model_name           = 'boltz'
    min_sequence_length  = 30
    msa_kmer_size        = 15
    folding_recycles     = 3
    reference_genome     = null
    blastdb              = null
}

profiles {
    test {
        // Supply input and output values from the user-selected test configuration.
    }
    docker      { docker.enabled      = true }
    singularity { singularity.enabled = true }
}

// Include the user's module configuration using its actual project reference.
```

And `modules.config` plumbs surfaced parameters into module invocations — the module itself does not read `params`:

```groovy
// modules.config
process {
    withName: 'RUN_MSA' {
        ext.args = { "--kmer-size ${params.msa_kmer_size}" }
    }
    withName: 'FOLDING_MODEL_INFERENCE' {
        ext.args = { "--recycles ${params.folding_recycles}" }
    }
}
```

Rules of thumb for this first draft:

- **`nextflow_schema.json` is the source of truth.** Every surfaced parameter lives there with a type, a default, and a description. `nextflow.config` mirrors the defaults; `modules.config` consumes them via `ext.args`.
- **Important tool flags are surfaced, not hidden.** If a value is discussed in the paper/README or if changing it measurably affects output, it is a top-level parameter — even though it ultimately flows through `ext.args`.
- **Defaults or `null` for every field.** Use `null` (and list in `required`) for inputs the pipeline cannot run without, so validation fails loudly at launch.
- **Pipeline-meaningful names** — `protein_sequence_to_be_folded`, not `fasta` (see `nf-pipeline-design` on naming). Align with nf-core conventions where they exist: `--input`, `--outdir`, `--genome`.
- **Group parameters in `$defs`** by user intent (I/O, pipeline logic, tool arguments, reference data). This grouping becomes the structure of the launch form and the docs.
- **Modules never read `params`.** Surfaced tool parameters reach modules only via `ext.args` configured in `modules.config`, as shown above.

### Critique the map

Once the data flow is drafted, **stop and critique it** — ideally out loud with the user. Ask:

- Is any step doing more than one thing? Can it be split?
- Are there steps that only exist because the original script was linear, but could run in parallel in a real pipeline?
- Are there steps that could be collapsed into a single module because they are trivially small and only make sense together?
- Is the samplesheet really the right unit of parallelism, or should fan-out happen at a different level?

This is the most important step. A wrong map here will produce a pipeline that costs more to refactor later than it cost to write the first time.

## Phase 2 — Design channels and metadata

Only once the data flow map is accepted should you start thinking in terms of channels. Write the channel design into the planning document's "Channels" section as you go — it is what Phase 3 reads when scoping subworkflows, and what Phase 5 reads when writing the actual `.map { }` / `.join(...)` calls.

Without writing any subworkflow yet, outline:

- Which channels are needed? What is the **shape** of each (`tuple val(meta), path(x)`, `path(x)`, `val(x)`)?
- Which channels need to **merge** (`join`, `combine`, `mix`)? On what key?
- Which steps run **per sample** (fan-out) vs. **across all samples** (fan-in, e.g. a joint report)?
- Where does metadata get attached, augmented, or consumed?

Draw this. A short ASCII diagram or bullet list of channel transitions is enough — the goal is that someone reading your plan can predict what `.map { ... }` / `.join(...)` calls will appear in the subworkflows.

### The cost-of-a-process rule

Submitting a Nextflow process has a real cost: scheduling overhead, possible queue wait, container pull, staging. For tiny read/write operations — reading a header, renaming a file, reshaping a CSV, merging metadata — a spawned task is almost always the wrong choice. Use this priority order:

1. **A Nextflow operator** (`map`, `filter`, `join`, `branch`, `combine`, `collect`, `flatMap`, `groupTuple`, …). First choice. No task submitted, runs in the head job.
2. **A Groovy function** at the subworkflow level. Use when one operator isn't quite enough but the logic is still cheap and deterministic.
3. **A single module that processes many files at once**, with fan-in at the input and fan-out at the output. Use when the operation is too awkward for operators/Groovy but amortizes across inputs.
4. **A module that processes each file individually.** Last resort for cheap ops — reserve this pattern for the expensive, genuinely-parallelizable work it was designed for.

The inverse rule matters too: for **large** operations (alignment, folding, variant calling, training loops — the actual bioinformatics workload), parallelize aggressively using channel operations and per-sample modules. Those are the jobs Nextflow exists to schedule. Do not batch them manually into a single monster task to "keep things simple."

## Phase 3 — Scope the subworkflows

With channel shapes pinned down, sketch the subworkflows. At this stage you are still not writing `.nf` code — you are deciding the *boundary* of each subworkflow.

Each subworkflow should:

- Express **one analysis step** (e.g. "run the MSA", "call variants", "build the report"), not a vague bucket.
- Execute **a single module per execution path**, even if multiple modules are interchangeable (choose between them with workflow-level `if/else`, as `nf-pipeline-design` describes).
- Be designed in clear **sections**: channel reshaping → module call → optional output reshaping. Readability matters.
- **Harmonize inputs** — much of the value of a subworkflow is turning whatever upstream shape you got into exactly what the underlying tool needs, so the module stays boring.

### Select tools without expanding scope

Use the tool requested by the user or required to preserve the source analysis. If none is specified, recommend one suitable default and resolve material license, compute or output-contract constraints before implementation.

Use `find-alternative-tools` only when the user requests a comparison or the selected tool is blocked. Research produces candidates; it does not authorize implementing them.

Add multiple tool branches only when the user explicitly requests or approves them. For approved branches:

1. Add only implemented choices to the schema's `enum`, keeping config defaults aligned.
2. Put tool selection in the subworkflow and reject unsupported choices explicitly.
3. Normalize inputs and outputs around the selected modules.
4. Verify every implemented branch's command, container and representative outputs.

Record the selected tool and rationale. A single-tool pipeline does not need an alternative-tool parameter or exhaustive rejected-tool inventory.

### Search for existing modules before writing one

For each selected tool, search for existing implementations in this priority order. This is not optional: reinventing an nf-core module costs days and produces worse code.

1. **nf-core modules.** Search https://nf-co.re/modules and the `nf-core/modules` GitHub repo. If a module exists, use it (or vendor it). If one exists but doesn't quite fit, prefer extending it over writing a new one.
2. **The inputted source.** Check the scripts, notebooks, or repository you were given — is there already a containerized version, a Dockerfile, a conda env? Reuse it.
3. **Online: papers, preprints, GitHub, biorxiv, arxiv.** Look for community implementations, reference pipelines, or tool authors' recommended invocations. A five-minute search often finds a canonical invocation that saves hours of trial and error.

**Delegate this search when available.** For each tool selected, use the `search-existing-modules` skill with the tool name, the path to the pipeline's source material, and — if the prior enumerate step flagged one — the nf-core module hint. The hint enables a quick verification instead of re-running the full priority-order search. These tasks can run in parallel with one subagent per tool. Each search produces a verdict (reuse / vendor / write new) plus the structured metadata needed by the next phase. Paste the verdicts into the planning document's "Module sourcing" section.

If a verdict reveals tool-argument parameters that were not captured in Phase 1's schema draft (e.g. the canonical invocation uses a `--kmer-size` flag that the source material omitted), add them to `nextflow_schema.json` now — this is the second expected schema augmentation. Keep `nextflow.config` in sync.

## Phase 4 — Build and test the container before writing module code

Before writing a single line of module code, **build the container and run the tool inside it end-to-end on a representative input.** Use [Wave](https://docs.seqera.io/wave) to build from a conda env or Dockerfile.

The loop looks like:

1. Declare the software environment (conda YAML, Dockerfile, or image reference).
2. Build with Wave.
3. Shell into the container (or run it interactively) and execute the exact command line you plan to put in the module's `script:` block, on a small but real input.
4. Iterate until the command works cleanly — correct outputs, correct exit code, no `stderr` surprises, no silent empty files.

Only then translate that working command into a Nextflow module. This flips the usual debugging loop: instead of discovering tool-level problems through the Nextflow harness (slow, noisy, hard to iterate on), you catch them at the shell level first. The module becomes a thin, predictable wrapper around a command you already know works.

Why this ordering matters: most "Nextflow bugs" when porting existing work are not Nextflow bugs at all — they're environment bugs (missing dependency, wrong version, missing reference data, different glibc) that happened to surface inside a process. Separating the two is the single biggest time-saver in the build loop.

**Delegate per-tool container builds when available.** For each tool that came back as "vendor" or "write new" from the module search, use the `create-container` skill with:

- The tool name and target version.
- The target command line (from the `CANONICAL CLI` field of the search-existing-modules verdict).
- The test-input row for that tool from the planning document's "Test inputs" section (produced by `nextflow-schema` in Phase 1). Include path, format, size, and expected output shape.
- Any constraints the main conversation has committed to (executor, GPU, container registry).

Run these **in parallel** — one per tool — since each is fully independent. Each subagent returns a verified image reference and the exact working command, which you paste into the planning document's "Container verification" section ready for the module-writing step.

For tools that came back as "reuse nf-core module", skip the container step — nf-core modules come with their containers already pinned. Only re-verify if the module is old enough that you suspect the image tag has rotted.

## Phase 5 — Write the code

At this point, hand off to `nf-pipeline-design`. You have:

- A data flow map.
- A channel-shape plan.
- A subworkflow scoping document with chosen tools.
- Working containers for each heavy step.

Writing the actual `main.nf`, subworkflows, and modules should now be a mechanical translation of the plan into Nextflow code that obeys the rules in `nf-pipeline-design`. If it isn't mechanical — if you find yourself making new design decisions in the middle of writing a subworkflow — stop, go back to the relevant planning phase, and update the plan first.

## Checklist before writing any `.nf` file

- [ ] the planning document exists and has sections filled for: data flow map, samplesheet design, channels, selected tools (alternatives only if requested), module sourcing, test inputs, container verification
- [ ] Data flow map has been critiqued, not just drafted
- [ ] "Sample" is defined and the sample-sheet schema is drafted
- [ ] `nextflow_schema.json` is available with defaults and descriptions; any tool-choice `enum` lists only approved implemented branches
- [ ] `nextflow.config` mirrors the schema defaults and declares profiles
- [ ] `modules.config` plumbs surfaced parameters into `ext.args` per module
- [ ] Channel shapes are listed for every major transition in the planning document
- [ ] Cheap operations are assigned to operators/Groovy; expensive ones to modules
- [ ] Each subworkflow has a one-sentence rationale and a named underlying tool (or set of alternative tools with branching)
- [ ] The selected tools match the requested analysis; any additional implemented branches were explicitly requested or approved
- [ ] Existing-module search has been done for every chosen tool, using `search-existing-modules`, with verdicts in the planning document
- [ ] Test inputs exist at the user's selected test data (or are specified to be derived) for every heavy tool
- [ ] Containers are built and tool commands verified, delegated to `create-container` subagents where possible, with verified image references in the planning document
- [ ] `nf-pipeline-design` has been re-read for the code-level rules

## Checklist before declaring the pipeline ready

- [ ] The pipeline has run end-to-end with representative data
- [ ] Primary report files exist in the published output tree
- [ ] `tower.yml` patterns match the user's actual published reports
- [ ] After a Platform launch, the workflow Reports tab lists the expected entries

## Registry composition playbook

When composing multiple Nextflow Registry modules, use the registry path instead of the source-material planning phases. Read [registry composition](references/create-workflow/README.md) before proceeding.

## Python conversion playbook

When converting a standalone Python script or Python analysis into Nextflow, preserve its algorithm and dependency contract. Read [Python conversion](references/convert-python-script/README.md) before proceeding.

## R conversion playbook

When converting R source into Nextflow, preserve statistical semantics, dependencies and meaningful task boundaries. Read [R conversion](references/convert-r-script/README.md) before proceeding.

## Notebook conversion playbook

When converting a Jupyter notebook, inventory cell dependencies and hidden state before extracting executable steps. Read [notebook conversion](references/convert-jupyter-notebook/README.md) before proceeding.
