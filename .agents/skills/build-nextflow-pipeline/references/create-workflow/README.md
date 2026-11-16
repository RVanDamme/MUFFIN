<!-- Modified for the Nextflow plugin: generic host tools and companion guidance. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../../../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Nextflow Workflow Writer

Create complete Nextflow workflows by composing validated modules from the [Nextflow Registry](https://registry.nextflow.io).

Modules are published under namespaces (e.g. `nf-core/fastqc`), and these skills compose modules from any of them.

**Requires Nextflow 26.04 or later** (for the `nextflow module` commands used during validation). Including workflow modules requires **Nextflow 26.10 or later**.

Registry modules come in two kinds (shown as `Kind:` by `nextflow module view`):
- **Process modules** — one tool, e.g. `nf-core/star/align`
- **Workflow modules** — a reusable subworkflow chaining several modules, e.g. `nf-core/fastq_align_star` (STAR alignment + samtools sort/index/stats)

**NEVER write a wrapper workflow just to run/test a single module.**

❌ **WRONG** - Writing a workflow to run one module:
```groovy
// DO NOT DO THIS - not even when a module run fails!
include { FASTQC } from 'nf-core/fastqc'
workflow { FASTQC(Channel.fromPath(params.input)) }
```

✅ **CORRECT** - Use the `run-module` skill:
```
run-module
```

**The `run-module` skill:**
- Uses `nextflow module search/view` to discover and get proper inputs/parameters
- Runs modules directly via `nextflow module run <namespace>/<module>`
- No wrapper workflow needed

**⚠️ When a module run fails due to missing args:**
- DO NOT write a wrapper workflow as a "fix"
- Instead: run `nextflow module view <module>` to get correct parameters
- Fix the command-line arguments and re-run directly

**Only write a workflow in Step 4** when composing multiple validated modules together.

## 4-Step Workflow Creation Process

**ALWAYS follow this structured process when creating a new workflow:**

### Step 1: Identify Modules and Propose Plan

1. Use `nextflow module search <term>` to find Registry modules for each processing step
2. **Look for workflow modules** that already cover a chain of steps. `nextflow module search` rarely ranks them high enough to appear, so query the registry API with the `kind=Workflow` filter:
   ```bash
   curl -s "https://registry.nextflow.io/api/v1/modules?query=align%20reads%20STAR&kind=Workflow&limit=10" \
     | python3 -c 'import json,sys; [print(r["name"], "-", r["description"]) for r in json.load(sys.stdin)["results"]]'
   ```
   When a workflow module matches a chain of steps in the plan, use it instead of hand-wiring the same process modules.
3. Use `nextflow module view <name>` to understand inputs/outputs of each module (for a workflow module: its `take:` inputs in call order and its `emit:` outputs)
4. **Present a plan to the user** with:
   - List of identified modules, marking each as process or workflow module
   - Processing sequence (which module runs first, second, etc.)
   - Data flow between modules (outputs → inputs)

**STOP and wait for user approval before proceeding.**

### Step 2: User Agreement

- Wait for user to review and approve the proposed plan
- Address any questions or modifications requested
- Only proceed when user explicitly agrees

### Step 3: Validate Modules ONE BY ONE (MANDATORY)

**NEVER skip this step.** After user agreement:

1. Determine appropriate **test data** for validation
2. For EACH module in the plan, sequentially:
   - **Install and run with test data**: Invoke `run-module`
   - **Verify outputs** - confirm expected data is produced
   - Only proceed to next module after current one succeeds
   - **Untyped workflow modules** (most nf-core ones) fail `nextflow module run` with "cannot be executed directly because it is not typed". This is expected, not a failure to fix: do not stop to ask the user. Record it in the validation record as "validated in Step 4" and move on. It gets validated by the end-to-end run in Step 4.
3. Record ALL module run commands and their outputs in a validation record selected by the user
4. If a command fails, stop and show the user the command used and the output generated before trying something else

> **Note**: The `run-module` skill uses `nextflow module` commands for discovery, configuration, and execution — modules are installed on-the-fly.

**DO NOT proceed to Step 4 until ALL modules have been individually validated.**

### Step 4: Compose Final Workflow

Only after ALL modules run successfully:

1. Configure `nextflow.config` with Wave + Conda:
   ```groovy
   wave.enabled = true
   wave.strategy = 'conda,container'
   docker.enabled = true
   ```

2. Write the workflow script composing all validated steps using **Nextflow managed modules**, which Nextflow downloads and installs from the Registry:
   ```groovy
   include { MODULE_A } from 'nf-core/module_a'
   include { MODULE_B } from 'nf-core/module_b'

   workflow {
       MODULE_A(input_ch)
       MODULE_B(MODULE_A.out.results)
   }
   ```

   Workflow modules use the same include syntax (`include { FASTQ_ALIGN_STAR } from 'nf-core/fastq_align_star'`) and are called with positional arguments in their `take:` order, as listed by `nextflow module view`. For untyped workflow modules, the input descriptions from `nextflow module view` are approximate. Read the installed `modules/<namespace>/<name>/main.nf` to see which process each input feeds, and build the channel shape that process expects (e.g. `[meta, index]` tuples). On first run, Nextflow installs each one under `modules/<namespace>/<name>/` together with its dependencies in a nested `modules/` directory. Commit the `modules/` directory with the pipeline.

3. **Run the complete workflow using the same test data** to validate end-to-end. This is the validation step for any untyped workflow module skipped in Step 3.

## Critical Guidelines

### Command Execution
- Use `-resume` flag to leverage cached results when appropriate
- Use user-confirmed input references that are readable by the intended execution environment

### File Handling
- When specifying multiple files, separate with comma and wrap in double quotes: `--input "file1.fq,file2.fq"`
- ALWAYS expand wildcards/globs to comma-separated file lists before running
- Reference task IDs from stdout to locate output files in work directories

### Module Include Syntax
- `include { MOD } from 'nf-core/module'` — **Nextflow managed module** (default). Nextflow automatically downloads and installs it from the Nextflow Registry. Always prefer this form. (`nf-core` here is the namespace; substitute the actual namespace of the module you are using.)
- For locally modified modules, use the actual module reference provided by the user's project rather than assuming its location.

### Module Selection
- Prefer an existing workflow module over hand-wiring the process modules it already chains together
- Find workflow modules via the registry API `kind=Workflow` filter (see Step 1); `nextflow module search` rarely surfaces them
- Do not write wrapper workflows to test single modules - use `run-module` instead
- Use `nextflow module search` to find modules, then `nextflow module view` for details

### Debugging Protocol
1. Check the task work directory using the task ID from stdout
2. Examine `.command.log`, `.command.err`, and `.command.out` files
3. Verify input files exist and are accessible
4. Check resource requirements (memory, CPUs) match available resources

## Skill Delegation (IMPORTANT)

**Read the relevant companion skill before starting its task:**

| Task | Invoke Skill |
|------|--------------|
| Install/run/test a module | `run-module` |

### When to Delegate

- **Step 3 (Validate Modules)**: Use `run-module` skill for each module validation

### Example: Step 3 Validation

For each module in your plan:
```
1. Invoke: run-module → install, run with test data, and verify outputs
2. Only proceed to next module after success
```

These skills contain detailed instructions for their specific tasks and ensure consistent execution patterns.

## Quick Reference

```
Step 1: Identify process + workflow modules → Propose plan to user
Step 2: Wait for user agreement
Step 3: Validate ALL modules ONE BY ONE with test data (untyped workflow modules: defer to Step 4)
Step 4: Compose final workflow (only after Step 3 succeeds)
```
