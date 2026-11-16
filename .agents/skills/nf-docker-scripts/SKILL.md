---
name: nf-docker-scripts
description: >
  Guide users on how to manage scripts in Nextflow processes and Docker images.
  Covers shared helper commands, template directive, inline scripts, and why module
  binaries should not be recommended by default. Use when the user asks about
  Docker images for Nextflow, where to put scripts, how to handle custom code
  in processes, or builds a container that bakes scripts in.
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


# Nextflow Script & Docker Strategy

Help the user choose the right approach for managing scripts in Nextflow processes,
especially when building Docker images. The guiding principle: **keep scripts out of
Docker images whenever possible.**

## When to Use

Load this skill when the user:

- Asks how to handle scripts in a Nextflow process
- Builds or writes a Dockerfile that COPYs pipeline scripts into the image
- Asks about shared helper commands, `template`, or module binaries
- Needs a container for a process that runs custom code
- Is unsure where to put a script that a Nextflow process calls

## Decision Tree

Ask the user which scenario matches their situation. If they don't know, walk through this tree:

```
Is the script a simple one-liner or short inline command?
  YES → Inline script: block (Strategy 1)
  NO ↓

Is the script reused across multiple processes or pipelines?
  YES → shared helper commands (Strategy 2) — shared scripts on $PATH
  NO ↓

Is this a multi-line, process-specific script?
  YES → template directive (Strategy 3) — safe, portable default for this Nextflow release
  NO → shared helper commands (Strategy 2) — reusable/testable fallback
```

**Recommended default for process-specific custom scripts:** `template` directive (Strategy 3).
Templates are portable across execution backends and avoid depending on Wave-only
module binary behavior.

**Recommended default for shared/reusable scripts:** shared pipeline helper commands (Strategy 2).

**Do not recommend module binaries as an approach.** Only discuss them if the user explicitly
asks about an existing module-bin setup or specifically requests `nextflow.enable.moduleBinaries`.
See Strategy 4 for details.

## Strategy 1: Inline `script:` Block

For short, simple commands. The script lives directly in the process definition.

```groovy
process SORT_BAM {
    container 'community.wave.seqera.io/library/samtools:1.21--abc123'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sorted

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam ${bam}
    """
}
```

Use when:
- The command is a single tool invocation or a few lines of Bash
- No custom logic beyond argument wiring

Do NOT use when:
- The script is >20 lines
- It contains complex logic, loops, or error handling
- It mixes languages (Python inside Bash)

## Strategy 2: Shared Helper Commands

Provide shared helper commands through the packaging mechanism supported by the
user's Nextflow execution environment. Verify that tasks can invoke them by name
before relying on them; local command availability is insufficient for remote execution.

Keep the relevant process definition and helper command under the user's project
organization, and verify their availability in the selected execution environment.

```groovy
process RUN_ANALYSIS {
    container 'community.wave.seqera.io/library/python_pandas:3.11--abc123'

    input:
    path input_csv

    output:
    path "*.csv", emit: results

    script:
    """
    run_analysis.py --input ${input_csv} --output output.csv
    """
}
```

Rules:
- Declare the interpreter supported by the execution environment
- Make the helper executable using the host's supported mechanism
- Use CLI arguments (argparse, click, getopts) — not hardcoded paths
- The Docker image needs only the runtime (Python, R, etc.) — **not the script itself**

Advantages:
- Scripts stay in version control with the pipeline
- Portable across any container that has the runtime installed
- Easy to test independently outside Nextflow
- Works with all execution backends (local, Slurm, AWS Batch, etc.)

## Strategy 3: `template` Directive

Use Nextflow's `template` directive for multi-line, process-specific scripts. For the current
Nextflow release, templates are the safest portable option when the script belongs to one process:
they work without Wave, keep pipeline logic out of Docker images, and keep the script beside the
process definition. Templates also allow `${var}` substitution by Nextflow before execution.

Use the template-loading convention supported by the selected Nextflow version
and the user's existing process definition:

Keep the relevant process definition and helper command under the user's project
organization, and verify their availability in the selected execution environment.

```groovy
process ANALYSIS {
    container 'community.wave.seqera.io/library/python_pandas:3.11--abc123'

    input:
    path input_file
    val threshold

    output:
    path "output.csv", emit: results

    script:
    template 'analysis.py'
}
```

The template file uses Nextflow variable interpolation:

```python
import pandas as pd

df = pd.read_csv("${input_file}")
filtered = df[df["score"] > ${threshold}]
filtered.to_csv("output.csv", index=False)
```

Use when:
- The script is multi-line and tightly coupled to one specific process
- Process script needs Nextflow variable interpolation (`${input_file}`, `${task.cpus}`, etc.)
- You want a portable alternative to module binaries for module-local custom code
- Escaping `\$` everywhere in an inline block would be painful or error-prone

Avoid when:
- The script is reusable across processes or pipelines (use shared pipeline helper commands instead)
- The script is a short one-liner (use inline `script:` instead)
- The script needs substantial independent unit testing as a standalone CLI (use shared helper commands instead)

## Strategy 4: Wave Module Binaries (Do Not Recommend)

**Do NOT recommend module binaries as the solution.** For this Nextflow release, prefer templates
for process-specific scripts and shared pipeline helper commands for shared scripts. Module binaries reduce
pipeline portability because they require Wave — pipelines using module binaries cannot run on
bare Nextflow, HPC schedulers without Wave, or environments where Wave is unavailable.

Only use this strategy when the user says something like:
- "use module binaries"
- "I want scripts in the module's helper packaging"
- "enable moduleBinaries"
- "use `nextflow.enable.moduleBinaries`"

If the user does not specifically request module binaries, always prefer Strategies 1-3
instead. Do not suggest module binaries merely because Wave is available.

### How it works

When the pipeline runs with Seqera Wave and `nextflow.enable.moduleBinaries = true`, Nextflow
bundles each module's helper commands into the container at runtime. This means scripts travel
with the module definition, not in the Docker image and not in the shared pipeline helper commands.

Keep the relevant process definition and helper command under the user's project
organization, and verify their availability in the selected execution environment.

```groovy
// nextflow.config
nextflow.enable.moduleBinaries = true
wave.enabled = true
docker.enabled = true
```

```groovy
process MY_TOOL {
    container 'community.wave.seqera.io/library/python:3.11--abc123'

    input:
    path input_file

    output:
    path "output.txt"

    script:
    """
    my_tool.py --input ${input_file} > output.txt
    """
}
```

Advantages (when Wave is available):
- Scripts are self-contained with the module — no shared pipeline helper commands coordination
- Different modules can have different scripts without name collisions
- Clean module packaging for nf-core-style reusable modules

**Cannot use when:**
- Wave is not enabled — module binaries require Wave to inject scripts into the container at runtime
- Running on bare Nextflow without Seqera Wave
- Pipeline needs to be portable to non-Wave environments

If the user asks for module binaries but is not using Wave, explain the Wave requirement and
guide them to Strategy 2 (shared helper commands) instead.

## Anti-Pattern: Baking Scripts into Docker Images

**Do NOT COPY pipeline scripts into Docker images.** This is the most common mistake.

For example, copying a changing pipeline helper into a runtime image ties every
helper edit to an image rebuild. Supply the helper through pipeline packaging instead.

Problems:
- Every script change requires rebuilding and re-pushing the image
- Script version drifts from pipeline version in git
- Cannot use the same image across different pipeline versions
- Defeats Nextflow's caching — image hash changes even when the tool hasn't changed
- Harder to review: scripts hidden inside image layers

Docker images should contain **tools and runtimes only** (Python, R, samtools, etc.).
Pipeline scripts should come from inline `script:`, shared pipeline helper commands, or templates — not
from Docker image layers, and not from module binaries unless the user explicitly requests
that Wave-only behavior.

**One exception:** third-party tool images (biocontainers, rocker, etc.) that ship their own
executables are fine — those aren't your pipeline scripts.

## Dockerfile Best Practices for Nextflow

When you do need a custom Dockerfile (because no Wave/conda/biocontainer image exists):

```dockerfile
FROM python:3.11-slim

RUN pip install --no-cache-dir pandas==2.1 numpy==1.26 biopython==1.82

# Supply pipeline helpers through the workflow's supported packaging mechanism.
```

Prefer this priority order for container images:
1. **Wave + conda** — `wave.strategy = 'conda,container'` resolves packages automatically
2. **Existing community images** — biocontainers, rocker, community.wave.seqera.io
3. **Custom Dockerfile** — only when packages aren't in conda or existing images

For more on container provisioning, load the `create-container` skill.

## Quick Reference

| Strategy | Where scripts live | Requires Wave? | Best for |
|----------|-------------------|---------------|----------|
| Inline `script:` | Process definition | No | Short commands, tool invocations |
| shared helper commands | Shared pipeline command packaging | No | Shared/reusable scripts, standalone CLIs |
| `template` | Process-specific template loading | No | **Default for process-specific multi-line scripts** |
| Module binaries | Module-owned command packaging | **Yes** | **Do not recommend; discuss only on explicit request** |
| ❌ COPY in Dockerfile | Docker image layer | No | **Never for pipeline scripts** |

## Interaction Pattern

When the user asks about Docker + scripts, follow this flow:

1. **Ask**: "Are you building a Docker image, or do you need to run a script inside a Nextflow process?"
2. **Recommend**: Based on the decision tree, suggest the appropriate strategy (default
   to `template` for process-specific multi-line scripts; shared helper commands for shared/reusable scripts)
3. **If they chose wrong**: If you see a Dockerfile that COPYs pipeline scripts, explain the anti-pattern and suggest the correct alternative
4. **If they ask for module binaries**: Only then discuss Strategy 4, and confirm Wave is available
5. **Cross-reference**: If they need container help beyond script placement, load `create-container`

## Related Skills

- `create-container` — Wave/conda container provisioning
- [registry composition](../build-nextflow-pipeline/references/create-workflow/README.md) — composing workflows from modules
- [Python conversion](../build-nextflow-pipeline/references/convert-python-script/README.md) — converting Python scripts to Nextflow processes
- `nf-pipeline-design` — understanding pipeline layout including shared helper commands
