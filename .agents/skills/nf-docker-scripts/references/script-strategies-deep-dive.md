<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Script Strategy Deep Dive

Extended guidance for each script management strategy in Nextflow processes.

## Shared Helper Commands — Advanced Patterns

### Multi-language shared helper commands

Shared helper commands can use any language. The shebang determines the interpreter.

Keep the relevant process definition and helper command under the user's project
organization, and verify their availability in the selected execution environment.

The container image must include the matching interpreter. For multi-language pipelines,
a combined conda environment or multi-tool Wave container works well.

### Testing shared helper commands outside Nextflow

Shared helper commands are standalone executables and can be tested independently:

Test the selected helper independently with representative user-provided inputs.
Run its unit tests, type checks, and shell lint checks with the tools available in
the user's environment.

This is a major advantage over inline scripts or templates, which are harder to test in isolation.

### Shared helpers with nf-core modules

Verify how the selected project makes shared commands available to tasks. Module-owned
commands are available only when supported module packaging is configured; Wave module
binaries additionally require `nextflow.enable.moduleBinaries = true` and Wave.
For environments without Wave, use a project-supported command or template mechanism.

## template Directive — Advanced Patterns

### Template with meta map

Templates can access complex Nextflow objects:

```python
import json

sample_id = "${meta.id}"
is_paired = ${meta.single_end ? 'False' : 'True'}
reads = "${reads}".split()

print(f"Processing {sample_id}, paired={is_paired}")
print(f"Input files: {reads}")
```

### Template escaping rules

In template files, Nextflow interpolates `${...}` before the script runs:

- `${variable}` — resolved by Nextflow
- `\${variable}` — escaped, passed literally to the script language
- `$variable` (no braces) — resolved by Nextflow only in Bash templates; in Python/R, it's literal

For Python f-strings in templates, avoid `${...}` syntax for Python variables — Nextflow will
try to resolve them. Use regular Python variables and f-string `{var}` (no dollar sign) instead:

```python
name = "${meta.id}"  # Nextflow resolves this to the actual sample ID
count = len(open("${input_file}").readlines())  # Nextflow resolves ${input_file}
result = f"Sample: {name}, count: {count}"  # Pure Python f-string — no dollar signs
```

If you must use a shell environment variable inside a Python template, escape it with `\$`:

```python
import os
home_dir = os.environ.get("HOME", "\${HOME}")  # \$ prevents Nextflow interpolation
```

### Template vs shared-command decision

For this Nextflow release, prefer templates for multi-line scripts that belong to one process.
Prefer shared pipeline helper commands when the script is a reusable CLI shared across processes or pipelines.
Do not recommend module-owned helper commands as the default alternative; it is Wave-only behavior and
should only be discussed when the user explicitly asks for it.

| Factor | Template | shared helpers |
|--------|----------|------|
| Best default use | Process-specific multi-line scripts | Shared/reusable standalone CLIs |
| Nextflow variable access | Direct `${var}` interpolation | Must pass as CLI arguments |
| Reuse across processes | One template per process | One script, many processes |
| Independent testing | Harder (needs NF var substitution) | Easy (standalone executable) |
| Version control visibility | Beside module/process definition | Shared command sources are easy to review |
| Language mixing | Single language per template | Each script is independent |

## Wave Module Binaries — Technical Details (Explicit Request Only)

### How it works internally

Only discuss this section when the user explicitly asks for module binaries or an existing
pipeline already uses them. Do not recommend this approach by default.

When `nextflow.enable.moduleBinaries = true` and Wave is active:

1. Nextflow discovers module-owned helper commands through the supported module packaging convention
2. Wave extends the container image to include those scripts at a known `$PATH` location
3. The scripts are available inside the container as if they were installed tools
4. Different modules can have different scripts without name collisions

### Module shared helpers vs pipeline shared helpers precedence

With module binaries enabled, both module-owned and shared pipeline helper commands can be on
`$PATH`. If there's a name collision, module-owned helper commands take precedence for that module's processes.

### Fallback when Wave is unavailable

If a user has module-owned helper commands but cannot use Wave:

1. Move the scripts to the shared pipeline helper commands
2. Or include them inline in the process `script:` block
3. Or use the `template` directive (preferred for process-specific multi-line scripts)

## Dockerfile Patterns — When Custom Images Are Necessary

### Minimal Dockerfile for a Nextflow process

Use a compatible Micromamba base image, install the packages from the user's
environment specification, clean build caches, and run tasks as the supported
nonroot user. Keep pipeline logic supplied by the workflow rather than the image.

### Multi-stage build to keep images small

Build dependencies in a separate build stage, then copy the required runtime
artifacts into the final image using destinations appropriate for the chosen
base image. Verify command discovery and interpreter availability in that image.

### What belongs in the image vs. what doesn't

**In the image (tools and runtimes):**
- Language runtimes (Python, R, Java)
- Package dependencies (pip, conda, CRAN)
- Compiled bioinformatics tools (samtools, bwa, GATK)
- System libraries required by tools

**NOT in the image (pipeline logic):**
- Custom analysis scripts
- Pipeline-specific wrappers
- Data transformation scripts
- Reporting scripts
- Any script that changes when pipeline logic changes

### Docker layer caching and Nextflow

Nextflow caches work based on task inputs + container hash. If the container image changes
(because you rebuilt it with updated scripts), all cached results are invalidated even if
the tool version is unchanged. Keeping scripts out of the image preserves cache validity
across script-only changes.
