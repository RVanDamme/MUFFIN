<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../../../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Convert Python Script to Nextflow

Convert standalone Python scripts into idiomatic Nextflow DSL2 processes and workflows.

## When to Use

Load this skill when the user wants to:
- Convert a Python script to a Nextflow process
- Migrate a Python-based analysis pipeline to Nextflow
- Wrap Python tools in Nextflow processes
- Port Python data processing logic to Nextflow DSL2

## Workflow

0. **Audit readiness** — `audit-conversion-readiness` for anything
   beyond a self-contained script: it finds hardcoded paths you cannot reach and
   imports or helper scripts that were never provided, before you design around
   them.
1. **Analyze the script** — Identify inputs, outputs, dependencies, and logic.
2. **Determine conversion strategy** — Wrap vs. rewrite (see below)
3. **Generate Nextflow process(es)** — One process per logical unit
4. **Generate workflow** — Wire processes with channels
5. **Add `publishDir`** — Expose final outputs to user-accessible directories
6. **Create conda/container environment** — Capture Python dependencies
7. **Write tests** — Generate nf-test stubs for each process

## Conversion Strategy

### Strategy 1: Wrap (Preferred)

Keep the Python script intact and wrap it in a Nextflow process. Best when:
- Script is complex or well-tested
- Uses libraries heavily (pandas, numpy, scikit-learn, etc.)
- Script is maintained externally

```groovy
process PYTHON_ANALYSIS {
    publishDir params.outdir, mode: 'copy'

    conda "python=3.11 pandas=2.1"
    container "biocontainers/python:3.11"

    input:
    path input_file
    path analysis_script

    output:
    path "output.csv", emit: results

    script:
    """
    python "${analysis_script}" \\
        --input "${input_file}" \\
        --output output.csv
    """
}
```

**Wrapping rules:**
- Make the user's script available as a staged input or a verified executable in
  the selected environment; do not assume a helper location.
- Derive arguments and output names from the script's actual interface.
- Use `argparse` or `click` for CLI args (refactor if using hardcoded paths)
- Pin dependencies in `environment.yml` or `requirements.txt`

### Strategy 2: Rewrite

Rewrite Python logic as Nextflow native code or simpler scripts. Best when:
- Script is short (<50 lines of logic)
- Logic is simple file manipulation, parsing, or reformatting
- No heavy library dependencies

```groovy
process PARSE_RESULTS {
    input:
    path csv_file

    output:
    path "filtered.csv", emit: filtered

    script:
    """
    python - <<'PYTHON'
    import csv
    with open("${csv_file}") as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader if float(r["score"]) > 0.5]
    with open("filtered.csv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    PYTHON
    """
}
```

## Input/Output Mapping

Map Python I/O patterns to Nextflow channels:

| Python Pattern | Nextflow Equivalent |
|---------------|---------------------|
| `open("input.txt")` | `path input_file` (input) |
| `sys.argv[1]` | `val arg` or `path file` (input) |
| `argparse` positional | `val`/`path`/`tuple` inputs |
| `argparse` flags | `val` inputs |
| `glob.glob("*.fastq")` | Channel.fromPath("*.fastq") |
| Listing the user's input selection | A channel constructed from the user-provided input parameter |
| `open("out.csv", "w")` | `path "out.csv"` (output) |
| Creating an output destination | Declare the actual task outputs and publish to the user-selected destination |
| `pd.read_csv()` | `path csv_file` (input) |
| `plt.savefig("plot.png")` | `path "plot.png"` (output) |
| Return value / print | `env` or `stdout` (output) |

## Dependency Management

### Extract Dependencies

Analyze imports and generate environment spec:

```yaml
# environment.yml
name: analysis
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - pandas=2.1
  - numpy=1.26
  - biopython=1.82
```

### Container Strategy

```groovy
process PYTHON_TOOL {
    conda "python=3.11 pandas=2.1"
    container "${ workflow.containerEngine == 'singularity' ?
        'oras://community.wave.seqera.io/library/pandas_numpy:2.1--abc123' :
        'community.wave.seqera.io/library/pandas_numpy:2.1--abc123' }"
    // ...
}
```

- Prefer Wave containers via `community.wave.seqera.io` for bioinformatics stacks
- Use BioContainers when available for standard tools
- Fall back to custom Dockerfile only when necessary

## Process Design Rules

1. **One logical task per process** — Don't combine unrelated steps
2. **Declare all outputs explicitly** — Use `path`, not shell redirects to undeclared files
3. **Use `task.cpus` and `task.memory`** — Don't hardcode resource values
4. **Stage inputs, don't copy** — Nextflow handles file staging; use `path` inputs
5. **No hardcoded paths** — Use `${projectDir}`, `${launchDir}`, or inputs
6. **Idempotent scripts** — Process must produce same output given same input
7. **Exit codes matter** — Ensure Python script returns non-zero on failure (`sys.exit(1)`)

## Workflow Wiring

Convert Python's sequential flow to Nextflow channels:

```python
# Python sequential
data = load_data("input.csv")
cleaned = clean(data)
results = analyze(cleaned)
plot(results)
```

```groovy
// Nextflow DSL2
workflow {
    input_ch = Channel.fromPath(params.input)

    LOAD_DATA(input_ch)
    CLEAN(LOAD_DATA.out.data)
    ANALYZE(CLEAN.out.cleaned)
    PLOT(ANALYZE.out.results)
}
```

### Channel Patterns

| Python Pattern | Nextflow Channel Pattern |
|---------------|-------------------------|
| `for f in files:` | `Channel.fromPath(...)` parallelizes automatically |
| `if condition:` | `.branch { }` or `.filter { }` |
| `zip(a, b)` | `.join()` or `.combine()` |
| `dict[key]` | Tuple channels with key: `[meta, file]` |
| `try/except` | `errorStrategy 'retry'` or `'ignore'` |
| Sequential loop | Chained processes (automatic parallelism) |
| `multiprocessing` | Preserve intentional within-task parallelism; let Nextflow schedule independent task instances |

## Publishing Outputs

Nextflow runs processes in isolated work directories. Use `publishDir` to copy/link results to user-visible locations:

```groovy
process ANALYZE {
    publishDir params.outdir, mode: 'copy'
    // ...
}
```

- `mode: 'copy'` — copies files (safe, portable)
- `mode: 'symlink'` — symlinks (fast, saves space, breaks if work dir cleaned)
- `mode: 'link'` — hard links (fast, same filesystem only)
- Use `pattern:` to selectively publish: `publishDir "${params.outdir}", pattern: "*.csv"`

Map Python's `shutil.copy(result, output_dir)` or `os.rename()` calls → `publishDir` instead.

## Snakemake Migration

For full Snakemake workflow conversion, load the dedicated `migrate-from-snakemake` skill.
This skill focuses on standalone Python script conversion.

## Common Pitfalls

1. **Separate task scheduling from tool threading** — Preserve intentional within-task parallelism and expose resource requirements via `task.cpus`. Nextflow schedules independent instances; it does not replace a tool's internal parallel algorithm.
2. **Don't read/write to absolute paths** — Nextflow work dirs are isolated. Use relative paths in scripts.
3. **Don't manage temp files** — Nextflow handles work directory cleanup.
4. **Don't install dependencies in script block** — Use `conda`/`container` directives.
5. **Watch for statefulness** — Python scripts that write to global state or databases need special handling (use `val` outputs or external storage).

## Nextflow Script Block Details

Nextflow `script:` blocks normally use a shell. Invoke the chosen interpreter
explicitly when embedding another language:

```groovy
process INLINE_PYTHON {
    script:
    """
    python - <<'PYTHON'
    print("Hello from Python")
    PYTHON
    """
}
```

**Variable interpolation**: Nextflow interpolates `${var}` before the script runs. Escape Python f-strings and `$` usage:
- Use `\$` for shell/Python variables: `\${PYTHONPATH}`, `\$(which python)`
- Nextflow variables resolve first: `${input_file}` → actual path

**Multi-language scripts**: Keep shell and Python execution explicit. Use a
verified interpreter invocation, separate processes, or a user-provided helper
made available through staged inputs or the selected environment.

## Params Mapping

Map Python CLI arguments to Nextflow `params`:

```python
# Python argparse
parser.add_argument("--threshold", type=float, default=0.5)
parser.add_argument("--output-dir", required=True)
parser.add_argument("input_file")
```

```groovy
// nextflow.config
params {
    input     = null      // required — no default
    threshold = 0.5
    outdir    = null
}
```

```groovy
// workflow
workflow {
    Channel.fromPath(params.input, checkIfExists: true)
    | PROCESS_DATA
}
```

## Checklist

Before delivering the converted pipeline:
- [ ] All Python inputs mapped to Nextflow inputs
- [ ] All Python outputs declared as Nextflow outputs
- [ ] `publishDir` configured for user-facing outputs
- [ ] Dependencies captured in `environment.yml` or container
- [ ] Wrapped Python script is available as an input or verified executable
- [ ] No hardcoded paths remain
- [ ] `task.cpus`/`task.memory` used for resource params
- [ ] `params` defined in `nextflow.config` for user-configurable values
- [ ] Workflow wires processes correctly via channels
- [ ] Error handling preserved (non-zero exit on failure)
- [ ] String interpolation escaped correctly (`\$` for shell/Python vars)
- [ ] Representative execution matches the source's expected scientific outputs, using justified numerical tolerances when needed
