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


# Convert Jupyter Notebook to Nextflow

## When to Use

- User asks to convert a Jupyter notebook to Nextflow
- User provides an `.ipynb` file and wants a pipeline
- User asks about migrating notebook-based analysis to Nextflow

## Instructions

Identify meaningful executable tasks rather than making each cell a process. Extract shared variables, cell-order dependencies, random state, working-directory changes and interactive assumptions explicitly. Wrap a single coherent analysis when appropriate; split independent tasks only when their file/value contracts are clear.

### Plan First

0. **Audit** what the notebook points at — run
   `audit-conversion-readiness` to find hardcoded paths you cannot
   reach and helper modules the notebook imports but that were not provided
1. **Analyze** the notebook: understand the analysis steps, inputs, outputs, and dependencies between cells
2. **Outline** the Nextflow process structure: identify logical groupings of cells that form distinct pipeline steps
3. **Identify** container/environment requirements for each process

### Execute

For each logical step:

1. Extract the relevant code cells into a standalone script
2. Create a Nextflow process wrapping that logic
3. Define input/output channels connecting processes
4. Specify the container environment (`container` directive)

### Verify

After generating the pipeline:

1. Validate compilation using the user's selected pipeline and installed runtime
2. If there are syntax errors, fix them and re-run
3. Run `nextflow lint` to check for best practice issues
4. Run representative input through the extracted source and the converted pipeline; compare expected scientific outputs and report any unresolved notebook state

Run validation with an available Nextflow installation and host execution tools.
If either is unavailable, report validation as unverified and read
`install-nextflow` for installation guidance.

### Common Patterns

- **Data loading cells** → first process with file inputs
- **Preprocessing/transformation cells** → intermediate processes
- **Analysis/modeling cells** → compute processes
- **Visualization/output cells** → final processes with `publishDir`
- **Imports at the top** → distribute to each process that needs them
- **Shared variables** → pass via channels between processes
