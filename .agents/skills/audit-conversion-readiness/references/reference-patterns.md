<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Finding External References in Source Material

Inventory every input, output, configuration, helper, reference bundle, and
temporary-storage dependency touched by the user's source. Record its literal
user-provided reference, source filename and line or notebook cell, and role.
Do not prescribe where those items must live.

## General sweep

Use the host's available search and file-reading tools to inspect the source
material the user supplied. Include literal locations, assembled references,
environment-variable roots, cloud references, and external URLs. Exclude
generated task artifacts, dependencies, notebook output cells, and vendored
material when they are irrelevant to the source contract.

A literal-string search is incomplete. Follow expressions that assemble names
or locations, and resolve configuration values before claiming the inventory is
complete.

## Language-specific signals

- Python: `open`, `Path`, `os.path.join`, `os.chdir`, dataframe readers and
  writers, environment lookups, import-search changes, globs, and subprocesses.
  Distinguish a fixed root from a variable component.
- R: `setwd`, `getwd`, table readers and writers, RDS operations, `file.path`,
  `source`, `library`, `require`, environment lookups, and `normalizePath`.
  A working-directory change can affect every subsequent relative reference.
- Shell and Makefiles: working-directory changes, sourced scripts, environment
  assignments, module loads, and exported variables. Module loads often provide
  tool names and exact versions.
- Snakemake: configuration keys, included configuration, helper scripts, wrappers,
  Conda specifications, containers, shell blocks, and input/output declarations.
  Missing configuration or helper code is a dependency gap.
- Nextflow: `params`, `publishDir`, `storeDir`, `workDir`, containers,
  `includeConfig`, and project/launch/base-directory expressions.
- WDL and CWL: file declarations, command sections, containers, and input
  locations.
- Notebooks: inspect code cells rather than stale output. Report cell indices
  where practical and use only user-approved locations if temporary extraction
  is necessary.
- Configuration: inspect configuration and lockfiles the source actually loads.
  Missing configuration can block enumeration of the underlying references.

## Classify each reference

| Role | Typical conversion treatment |
|---|---|
| Sample or input data | Sample-sheet field or explicit pipeline input |
| Reference or index bundle | User-selected pipeline parameter or data link |
| Output destination | User-selected output parameter and publishing behavior |
| Configuration | Parameter defaults or execution profile |
| Helper code | Reproducibly available executable or module dependency |
| Scratch or temporary storage | Runtime-managed task storage |
| Logs | Runtime-managed execution evidence |

Temporary storage and logs may disappear as source requirements during
conversion. Explain that treatment instead of guessing replacement locations.
