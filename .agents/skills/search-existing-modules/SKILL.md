---
name: search-existing-modules
description: >
  Given a bioinformatics tool name (e.g. MAFFT, Boltz, bcftools, salmon),
  search for existing Nextflow module implementations in priority order —
  nf-core/modules first, then a supplied source repository, then the
  broader community (GitHub, nf-core pipelines, biorxiv/arxiv
  supplementary pipelines). Return a verdict — reuse the nf-core module
  at path X, vendor/adapt the module from repo Y, or no existing
  implementation (in which case provide the canonical invocation for a
  new module). Use whenever a subworkflow is being scoped and needs to
  decide whether to reuse an existing Nextflow module or write a new one,
  or as a subagent task during the `nf-pipeline-design` /
  build-nextflow-pipeline flow. Self-contained — needs only the tool name and
  optionally the source repo path, no prior conversation context.
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


# Search for an Existing Nextflow Module

Given a bioinformatics tool and optional source material, find a reusable Nextflow implementation before proposing a new module. Do not assume a local repository or module directory.

## Search Order

1. Verify a supplied implementation hint when available. Read the referenced module source and confirm its inputs, outputs, and dependencies.
2. Search the Nextflow registry and nf-core module collection.
3. Inspect source material supplied by the user, including existing Nextflow code, container definitions, environment specifications, and helper scripts.
4. Search the broader community and tool authors' public examples.
5. Compare a candidate with the planned input/output contract before recommending reuse.

For native registry discovery when supported:

~~~bash
nextflow module search "quality control"
nextflow module view nf-core/fastqc
~~~

Use a candidate's current module template for direct validation on a small user-provided input. Do not create a wrapper solely to test a registry module.

## Evaluate the Candidate

Record the verified public source or user-supplied implementation reference, module identifier, process name, input tuple shape, output channels, required configuration, and declared container or Conda environment.

Check whether the module matches the planned upstream and downstream contracts, whether its dependency version suits the analysis, and whether the chosen execution environment can access its container and inputs.

For tools with multiple subcommands, inspect the relevant module family rather than assuming a single process covers every operation.

## User-Supplied Source Material

Look for reusable process definitions, a working container/environment, test fixtures, and helper scripts. Identify how these are referenced in the supplied material rather than prescribing folders.

A container and wrapper script can be valuable even when no Nextflow module exists. Record that partial match so subsequent module development can reuse verified behavior.

## Return a Verdict

### Reuse

Report the module identifier, verified implementation reference, process name, input and output contracts, required settings, dependency pinning, and validation result.

### Vendor and Adapt

Report the source, why it does not fit as-is, the minimum required adaptations, and reusable environment or helper components.

### Write New

Report the searches performed, the tool authors' canonical invocation, known input/output formats, container basis, and unresolved prerequisites. Do not claim no implementation exists after an incomplete search.

Keep the verdict concise and suitable for the caller's planning document. Distinguish a source inspection from an executed validation.

## Public References

- https://registry.nextflow.io
- https://nf-co.re/modules
- https://github.com/nf-core/modules
- https://nf-co.re/pipelines

Use create-container when a new module needs a verified execution environment.
