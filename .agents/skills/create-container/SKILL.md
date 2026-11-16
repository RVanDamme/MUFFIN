---
name: create-container
description: >
  Build, claim, or recommend a verified container for a bioinformatics tool —
  and prove it runs the intended command on representative test data before
  handing it off for module authoring. Make sure to use this skill whenever
  the user wants to containerize a tool, wrap a tool for a Nextflow module,
  set up Wave for a custom package stack, pin a Bioconda/conda-forge
  environment, or unblock a "this tool's dependencies are weird" problem.
  Trigger phrases include "containerize X", "I need an image for Y",
  "set up Wave for Z", "build a conda env for this tool", "verify this
  container works on test input", "what container does <tool> need",
  "wrap <tool> for a Nextflow module", "fix this conda solve". Follows a
  build → run-test-input → iterate loop and returns a structured handoff
  (verified image reference + working command). Do not author the Nextflow
  module here — that belongs in `build-nextflow-pipeline` after verification.
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


# Container Provisioner

Provision or recommend containers for conda-based bioinformatics tools, preferably with Seqera Wave.

## When to Use

Use this skill when the user:

- needs a container for a tool or package set
- wants a reproducible container without hand-writing a Dockerfile
- needs a container for validating a module or workflow step

## Build → run test input → iterate loop

This skill is not just a package-name recommender. For a new module or tool
wrapper, the useful output is a container **and** a command that has already run
successfully on representative data.

1. **Resolve the package/image spec.** Prefer pinned Bioconda/conda-forge
   packages, an official image, or a minimal Dockerfile when Conda cannot model
   the environment.
2. **Build or claim the container.** Use Wave when available; otherwise provide
   the exact conda YAML, Dockerfile, or image reference.
3. **Run the real command inside the container** against the supplied test input
   using Docker/Singularity or the available runtime.
4. **Check real outputs.** Success means exit code 0, expected files exist and
   are non-empty, output format looks right, and stderr has no unexplained red
   flags.
5. **Iterate** on dependencies, versions, paths, and flags until the command is
   genuinely runnable. Only then hand the image/command to module authoring.

Do not write a Nextflow module here. Environment bugs are faster to catch at the
shell/container layer than through a Nextflow process wrapper.

## Test input requirements

Ask for or derive a small representative input before declaring the container
usable:

- tool name and version
- input path and format
- expected output shape
- command the eventual module should run
- any reference data/model/GPU requirement

If any of these are missing, return the gap instead of claiming the image is
ready.

## Step 1: Resolve the Package Spec

First determine the exact package names and versions the tool needs.

Preferred approach:

- use the available Seqera package-search capabilities if they exist
- otherwise identify the package name and a concrete version from trusted package metadata before recommending a container

Use channels in this order unless the user says otherwise:

- `bioconda`
- `conda-forge`

## Step 2: Claim or Describe the Container

If Wave container-claim capabilities are available in the current environment:

- request a container from the resolved package set
- return the resulting image reference

If not:

- provide the exact package entries needed for Wave/conda resolution
- show how to use them from Nextflow config

Example package entries:

```text
bwa=0.7.19
samtools=1.17
multiqc=1.19
```

## Nextflow Usage

Recommended Wave configuration:

```groovy
wave.enabled = true
wave.strategy = 'conda,container'
docker.enabled = true
```

Fixed image usage:

```groovy
process {
    container = 'community.wave.seqera.io/library/bwa:0.7.19--example'
}
```

## Frozen vs Ephemeral Containers

- use ephemeral containers for quick testing
- use frozen or permanent image references for reproducible production workflows

If the environment supports a `freeze`-style option, use it when the user needs:

- stable image names
- long-lived reproducibility
- published workflow references


## Return format

Return a structured block that can feed the module-writing step:

```text
TOOL: <name/version>
IMAGE: <verified image ref or exact conda/Dockerfile spec>
COMMAND (working): <literal command tested in the container>
INPUT FORMAT: <what was tested>
OUTPUT FORMAT: <what was produced>
VERIFIED ON: <test input path/description>
NOTES: <warnings, GPU/platform/reference-data caveats>
FAILED ATTEMPTS: <brief list so the next maintainer does not repeat them>
```

## Practical Rules

1. Resolve exact package names before asking for a container
2. Prefer reproducible pinned versions
3. Use Wave/conda before custom Dockerfiles when possible
4. When multiple tools belong together, combine them into one container only if they are actually used together
5. If the environment lacks live Wave APIs, still provide the concrete package spec and Nextflow config
6. Do not call a container verified until the command has run on test data or you clearly label it as unverified
