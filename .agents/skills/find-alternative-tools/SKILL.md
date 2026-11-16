---
name: find-alternative-tools
description: Given a named bioinformatics analysis step (e.g. "multiple sequence alignment", "somatic variant calling", "protein structure prediction", "single-cell demultiplexing"), find the credible alternative tools that could implement it, with pros/cons, license notes, typical compute requirements, maintenance status, and a ranked recommendation. Output is a structured comparison to inform tool selection; research does not authorize implementation of additional branches. Use whenever a subworkflow is being scoped and needs multiple tool options, when asked to "list alternatives for X", "what could we use instead of Y", "what tools exist for Z", or as a subagent task inside the `build-nextflow-pipeline` flow at the alternative-tools step. Self-contained — needs only the analysis step and any constraints, no prior conversation context.
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


# Find alternative tools for an analysis step

You have been given the name of an analysis step. Your job is to produce a ranked list of credible tools that can implement it, enough for the caller to decide which alternatives to implement as parallel branches in a Nextflow subworkflow.

## Scope

Compare credible tools for the requested analysis and constraints. Recommend one default and explain when alternatives are justified. The caller decides which candidate to use; implementing additional branches requires explicit user approval. An alternatives report is research, not a pipeline implementation plan.

Err on the side of inclusion. A tool that is slightly dated but still published-in-papers might still belong on the list as a reproducibility option. Skip only tools that are:

- Abandoned and unmaintained for >5 years with no clear successor
- Proprietary in a way that blocks distribution (unavailable under any license the pipeline could ship)
- Architecturally impossible in a Nextflow pipeline (e.g. require interactive GUI, cannot be containerized)

## How to search

Work in this order:

1. **Recent review papers and benchmarks** — search biorxiv, arxiv, pubmed for "<step> benchmark" or "<step> comparison" from the last 3–5 years. These give you the current tool landscape with pros/cons already weighed.
2. **nf-core pipelines that already do this step** — look at https://nf-co.re/pipelines and see which tools the canonical pipelines wrap. This is a strong signal of community acceptance.
3. **GitHub topic and keyword search** — look for active repositories with recent commits, sensible star counts, and clear CLI interfaces.
4. **Tool authors' own recommendations** — READMEs often list contemporaries or alternatives honestly.

If you have web search available, use it — the landscape shifts every 6–12 months and cached knowledge goes stale fast.

## What to report per tool

For each tool, a short block:

- **Name and canonical URL** (repo or tool homepage)
- **One-sentence description**
- **License** — critical for pipeline distribution; flag anything non-OSI or non-commercial
- **Maintenance status** — last release date, commit activity, responsive issue tracker
- **Compute profile** — CPU-only / GPU-required / memory footprint / typical runtime per sample
- **Container availability** — official Docker image, bioconda package, or needs-build-from-source
- **nf-core module hint** — if you happen to notice an existing nf-core module for this tool during your research on a public module page or in the user's source, record it here. This is not the job of this skill (that is `search-existing-modules`), but free signal is free signal — the next phase can verify and skip redundant discovery.
- **Output format** — what downstream tools need to consume
- **Strengths** — 1–3 bullets, concrete not marketing
- **Weaknesses / caveats** — 1–3 bullets, honest; include known bugs, quality issues, or sharp edges

## Ranked output

After the per-tool blocks, produce a ranked recommendation with reasoning, structured as:

- **Default branch** — the tool most users should use most of the time. Usually the best balance of quality, speed, and maintenance.
- **Accuracy branch** — if there's a clear "slower but better" option, name it.
- **Speed/budget branch** — if there's a "faster, slightly lower quality" option for large-scale runs, name it.
- **Niche / specialized branches** — tools that shine for specific organisms, data types, or hardware.
- **Do not implement** — tools you evaluated and rejected, with a one-line reason for each. This list matters as much as the recommended set — it prevents the caller from re-evaluating the same alternatives later.

## Parameter naming hint

Suggest schema `enum` values only for branches explicitly approved for implementation. A research shortlist must not become the pipeline's advertised supported choices. Use short, lowercase, conventional tool names for approved choices.

## What to return

Return the ranked list as markdown. Keep it terse — the caller will read every block, so padding hurts. A complete answer for one analysis step is usually 300–800 words.
