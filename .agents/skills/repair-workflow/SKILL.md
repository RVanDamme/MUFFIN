---
name: repair-workflow
description: >
  Repair or debug an existing Nextflow workflow or pipeline. Use when the user
  asks to fix, debug, or improve an existing workflow, especially when choosing
  the right validation loop matters more than writing new workflow structure.
  Also use for pipeline lint, config validation and preview/compilation diagnostics.
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


# Repair Existing Nextflow Workflows

Diagnose or repair the user's existing workflow in place. Work with the supplied files, modules, tests and validation commands rather than assuming a repository layout. A diagnosis or lint request authorizes inspection and reporting, not automatic source edits.

## Choose the Validation Loop

For syntax/config/compilation failures, read [static diagnostics](references/nf-debug/README.md). For executed failures, use `debug-local-run` or `debug-seqera-failed-run`; for test expectations, use `nf-test`. Static diagnostics supplement representative execution and cannot prove output correctness.

Run the existing verifier or affected nf-test target early. If the project explicitly requires nf-core linting, use its required lint command and metadata. For channel-shape bugs, inspect representative tuples before changing operators.

Check current official Nextflow documentation before changing operator, channel, or output semantics:
https://docs.seqera.io/nextflow/

## Make the Smallest Correct Edit

Fix the intended workflow, module, or wiring. Do not introduce a shadow implementation, wrapper, or broad redesign to avoid repairing the actual surface.

Do not rewrite fixtures, expected outputs, or snapshots merely to make validation pass.

## Deterministic Lint Repairs

When failures concern configuration or metadata drift, inspect the actual messages and repair the relevant setting. Examples include stale nf-core metadata, missing custom configuration settings, missing profile inclusion, or an incorrect DAG output format.

Classify failures as configuration drift, metadata drift, missing validation setup, structural issues, or stylistic warnings. Rerun after a mapped repair. If the same failure persists twice, inspect the verifier expectations and produced artifacts before further editing.

## Reason About Channels Explicitly

1. Write representative input tuples.
2. Write expected output tuples.
3. Identify the join key.
4. Identify the authoritative stream.
5. Then choose or modify operators.

Optional channels may enrich primary rows without manufacturing new primary rows. Be cautious with remainder joins and null padding.

## Avoid Fragile File Handling

Reject patches that rely on unverified project-relative helper locations, manual copy-back into an assumed output folder, trivial wrappers, or duplicate local implementations of the intended module.

Use helpers already available in the user's environment and preserve their staging contract. Do not introduce broad tool-version bumps, module patch churn, or snapshot rewrites without evidence.

## Static and Runtime Validation

Use Nextflow's built-in lint command when supported by the installed version:

~~~bash
nextflow lint
~~~

When available, use the bundled nflint structural linter described by nf-pipeline-design. Let the host locate the helper; do not prescribe its installation or project location.

Static validation complements semantic review and representative runtime tests.

## Staging Collisions

When inputs may share the same basename, inspect the process staging contract. Use Nextflow's stageAs facility or typed-process staging equivalent to assign distinct staged names while preserving the downstream tool's expectations.

Choose names from the actual input roles and cardinality. Do not prescribe fixed staging subdirectories. Validate that every tool reference agrees with the staged names and that multiple input files remain distinguishable.

## Escalate After Two Unproductive Cycles

Reread the test output, inspect the actual produced artifacts, compare expected and observed tuple/file shapes, and state a concrete failure hypothesis before editing again.

## Completion Criteria

- The intended validation ran.
- The defect is fixed in the existing surface.
- Fixtures were preserved unless an intentional contract change justified updating them.
- No shadow implementation or copy-back workaround was introduced.
- The representative test passes, or any remaining limitation is reported precisely.

Use [registry composition](../build-nextflow-pipeline/references/create-workflow/README.md) for new workflows, run-module for a single registry module, `migrate-nextflow-code` for strict-syntax changes, and [nf-test failure repair](../nf-test/references/repair-nf-test/README.md) for test-specific failures. Consult repair-checklists.md when available.

## Static diagnostics playbook

When the request concerns lint, config validation, compilation or why a pipeline will not start, use static diagnostics and honor whether the user requested diagnosis or repair. Read [static diagnostics](references/nf-debug/README.md) before proceeding.
