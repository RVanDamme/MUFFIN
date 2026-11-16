---
name: nf-test
description: Testing Nextflow pipelines with nf-test. Use when setting up nf-test, writing pipeline tests, testing modules/subworkflows, using nft-utils assertions, or debugging test failures. Trigger phrases include "test nextflow", "nf-test", "test pipeline", "write test for workflow", "snapshot testing".
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


# Testing Nextflow with nf-test

Use nf-test for pipeline, module, subworkflow, and function tests without assuming a project layout.

## Prerequisites and Configuration

Check whether nf-test and the user's intended Nextflow/container runtime are available. Use nf-test version and the current official installation instructions when setup is needed. Do not assume an installation location.

~~~bash
nf-test version
~~~

For a new test setup, nf-test init can generate configuration. Adapt it to the user's existing project and chosen test/work locations. Preserve existing profiles, fixtures, and test organization.

Configure required nf-test plugins in nf-test.config:

~~~groovy
config {
    plugins {
        load "nft-utils@0.0.9"
        load "nft-vcf@1.0.7"
        load "nft-bam@0.6.1"
    }
}
~~~

These versions describe the bundled reference examples. Confirm compatibility with the installed nf-test version when authoring a real test.

## Select the Test Surface

| Surface | nf-test declaration | What to validate |
|---|---|---|
| Pipeline | nextflow_pipeline | End-to-end success and required outputs |
| Module/process | nextflow_process | Input contract, outputs, versions, and failure behavior |
| Subworkflow | nextflow_workflow | Channel wiring and combined outputs |
| Function | nextflow_function | Return values, types, and edge cases |

Use the source file and target name actually supplied by the user. Test filenames use the .nf.test extension; do not prescribe where they must be stored.

## Write a Meaningful Test

Arrange the smallest representative input, invoke the intended surface, and assert the behavior the user cares about. Preserve the module's metadata and tuple shape.

A then block can combine independent checks:

~~~groovy
then {
    assertAll(
        { assert process.success },
        { assert process.out.bam.size() == 1 },
        { assert process.out.bai.size() == 1 },
        { assert snapshot(process.out.versions).match("versions") }
    )
}
~~~

For workflow tests, use workflow.success and the actual output channels. For functions, check function.success and function.result. Expected failures should assert the intended error and exit behavior.

## File Assertions

Prefer format-aware nf-test plugins to manual parsing:

- VCF: nft-vcf, including variantCount, sampleCount, chromosomes, variantsMD5, and targeted variant queries.
- BAM/CRAM: nft-bam.
- FASTA/FASTQ: the applicable format plugin.
- General text, JSON, YAML, checksums, and snapshots: nft-utils.

For example:

~~~groovy
with(path(process.out.vcf[0][1]).vcf) {
    assert variantCount > 0
    assert sampleCount == expectedSampleCount
    assert chromosomes == expectedChromosomes
}
~~~

The expected-value variables stand for values established from the user's fixture and test contract.

## Snapshots

Snapshot stable, meaningful output channels. Name separate snapshots when a test checks several surfaces. Sort order-independent values and filter timestamps or other unstable metadata deliberately.

Update a snapshot only after reviewing why the output changed and confirming the change is intentional. Do not update fixtures or snapshots to conceal a regression.

~~~bash
nf-test test --verbose
nf-test test --update-snapshot
~~~

Before updating, select the affected test using the project's supported test-selection mechanism. Do not indiscriminately update the whole suite.

## Running and Debugging

~~~bash
nf-test test
nf-test test --profile docker
nf-test test --verbose
nf-test test --debug
~~~

Use only profiles available in the user's configuration. Narrow execution to the affected test or tag before broadening to the relevant suite.

Inspect the work artifacts and logs reported by nf-test using available host tools. For a failed task, inspect .command.sh, .command.log, and .command.err from the actual task reported by the test. Do not assume their locations.

Treat failures as evidence. Read [nf-test failure repair](references/repair-nf-test/README.md) to distinguish regression from an intentional contract change. Check missing inputs, tuple shape, profile availability, container access, resource requirements, and actual assertion differences before changing pipeline logic.

## Setup, Cleanup, and CI

Use nf-test setup for genuine prerequisites and cleanup for temporary artifacts created by the test. Respect the user's chosen locations and preserve unrelated files.

In CI, reproduce the required Nextflow, nf-test, plugin, container, and profile environment. Pre-download plugins before parallel tests when required by the installed plugin mechanism. Use small representative inputs and deterministic assertions.

## Bundled References

Consult plugins.md, nft-utils.md, assertions.md, test-data.md, and troubleshooting.md by name when the host makes them available.

For failing suites, assertion errors or snapshot mismatches, read [nf-test failure repair](references/repair-nf-test/README.md) before changing code or expectations. Preserve fixtures unless a reviewed intentional contract change justifies updating them.

## Failure repair playbook

When an existing nf-test fails, classify intentional contract change, regression, nondeterminism or infrastructure before editing code, assertions or snapshots. Read [nf-test failure repair](references/repair-nf-test/README.md) before proceeding.
