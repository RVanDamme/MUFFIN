<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Audit Report Shape

Lead with the verdict and blockers, then provide the inventory needed to audit
your findings. Use references taken from the user's material, never invented
locations or a prescribed directory layout.

## Verdict

State `READY`, `READY WITH GAPS`, or `BLOCKED`, the number and nature of
blockers, the source material examined, and the execution access available.

`READY` means no observed blockers. `READY WITH GAPS` means conversion can
start but named items are still required before end-to-end execution.
`BLOCKED` means proceeding would require guessing important dependencies.

## Blockers

For each blocker, identify the user-provided reference and source filename,
line, or notebook cell. Give one concrete clearing action:

- Confirm the reference genome build and index-producing tool version.
- Supply missing helper code or approve a defined replacement.
- Explain how a required license or restricted download is made available.
- Confirm the intended input or execution environment.

Do not substitute a plausible reference or drop a missing scientific step.

## Inventory

Use a table with columns for the user's reference, role, provenance, status,
and intended conversion treatment. Sort blockers first, then unknown items,
then resolved references.

Distinguish `reachable`, `missing`, `unreachable-from-here`,
`needs-credentials`, and `unknown`. Report probe results separately from
assumptions. Include a dependency table from the tool-availability guidance
when external tools or helpers are involved.

## Consolidated questions

Collect all unresolved questions into one numbered list. Ask for the minimum
evidence needed to clear each item: data access, reference provenance, missing
code, validated tool versions, license availability, or permission to replace
a specific tool.

Report status without blame. Summarize irrelevant generated artifacts instead
of padding the report. A clean audit should conclude briefly and hand off to
conversion.
