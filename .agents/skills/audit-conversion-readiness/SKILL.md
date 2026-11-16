---
name: audit-conversion-readiness
description: >
  Read-only first-pass audit of source material (scripts, notebooks,
  Snakemake/WDL/CWL workflows, legacy in-house pipelines) before converting or
  migrating it to Nextflow. Inventories every hardcoded or external file path,
  probes which of them are actually reachable from this session, and flags tool
  code and dependencies that are referenced but were never provided. Run this
  before any conversion work, and whenever the user asks to "audit this before
  converting", "find the hardcoded paths", "can you actually access these files
  or reference data", "what's missing before we convert this", "what will block
  this port to Nextflow", or "why can't you read this path". Output is a
  blockers-first report plus one consolidated list of what the user needs to
  supply — no pipeline code.
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


# Audit conversion readiness

A conversion fails slowly when the source material points at data and code that
is not actually available. The expensive failure mode is discovering this halfway
through: a pipeline gets written against paths nobody can read, or against a
helper script that was never in the handover.

This skill is the cheap first pass. It is read-only: inventory, probe, report.
Write no `.nf` files and edit nothing during the audit. What happens next
depends on the verdict:

- **BLOCKED** — do not convert; the blockers would make the pipeline guesswork.
- **READY WITH GAPS** — conversion can start; named items must land before an
  end-to-end run.
- **READY** — nothing blocking; hand off to the conversion skill.

## Step 0 — establish what access you actually have

Every reachability claim you make later depends on this, so settle it first and
name it in the report.

| Session type | Probe with | Can reach |
|---|---|---|
| Local execution (CLI) | the host's local shell tool, the host's local file-reading tool, a directory listing through the host shell | The user's machine and whatever it mounts |
| Cloud sandbox (web) | the host's sandbox shell tool, if available, the host's sandbox file tool, if available | Only the sandbox VM, cloned repos, and downloads |
| Neither | — | Only what was pasted into the conversation |

Two consequences worth stating out loud:

- Access to a remote execution environment does not imply access to the user's
  local files, lab storage, or cluster mounts. Verify the actual access available.
- Cloud URIs (`s3://`, `gs://`, `az://`) may be reachable, but only with
  credentials. Absence of credentials is not absence of data.

Never report a path as verified unless you ran a probe that returned. When you
could not probe, the status is `unknown` — that is a useful answer.

## Step 1 — inventory every external reference

Collect references exhaustively before judging any of them. Read
the bundled reference-patterns guidance, when available, for per-language
signals covering user-provided locations, cloud references, URLs, working-directory
changes, configuration loads, and environment-derived references.

Record for each hit: the literal reference, `file:line`, and what it is for
(input data, reference/index bundle, output location, config, helper code,
scratch). Provenance matters — the user needs to audit your calls, and the path's
role decides its fate in the pipeline.

Exclude irrelevant generated task artifacts, dependency environments, vendored
material, and notebook output cells from the report.

## Step 2 — probe reachability

For each path, use the tools from Step 0 to establish existence, readability, and
rough size or file count. Batch the probes into one command per surface rather
than one call per path.

Assign exactly one status per reference:

| Status | Meaning |
|---|---|
| `reachable` | Probed and readable from here |
| `missing` | Probed here, and it does not exist |
| `unreachable-from-here` | Real for the user, structurally out of reach for this session |
| `needs-credentials` | Cloud URI with no working credentials in this session |
| `unknown` | Not probed — say why |

For cloud URIs, check whether the workspace already has a data link covering the
bucket before asking the user for keys. `search_seqera_api` only returns API
names — you must then `call_seqera_api` to list them (the authenticated MCP connection handles credentials; do not curl Platform). Typical sequence:

1. `search_seqera_api(query="list data links in the workspace")`
2. `call_seqera_api(service="platform", api_name=<the list-data-links name from
   step 1>, parameters={...})` with the workspace id identified by the user or
   returned by an authenticated workspace lookup
3. Match each `s3://` / `gs://` / `az://` URI against the returned links

Load `seqera-data-links` only if you then need to *create* a link.
An existing data link is usually the answer to "how will the pipeline read this".

## Step 3 — audit tool and helper code availability

The second silent blocker: code the source calls but the handover omitted. Every
`import`, `source()`, `library()`, `script:`, wrapper reference, jar, and bare
shell command is a dependency claim to resolve against what you were actually
given.

Read `tool-availability.md` for the resolution procedure and the checks that
matter: bioconda/biocontainers/Wave availability, existing nf-core
modules (`search_nfcore_module`), pinned versions, and the categories that block
a conversion outright — license-gated tools, registration-gated downloads,
GPU-only tools, and tools that need a separate reference bundle.

Resolve unpinned versions into an explicit question rather than a guess. "Which
`samtools` version was this validated against?" is cheap now and expensive after
the modules exist.

## Step 4 — report, blockers first

Read `report-template.md` for the exact shape. In summary:

1. **Verdict** — ready to convert / ready with gaps / blocked, and the surface
   you audited from.
2. **Blockers** — the short list that stops work, each with the one thing that
   would clear it.
3. **Full inventory table** — every reference, its role, status, and its intended
   fate in the pipeline (parameter, samplesheet column, container, staged input).
4. **Consolidated asks** — every question in one numbered list.

Ask once. A drip-feed of one question per turn is the exact time sink this audit
exists to prevent.

## Step 5 — hand off

With the audit answered, route onward:

- [parameter triage](../nextflow-schema/references/triage-pipeline-parameters/README.md) — turn the surviving paths and
  constants into a parameter surface and schema
- `create-container` or `nf-docker-scripts` — tools
  with no usable image
- `search-existing-modules` — tools that may already have a module
- `find-alternative-tools` — tools blocked by license, GPU
  requirements, or abandonment
- `seqera-data-links` — cloud data the pipeline must read
- Then the conversion skill itself: [Python conversion](../build-nextflow-pipeline/references/convert-python-script/README.md), [R conversion](../build-nextflow-pipeline/references/convert-r-script/README.md),
  [notebook conversion](../build-nextflow-pipeline/references/convert-jupyter-notebook/README.md), `migrate-from-snakemake`, or [registry composition](../build-nextflow-pipeline/references/create-workflow/README.md)

## Guardrails

- Report what you probed, not what you assume. An unprobed path is `unknown`.
- Never substitute a plausible-looking reference file, genome build, or container
  for one you could not find. A wrong reference genome produces a pipeline that
  runs and is silently incorrect.
- Never treat a missing file as a reason to drop the step from the conversion —
  it is a question for the user.
- Hardcoded paths are findings, not bugs to fix in place. The pipeline's
  parameter surface is where they get resolved, in [parameter triage](../nextflow-schema/references/triage-pipeline-parameters/README.md).
- Keep the audit read-only. No edits to the source material, no scaffolding.
