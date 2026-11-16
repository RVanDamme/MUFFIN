---
name: maintain-nf-core-pipeline
description: >
  Maintain an existing nf-core pipeline: work on the `dev` branch, resolve
  outstanding template-sync PRs, update shared modules from upstream, then
  patch any remaining code. Use whenever the user asks to update, maintain,
  sync, modernize, or contribute to an nf-core pipeline; whenever they
  mention resolving a template sync, a `tools X.Y.Z` merge, or stacked
  template PRs; or when they want to bring an nf-core pipeline up to a
  current tools version. Trigger phrases include "sync the template",
  "resolve the template merge", "update nf-core modules", "modernize
  atacseq", "fix this nf-core pipeline", "rebase against tools 4".
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


# Maintain an nf-core pipeline

Order of operations to bring an existing nf-core pipeline current. Do
these *in order* — most fixes higher up vanish once an earlier step
lands, so skipping ahead usually creates redundant code changes that
later conflict with a template sync or a module update.

## 1. Work on `dev`, not `main`/`master`

The default branch on every nf-core repo is the latest release, which
lags active code by months. Check out `dev` and target any PR against
`dev`. Maintenance work on `main` will be invisible to the next release
and will conflict at sync time.

For the full set of nf-core conventions (branch model, where edits
go, CLI commands, common "don't do" rules), read
`nf-core-conventions.md` once at the start of any
unfamiliar pipeline.

## 2. Resolve the outstanding template-sync PR

nf-core auto-opens a template-sync PR each time `nf-core/tools` cuts a
release. They stack up. Most maintenance fixes already exist in the
template — landing the oldest unresolved sync PR usually fixes them
for free, so this is the cheapest first move.

1. List open PRs labeled `template`, oldest first.
2. Take the oldest. Read the tools version from the PR title
   (e.g. *"Template update for nf-core/tools v4.0.0"*), then:
   - If a per-release guide exists at
     the bundled template merge guides, read it before touching conflicts.
   - Otherwise fetch the matching blog post at
     `https://nf-co.re/blog` and resolve conflicts from there.
3. Land the sync on `dev` before doing anything else.

## 3. Update shared modules

```bash
nf-core modules update --all
```

Many upstream module fixes land for free here. **Never edit a shared
module in-place** — those files are owned upstream and edits will be
overwritten by the next `update`.

## 4. Fix shared modules upstream, not in the pipeline

When a real bug lives inside a shared module:

1. Make the edit in the `nf-core/modules` repo.
2. Open a PR, wait for review and merge.
3. Back in the pipeline: `nf-core modules update <name>`.

If the fix is urgent and upstream review is slow, vendor the module to
a user-selected project component with a `// vendored from nf-core/modules#<PR>`
note and a TODO to delete the vendored copy after upstream merges.
Don't fork silently.

## 5. Patch remaining pipeline code

After 1-4, address whatever lint or test failures are left in
the entry workflow and project-owned workflow definitions. By this point the
residue is genuinely pipeline-specific.

If the residual work is a syntax migration:
- `migrate-nextflow-code` — strict syntax and other language migrations
- [schema migration](../nextflow-schema/references/nf-schema-migration/README.md) — schema updates required by v2

For structural analysis of an unfamiliar pipeline before editing:
- `nf-pipeline-design`

## References

- `nf-core-conventions.md` — branch model, file ownership, CLI commands, "don't do" rules.
- Template merge guides — per-release merge-conflict guidance.
- All `nf-core/tools` release blog posts: https://nf-co.re/blog
- Pipeline sync tutorial: https://nf-co.re/docs/tutorials/sync_a_pipeline
