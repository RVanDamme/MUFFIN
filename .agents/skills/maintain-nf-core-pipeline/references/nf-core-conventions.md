<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nf-core Pipeline Conventions

Quick-reference for the conventions that apply across all nf-core
pipelines. Source: nf-co.re/docs and `nf-core/tools` defaults.

## Branch model

Three required branches on every nf-core pipeline:

| Branch | Purpose |
|---|---|
| `dev` | Active development. **All maintenance PRs target this.** |
| `master` / `main` | Stable releases. Lags `dev` by months. |
| `TEMPLATE` | Holds the synced nf-core template state. Updated by `nf-core pipelines sync`. Never edit by hand. |

Releases are cut from `dev` → `master` via a release PR. The default
branch shown on GitHub is `master`/`main` (the latest release), which
is why `git clone` lands you behind the active code unless you check
out `dev` immediately.

## Where edits go

| Path | Owner | Edit policy |
|---|---|---|
| Project-owned workflow and process definitions | Pipeline | Edit freely. |
| Imported nf-core components | nf-core/modules (upstream) | **Never edit in place.** Fix upstream and pull update. |
| `nextflow.config`, `nextflow_schema.json` | Pipeline (with template defaults) | Edit pipeline-specific blocks; preserve template-managed blocks during sync. |
| `.nf-core.yml`, automation configuration, `.pre-commit-config.yaml`, template documentation, `ro-crate-metadata.json` | Template | Synced from `nf-core/tools`. Resolve conflicts during template sync; don't pre-emptively edit. |
| `multiqc_config.yml` | Mixed | Pipeline-specific content blocks are pipeline-owned; structural defaults are template-owned. |

The litmus test: if a file is regenerated/updated by `nf-core pipelines
sync` or `nf-core modules update`, treat it as upstream-owned.

## CLI commands for maintenance

```bash
# Sync the latest template state from nf-core/tools into TEMPLATE branch,
# then auto-open a PR into dev with the changes.
nf-core pipelines sync

# Update every shared module to its latest upstream version.
nf-core modules update --all

# Update a specific shared module.
nf-core modules update <module-name>

# Install a new shared module from nf-core/modules.
nf-core modules install <module-name>

# Patch a shared module locally (creates a tracked diff against upstream).
# Use this when you need a temporary fix and can't wait for upstream PR.
nf-core modules patch <module-name>

# Run lint against the pipeline.
nf-core pipelines lint
```

## Common "don't do" rules

- Don't edit anything under imported nf-core components directly. Fix upstream and pull the update.
- Don't merge to `main`/`master` directly. Maintenance PRs go to `dev`; releases roll from `dev` to `master` via a release PR.
- Don't ignore the template-sync PR for more than one release cycle. Conflicts compound across versions.
- Don't hand-edit files under automation configuration, `.nf-core.yml`, or `ro-crate-metadata.json` unless you're explicitly resolving a template-sync conflict — they'll be overwritten by the next sync.
- Don't pin tool versions outside containers/conda specs. The container directive owns the version.
- Don't use non-MIT licenses. nf-core requires MIT.
- Don't add a CI requirement that depends on manual intervention — pipelines must run end-to-end via a single `nextflow run` command.

## Tooling notes

- nf-core ships with `prek` (Rust-based pre-commit replacement) as of tools 4.0. Run `prek install --overwrite` to migrate from `pre-commit`.
- Strict syntax is required for pipelines targeting tools ≥4.0. Use `migrate-nextflow-code` to migrate.
- `nf-core.yml` carries the `nf_core_version` field — that's the source of truth for which tools version the pipeline was last synced against.
