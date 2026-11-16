<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Tools 4.0 template merge guide

Source: https://nf-co.re/blog/2026/tools-4_0_0
Video: https://youtu.be/LBB9NueK-vk

## What changed in the template

- **Nextflow minimum bumped to 25.10.4**, strict syntax enforced via pre-commit hook.
- **`prek` replaces `pre-commit`** (Rust rewrite, same `.pre-commit-config.yaml` format). Existing repos need `prek install --overwrite`.
- **Container configs auto-generated**: new container configuration files for Docker, Singularity/Apptainer, and per-arch Conda.
- **Apptainer added** alongside Singularity in module templates.
- **Webhook notifications removed**, replaced by an `nf-slack` plugin invoked from GitHub Actions. Existing `hook_url` / `slackreport` / `adaptivecard` settings must migrate.
- **Contribution guidance updated** with new AI/LLM instructions; compare the user's existing documentation with the new template.
- **RO-Crate metadata** now reads contributors from `manifest.contributors` automatically.
- **Deprecated CLI flags removed**: `--migrate-pytest` gone; prefix-free commands (`nf-core create`) no longer work — must be `nf-core pipelines create`.

## Files that commonly conflict

1. `.nf-core.yml` — `nf_core_version` field
2. `README.md` — Nextflow / template version badges
3. `nextflow.config` — Nextflow version requirement, new parameter blocks, plugin section
4. `.gitignore` — lineage metadata entries added
5. `nf-test.yml` — NFT version bump
6. `awsfulltest.yml` — action ref + nf-slack plugin
7. `multiqc_config.yml` — YAML reformatting (multiline → single-line)
8. `modules.json` — FastQC/MultiQC SHA updates
9. `ro-crate-metadata.json` — auto-regenerated fields
10. Workflow definitions — strict-syntax formatting of process calls
11. Contribution guidance — updated template guidance
12. `.pre-commit-config.yaml` — replaced/updated by prek migration

## Resolution patterns

| File | Strategy |
|---|---|
| `.nf-core.yml` | Take the higher `nf_core_version`. Remove any duplicate keys introduced by the merge. |
| `README.md` | Accept newer version badges unless the pipeline genuinely requires a higher Nextflow than the template offers. |
| `nextflow.config` | Take higher `nextflowVersion`. **Keep** pipeline-specific `params`, `doi`, plugin configs, and profile blocks. Take new template plugin section as-is. |
| `.gitignore` | Accept all template additions; move custom entries to end of file. |
| `nf-test.yml` | Accept template version bumps unless tests require a pinned older NFT. |
| `awsfulltest.yml` | Accept action-ref bump. Migrate any `hook_url` / `slackreport` config to the new `nf-slack` plugin invocation. |
| `multiqc_config.yml` | Accept single-line YAML format — content is identical, only formatting changed. |
| `modules.json` | Accept updated FastQC / MultiQC SHAs from template. Preserve entries for custom or non-template modules. |
| `ro-crate-metadata.json` | Take all template changes — file is auto-generated, no manual content to preserve. |
| Workflow definitions | Accept strict-syntax reformatting. If pipeline has custom logic, ensure it still passes strict syntax — chain `migrate-nextflow-code` if not. |
| Contribution guidance | Merge pipeline-specific guidance into the new template version and remove obsolete duplicates. |
| `.pre-commit-config.yaml` | Accept template's `prek`-compatible config. Run `prek install --overwrite` after merging. |

## Breaking changes / gotchas

- **Run `prek install --overwrite` after merge** to switch the repo's git hooks from `pre-commit` to `prek`.
- **Webhook config migration is mandatory** — existing pipelines with Slack webhooks need to switch to the `nf-slack` plugin. The template removes the webhook section entirely; don't try to keep it.
- **Strict syntax is required**, not optional, for tools 4.0+. After merge, run `nextflow lint` and chain `migrate-nextflow-code` / [schema migration](../../../nextflow-schema/references/nf-schema-migration/README.md) if there are residual failures.
- **No prefix-free `nf-core` commands** — any docs, scripts, or CI invoking `nf-core create`, `nf-core lint`, etc. must be updated to `nf-core pipelines create`, `nf-core pipelines lint`, etc.
- **`--migrate-pytest` is gone**. If you still have pytest-based tests waiting to migrate to `nf-test`, do that conversion outside the sync PR.
