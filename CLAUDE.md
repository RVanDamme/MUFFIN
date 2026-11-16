# MUFFIN — nf-core-style rebuild (branch `Upgrade`)

MUFFIN is a metagenomics pipeline (hybrid assembly, differential-coverage binning, bin QC, taxonomy, annotation, optional metatranscriptomics). The legacy DSL2 code on this branch does not parse on current Nextflow. We are **rebuilding it on the nf-core template**, reusing nf-core modules, on Nextflow ≥ 26.04 with strict syntax. It is an nf-core-*style* pipeline (`RVanDamme/muffin`), not an official nf-core pipeline.

## Read first, every session

1. `docs/migration/PROGRESS.md` — the **Now** block (stage, task, next action), task checklist, decisions. **Start here.**
2. `docs/migration/STAGES.md` — what each stage must deliver and how it is accepted.
3. `docs/migration/REPORT_ADDENDUM.md` — verified tool versions, real module signatures, corrections.
4. `docs/migration/MUFFIN_nfcore_migration_report.md` — the full audit and target design (the spec).

Precedence when sources disagree: the user's latest instruction → `PROGRESS.md` decisions → `REPORT_ADDENDUM.md` → the report → legacy code. Legacy code is evidence of intent, never the spec. If you find a new contradiction, stop, show file and lines, propose a fix, and record the outcome in `PROGRESS.md`.

## Toolchain

| Tool | Version | Notes |
|---|---|---|
| Nextflow | 26.10.6 (`NXF_VER=26.10.6`) | `manifest.nextflowVersion = '!>=26.04.0'` unless `PROGRESS.md` says otherwise |
| nf-core/tools | 4.2.2 | pinned in `.nf-core.yml` |
| nf-test | current | plugins as set by the template |
| Java | 17+ | |
| Containers | Docker or Apptainer/Singularity | use whichever the machine has; never assume Docker |

Run `.claude/scripts/setup-toolchain.sh` if anything is missing, then print the versions. If there is no working container engine (Claude Code cloud sessions have none), say so at the start: you can still edit, lint and install modules, but you must not claim that a run or test passed.

## Commands

```bash
nextflow lint .                                   # after every batch of edits
nf-core pipelines lint                            # end of each stage
nf-core modules lint --all                        # end of each stage
nf-core modules install <tool/sub>                # the only way to add an nf-core module
nf-core modules create <tool/sub>                 # local modules (writes to modules/local/)
nf-core pipelines schema build --no-prompts       # after changing params
nextflow run . -profile test,<engine> -stub --outdir results_stub
nextflow run . -profile test,<engine> --outdir results
nf-test test --profile +<engine>
```

## Rules

**Process**
- Claude does every stage and its review; there is no other agent and no file-ownership split.
- One stage at a time, split into tasks (see "Traceability and resuming"). Print the plan (goal, tasks, files, commands, verification) before editing. Wait for "go" on Stage 1 and on anything `STAGES.md` marks ⚠️.
- Never report success without the command and its output. Stub runs and lint are not a substitute for a real run; say which one you did.
- Stop and report instead of looping when: the same failure survives two targeted fixes; a module signature contradicts the docs; a tool flag does not exist in the pinned version; data or databases are unreachable.
- End every stage with the review pass, then the stage report (format in `STAGES.md`), saved to `docs/migration/reports/stage-<N>.md`.

**nf-core**
- Modules come from `nf-core modules install` and are tracked in `modules.json`. Never hand-edit `modules/nf-core/**` or `subworkflows/nf-core/**`; if a change is unavoidable use `nf-core modules patch` and explain why. No copies of nf-core modules under `modules/local/`.
- Before calling any module, open its installed `main.nf` and `meta.yml` and use the real signature. Never write a call from memory or from the tables in the docs alone.
- Conventions: uppercase process names; `tag "$meta.id"`; `process_*` resource labels; tool options only via `task.ext.args*` / `ext.prefix` in `conf/modules.config`, with `withName` selectors scoped by subworkflow where a module is reused; versions via `topic: versions`; `stub:` in every local module; no `params.*` inside modules; optional inputs are `[]`.
- Containers: whatever each nf-core module pins. Local modules pin an exact version (biocontainers or `community.wave.seqera.io`), never `latest`. Do not reuse `legacy/configs/container.config`.
- Publishing: `publishDir` in `conf/modules.config`, as the template generates. No `output {}` block yet.

**Nextflow language**
- Strict syntax only: no implicit `it`, no `for`/`while`/`switch`, no top-level statements mixed with declarations, no conditional `include`, no `include … params(…)`, no `shell:` blocks, lowercase `channel.` factories, no `nextflow.enable.dsl`.
- Static types, records and typed processes are **deferred**. Do not enable any static-typing feature flag or migrate to workflow outputs.
- `0 errors` from `nextflow lint .` everywhere; `0 warnings` in code we own (`when:` deprecation warnings inside `modules/nf-core/**` are expected and not ours to fix).

**Channels**
- `meta.id` is the join key. `meta` carries `id`, `mode` (`short|long|hybrid`), `has_rna`, `single_end`.
- Before any `join`, `combine`, `groupTuple`, `branch` or `transpose`, write the input shapes, the key and the output shape as a comment.
- Bins travel as `[meta, List<Path>]`; `transpose()` only where a tool is per-bin, adding `meta.bin_id`.
- No `join(remainder: true)` with null padding, no dummy placeholder files, no `collect()` that mixes samples.
- No processes that only glob, rename, `cat` or filter rows; use operators or `CAT_CAT`.
- Scripts in `bin/` are on `PATH`; never `python ${projectDir}/bin/…`. Nothing writes to `params.outdir` from a script.
- The legacy pankegg scripts in `bin/` are not ported. The report DB comes from the standalone PANKEGG tool (`pankegg_make_db`, addendum §4).

**Behaviour**
- Implement the *intended* behaviour (report §3), not the legacy bugs. Every deviation from legacy behaviour gets a `CHANGELOG.md` line under `### Changed` or `### Fixed`, including the ones already decided.
- Legacy code moves to `legacy/` in Stage 1 and stays until the user confirms the comparison is done.

**Git**
- Work on `Upgrade`. Small conventional commits (`feat(binning): …`, `fix(bin_qc): …`), one logical change each, with the task ID in the body (`Task: 4.3`). Never commit `work/`, `.nextflow*`, `results*/`, `.nf-test/`.
- **Push to `origin Upgrade` after every checkpoint commit** (standing permission from the user, so that an interrupted session loses nothing). Never force-push, never rewrite pushed history, never touch `master`, no PRs unless asked.

## Traceability and resuming

Sessions can end without warning (usage limit, timeout, lost sandbox). Anything not committed and pushed is lost, and the next session knows only what is in the repo. Work so that an interruption costs at most one task.

**Tasks.** At plan time, split the stage into tasks of roughly one subworkflow, one module group or one config file each, numbered `<stage>.<n>`. Write them as a checklist in `PROGRESS.md` and commit that before editing anything else.

**Checkpoint = after every task**, in this order:
1. `nextflow lint .` (and the task's own check); record the result.
2. Tick the task in `PROGRESS.md`; update the **Now** block (current task, next action, anything half-decided); add decisions to the decision table with the date.
3. Append to `docs/migration/log/stage-<N>.md`: task ID, what changed, each verification command with its result, anything not verified, and the commit hash of the previous entry.
4. Commit code + `PROGRESS.md` + log together, then push.

**Write before you run.** Before a long command (pipeline run, nf-test, many module installs), update the **Now** block with the command and what you expect, and commit. If the session dies mid-run, the next one knows what was in flight.

**Work in progress.** If a task is too big to finish cleanly, commit what parses with a `wip(<scope>): …` message and say exactly what is missing in the **Now** block. A pushed WIP commit beats a perfect uncommitted one. Never leave more than about 30 minutes of work uncommitted.

**Resuming.** On "resume" or `/resume`:
1. `git fetch && git status && git log --oneline -15`; read the **Now** block and the tail of the current stage log.
2. If there are uncommitted changes, show the diff summary. Keep them only if they belong to the current task and lint passes; otherwise stash them (`git stash push -u -m "orphan <date>"`) and say so. Never discard silently.
3. If the last commit is `wip(...)`, finish that task first.
4. State in three lines where things stand and what you will do next, then continue without waiting unless a ⚠️ item or an open decision blocks.

**Evidence.** A claim of "passes" in a log, report or commit message must be backed by the command and the relevant output lines in the stage log. Things read from source but not executed are labelled "source-checked, not executed".

**Review.** Before a stage is marked DONE, launch the `stage-reviewer` subagent (fresh context, read-only) on the stage's commit range. Put its findings in the stage report, fix everything it marks blocking, and re-run the affected checks. If subagents are unavailable, do the review yourself as a separate pass against the checklist in `.claude/agents/stage-reviewer.md` and say that it was not independent.

**Milestones.** When a stage is DONE and pushed, tag it `stage-<N>-done` and push the tag.

## Skills (plugin `nextflow@nextflow-plugin`, repo `seqeralabs/nextflow-plugin`)

Enabled through `.claude/settings.json`. If they are not listed, install with `npx skills add seqeralabs/nextflow-plugin`. The older `nextflow-io/agent-skills` repo is archived; do not use it.

| Use for | Skill |
|---|---|
| Installing/upgrading Nextflow and Java | `install-nextflow` |
| Strict-syntax and topic-channel questions (stop before its static-typing and workflow-output migrations) | `migrate-nextflow-code` |
| Deciding reuse vs. new module; reading a module's interface | `search-existing-modules`, `run-module` (`nextflow module view nf-core/<name>`) |
| Template, lint and nf-core conventions | `maintain-nf-core-pipeline` |
| `conf/*.config`, profiles, resource labels | `nextflow-config` |
| `nextflow_schema.json` and samplesheet schemas | `nextflow-schema` |
| `join` / `groupTuple` / `collectFile` correctness | `nextflow-output-patterns` |
| Local module containers (Wave) | `create-container` |
| nf-test authoring and failure repair | `nf-test` |
| Failed or hanging runs | `debug-local-run`, `repair-workflow` |

Where a skill and this file disagree, this file wins. Two known conflicts:
- Several skills prefer **registry-managed includes** (`include { X } from 'nf-core/x'`, `nextflow module install`). Do **not** use them here: nf-core lint needs `nf-core modules install`, `modules.json` and relative includes (`'../../modules/nf-core/x/main'`). `nextflow module search/view` is fine for discovery.
- Platform skills (`launch-workflow`, `seqerakit`, `debug-seqera-failed-run`, …) need a Seqera account. Do not use them or the Seqera MCP unless the user asks.
- `migrate-nextflow-code` says to keep behaviour identical. Here behaviour changes are allowed when they are listed in the report, the addendum or `PROGRESS.md`, and recorded in `CHANGELOG.md`.
