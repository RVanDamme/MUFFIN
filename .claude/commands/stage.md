---
description: Start or continue a stage of the MUFFIN nf-core rebuild
argument-hint: <stage number>
---

Continue the MUFFIN rebuild with Stage $ARGUMENTS.

1. `git pull`. Read `docs/migration/PROGRESS.md`, then the Stage $ARGUMENTS section of `docs/migration/STAGES.md` and the parts of `docs/migration/REPORT_ADDENDUM.md` and `docs/migration/MUFFIN_nfcore_migration_report.md` it refers to.
2. If the **Now** block shows this stage already in progress, switch to the "Resuming" steps in `CLAUDE.md` instead of re-planning.
3. Print the versions of `nextflow`, `nf-core`, `nf-test` and which container engine works (run `.claude/scripts/setup-toolchain.sh` if something is missing). If no engine works, list the acceptance checks you cannot run.
4. If earlier stages are not marked DONE in `PROGRESS.md`, say so and stop.
5. Ask any open decision this stage needs, each with a recommendation.
6. Print the stage plan with its task list. Wait for "go" if this is Stage 1 or the stage has a ⚠️ item; otherwise proceed.
7. Write the task checklist to `PROGRESS.md`, commit and push it, then work task by task with a checkpoint (lint, PROGRESS, log, commit, push) after each.
8. Finish with the `stage-reviewer` pass, the stage report saved to `docs/migration/reports/stage-$ARGUMENTS.md`, and the `stage-$ARGUMENTS-done` tag.
