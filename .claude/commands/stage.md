---
description: Start or continue a stage of the MUFFIN nf-core rebuild
argument-hint: <stage number>
---

Continue the MUFFIN rebuild with Stage $ARGUMENTS.

1. Read `docs/migration/PROGRESS.md`, then the Stage $ARGUMENTS section of `docs/migration/STAGES.md` and the parts of `docs/migration/REPORT_ADDENDUM.md` and `docs/migration/MUFFIN_nfcore_migration_report.md` it refers to.
2. Print the versions of `nextflow`, `nf-core`, `nf-test` and which container engine works (run `.claude/scripts/setup-toolchain.sh` if something is missing). If no engine works, list the acceptance checks you cannot run.
3. If earlier stages are not marked DONE in `PROGRESS.md`, say so and stop.
4. Ask any open decision this stage needs, each with a recommendation.
5. Print the stage plan: goal, files to create or modify, commands, verification per acceptance line. Wait for "go" if this is Stage 1 or the stage has a ⚠️ item; otherwise proceed.
6. Do this stage only. Finish with the stage report from `STAGES.md` and update `PROGRESS.md`.
