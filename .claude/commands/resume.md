---
description: Resume the MUFFIN rebuild after an interrupted session
---

Resume the MUFFIN rebuild following the "Resuming" steps in `CLAUDE.md`:

1. `git fetch && git status && git log --oneline -15`.
2. Read the **Now** block and task checklist in `docs/migration/PROGRESS.md` and the tail of the current `docs/migration/log/stage-<N>.md`.
3. Handle uncommitted changes as `CLAUDE.md` says (keep if they belong to the current task and lint passes, otherwise stash and report). Finish any `wip(...)` commit first.
4. Print the toolchain versions and the working container engine.
5. Say in three lines where things stand and what you do next, then continue from the first unticked task, checkpointing after every task.
