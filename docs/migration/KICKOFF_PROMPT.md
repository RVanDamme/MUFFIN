# Prompts for Claude Code

The standing rules are in `CLAUDE.md` and load automatically, so these prompts stay short. Open Claude Code at the repository root on branch `Upgrade` (pull first). Real pipeline runs need a machine with Docker or Apptainer; planning, scaffolding and linting work anywhere.

## 1. First session (Stage 1)

```
You are rebuilding MUFFIN as an nf-core-style pipeline on branch Upgrade.
CLAUDE.md has the rules; you do every stage and its review yourself.

Before anything else:
1. git pull, then read docs/migration/PROGRESS.md, STAGES.md,
   REPORT_ADDENDUM.md and MUFFIN_nfcore_migration_report.md completely.
2. Run .claude/scripts/setup-toolchain.sh if needed, then print the versions of
   nextflow, nf-core, nf-test and which container engine works. If none works,
   list the Stage 1 acceptance checks you cannot run here.
3. List the nextflow plugin skills you can see, and confirm you can push to
   origin Upgrade (git push --dry-run).

Then:
4. Ask me the open decisions for Stage 1 (N1, N2, N10, the template features)
   and N3, which blocks Stage 2, in one message, each with your recommendation.
5. Print the Stage 1 plan: the task list (1.1, 1.2, ...), every file you will
   create, move or modify, how you resolve collisions between the template and
   existing files, the commands, and the check for each acceptance line.

Do not modify any file until I reply "go". After "go": record my answers and
the task checklist in PROGRESS.md and push that first. Then do Stage 1 only,
checkpointing after every task exactly as CLAUDE.md describes (lint, PROGRESS,
log, commit, push). Finish with the stage-reviewer pass, the saved stage
report and the stage-1-done tag.
```

## 2. Resume after an interruption (usage limit, timeout, new machine)

```
Resume the MUFFIN rebuild. Follow the "Resuming" steps in CLAUDE.md: fetch,
check git status and the last 15 commits, read the Now block in
docs/migration/PROGRESS.md and the tail of the current stage log. Tell me in
three lines where things stand and what you will do next, then continue from
the first unticked task. Keep checkpointing after every task.
```

Slash command: `/resume`. If the previous session ran in the cloud and its sandbox is gone, nothing is lost beyond the task in progress, because every checkpoint was pushed.

## 3. Start the next stage

```
Continue the MUFFIN rebuild with Stage <N>. Read docs/migration/PROGRESS.md
first, then the Stage <N> section of STAGES.md and what it refers to. Print the
toolchain versions and the working container engine, ask any open decision
this stage needs, print the plan with its task list, and write the checklist
to PROGRESS.md. Proceed unless the stage has a ⚠️ item. Checkpoint after every
task; end with the review, the saved report and the tag.
```

Slash command: `/stage <N>`.

## If you see the usage limit coming

Send: `Checkpoint now and stop.` Claude finishes or WIP-commits the current task, updates the Now block, pushes, and stops. Next time, use prompt 2.

## Where to look

| Question | File |
|---|---|
| Where are we, what is next? | `docs/migration/PROGRESS.md` (Now block, checklist) |
| What was run and what did it show? | `docs/migration/log/stage-<N>.md` |
| What did a stage deliver, what did the review find? | `docs/migration/reports/stage-<N>.md` |
| Why was something decided? | decision tables in `PROGRESS.md`, options in `REPORT_ADDENDUM.md` §6 |
| What changed for users? | `CHANGELOG.md` |
| Restore points | git tags `stage-<N>-done` |
