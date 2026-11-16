# Prompts for Claude Code

The standing rules are in `CLAUDE.md` and load automatically, so these prompts stay short. Open Claude Code at the repository root on branch `Upgrade`, on a machine with Docker or Apptainer if you want real runs.

## First session (Stage 1)

```
You are rebuilding MUFFIN as an nf-core-style pipeline. CLAUDE.md has the rules.

Before anything else:
1. Read docs/migration/PROGRESS.md, STAGES.md, REPORT_ADDENDUM.md and
   MUFFIN_nfcore_migration_report.md completely.
2. Run .claude/scripts/setup-toolchain.sh if needed, then print the versions of
   nextflow, nf-core, nf-test and which container engine works. If none works,
   say which Stage 1 acceptance checks you will not be able to run.
3. Confirm that the nextflow plugin skills are available (list them).

Then:
4. Ask me the open decisions needed for Stage 1 (N1, N2, N5, N10 and the
   template features) in one message, each with your recommendation.
5. Print the Stage 1 plan: goal, every file you will create, move or modify,
   how you will resolve collisions between the template and existing files,
   the commands you will run, and how you will verify each acceptance line.

Do not modify any file until I reply "go". After "go", do Stage 1 only, record
my answers in PROGRESS.md, and finish with the stage report.
```

## Every later session

```
Continue the MUFFIN rebuild with Stage <N>.
Read docs/migration/PROGRESS.md first, then the Stage <N> section of STAGES.md
and the parts of REPORT_ADDENDUM.md and the report it refers to.
Print the toolchain versions and the working container engine, then the stage
plan. Ask me any open decision this stage needs. Proceed after the plan unless
the stage has a ⚠️ item. Do this stage only; end with the stage report and an
updated PROGRESS.md.
```

The same thing is available as the slash command `/stage <N>`.

## Review pass (for Seqera Co-Scientist or a second model)

```
Review the last stage of the MUFFIN nf-core rebuild on branch Upgrade. Do not
change code. Check the diff of the stage's commits against CLAUDE.md,
docs/migration/STAGES.md (acceptance criteria) and REPORT_ADDENDUM.md.
Report, with file and line: module calls whose arguments do not match the
installed module's input block; channel joins whose keys or shapes are wrong;
ext.args that do not exist in the pinned tool version; behaviour changes
missing from CHANGELOG.md; anything claimed as verified without evidence.
Rank findings by severity and say which ones block the next stage.
```
