# Stage 1 log — Scaffold

Environment: Claude Code cloud session; nextflow 26.10.6, nf-core/tools 4.2.2, Java 21; nf-test NOT installed (download 403); Docker client present, no daemon → no container runs here.

## 1.0 — decisions, checklist, skills (2026-11-16)
- Recorded N1, N2, N10, N11, N13–N18 in `PROGRESS.md`; version changed 3.0.0 → 2.0.0 in `STAGES.md`, `PROGRESS.md`, `CLAUDE.md`.
- `npx --yes skills add seqeralabs/nextflow-plugin --yes` → "Done!", 28 skills under `.agents/skills/`, symlinks in `.claude/skills/`, `skills-lock.json`.
- Not verified: skills are only loaded by a new session; this session worked from `CLAUDE.md`.
- Previous entry commit: 154ede9

## 1.1 (2026-11-16)
- `git mv main.nf nextflow.config modules readflows configs obsolete_process MUFFIN_params.yml bin legacy/` (directory names preserved).
- Kept in place: `.figure/`, `LICENSE`, `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, `README.md`, `CLAUDE.md`, `.claude/`, `.agents/`, `docs/migration/`.
- `.agents/` (installed skills, contain example .nf files with 1 lint error) and `legacy/` are excluded from lint; not our code.
- `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → Linting Nextflow code.. Nextflow linting complete!  No files found to process 
- Previous entry commit: be238b2

## 1.2 (2026-11-16)
- `nf-core pipelines create --template-yaml template.yml` in /tmp/scratch (yaml saved as `docs/migration/template.yml`); files copied with tar, without `.git` and without the template `README.md`.
- Collisions: `README.md` kept (v1 docs), new `README_v2.md` holds v2 status; `LICENSE`, `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md` kept (license feature skipped; template has only `docs/CONTRIBUTING.md`); `.gitignore` = template + `results*/`; `docs/` no collision.
- `.prettierignore` extended with legacy/, .agents/, .claude/, docs/migration/ and the kept root files.
- `git push origin TEMPLATE` from the scratch repo → `* [new branch] TEMPLATE -> TEMPLATE` (98f9ae6).
- Note: nf-core create offered an ORCID found online for the author; not used (template leaves `orcid: ''`).
- `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → Nextflow linting complete!  ??  2 files had 2 warnings  ? 27 files had no errors 
- Previous entry commit: 291707c

## 1.3 (2026-11-16)
- `manifest.nextflowVersion = '!>=26.04.0'`; CI minimum version in `.github/workflows/{nf-test,fix-linting}.yml` set to 26.04.0 to match.
- Contributor entry filled (github, author/maintainer); removed `validation.defaultIgnoreParams = ["genomes"]` (igenomes skipped). No `params.fasta`/`getGenomeAttribute` was generated.
- `workflows/muffin.nf` emits the samplesheet channel (removes the unused-parameter warning; placeholder until Stage 4).
- Remaining warning: `when:` in `modules/nf-core/multiqc` (expected, not ours).
- `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → Nextflow linting complete!  ??  1 file had 1 warning  ? 27 files had no errors 
- Previous entry commit: c7cef5c

## 1.4 (2026-11-16)
- `assets/schema_input.json`: `sample, short_reads_1, short_reads_2, long_reads, rna_reads_1, rna_reads_2`; at least one of short pair / long reads (`anyOf`); pairs complete (`dependentRequired`); unique sample names.
- `assets/schema_coverage.json` (`sample, short_reads_1, short_reads_2, long_reads`), `assets/schema_bins.json` (`sample, bins_dir`); `assets/samplesheet.csv` example updated.
- `PIPELINE_INITIALISATION` takes `coverage_input` and `input_bins`, derives `meta.mode` (short|long|hybrid), `meta.has_rna`, `meta.single_end=false`, emits `samplesheet`, `coverage`, `bins`; fails early if neither `--input` nor `--input_bins` is given. Template single-end/multi-run grouping removed (one row per sample).
- Not verified yet: schema behaviour on real samplesheets (checked in 1.8 with a run without containers).
- `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → Nextflow linting complete!  ??  1 file had 1 warning  ? 27 files had no errors 
- Previous entry commit: 539b978

## 1.5 (2026-11-16)
- 31 pipeline parameters added to `nextflow.config` (report §5.3 + addendum §4a: `bin_taxonomy`, `sourmash_lca_db`, `gtdbtk_db`). Defaults follow legacy (`hybrid_assembler=metaspades`, `binner=metabat2`, `semibin_mode=pretrained`, `semibin_environment=global`, `pilon_rounds=2`, `bin_min_score=50`). `medaka_model` defaults to unset = medaka auto-selection (still *verify*, Stage 2).
- `nf-core pipelines schema build --no-prompts` → "Pipeline schema looks valid (found 53 params)"; groups, descriptions, enums and the GTDB-Tk cost text filled by editing the JSON. `input` no longer required (alternative `--input_bins`; checked in `validateInputParameters`).
- `nf-core pipelines schema lint` → "[✓] Default parameters match schema validation", "[✓] Pipeline schema looks valid (found 53 params)".
- GTDB-Tk figures in the help text are orders of magnitude, not yet checked against the pinned release (Stage 2 item).
- `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → Nextflow linting complete!  ??  1 file had 1 warning  ? 27 files had no errors 
- Previous entry commit: 9a9adce

## 1.6 (2026-11-16)
- `assets/samplesheet_test.csv`: minigut_hybrid (short pair 3.3+3.4 MB + ONT 10.3 MB), minigut_short (1.6+1.8 MB), minigut_long (same ONT file); `conf/test.config` and `conf/test_full.config` (placeholder, TODO Stage 6) point at it.
- Executed WITHOUT a container engine (stub only, not the acceptance command): `nextflow run . -profile test -stub-run --outdir results_stub` → `[PROCESS] RVANDAMME_MUFFIN:MUFFIN:MULTIQC (muffin)`, `[SUCCESS] completed=1 failed=0 cached=0`.
- Channel check (temporary `.view`, reverted): rows came out as `[[id:minigut_hybrid, mode:hybrid, has_rna:false, single_end:false], [R1, R2], ont, []]`, `[[id:minigut_short, mode:short, …], [R1, R2], [], []]`, `[[id:minigut_long, mode:long, …], [], ont, []]`.
- Negative checks: samplesheet with R1 only and with no reads → "Validation of file failed … Object does not have some of the required properties [short_reads_2]" / "Value does not match against any of the schemas"; no `--input`/`--input_bins` → "Please provide a samplesheet with --input (reads) or --input_bins (existing bins)."
- Fixes found by these runs: `bin_min_score` schema type integer/default 50 (was shown as non-default); samplesheet `errorMessage` as plain strings.
- Observed, not fixed: with `-stub-run` on Nextflow 26.10.6 the template check warns "positional argument `true` has been detected" (no warning with `-preview`); comes from `subworkflows/nf-core/utils_nfcore_pipeline`, harmless.
- Not verified: RNA columns and `--coverage_input`/`--input_bins` parsing with real files (no test data yet).
- `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → Nextflow linting complete!  ??  1 file had 1 warning  ? 27 files had no errors 
- Previous entry commit: 11ea817

## 1.7 (2026-11-16)
- `CHANGELOG.md` rewritten for v2.0.0dev: breaking samplesheet interface (D2), renamed/removed parameters (report §5.3), taxonomy methods, minimum Nextflow. Behaviour fixes of later stages get their lines when implemented.
- `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → Nextflow linting complete!  ??  1 file had 1 warning  ? 27 files had no errors 
- Previous entry commit: 276e96c
