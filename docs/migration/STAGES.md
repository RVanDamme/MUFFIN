# Stages and acceptance criteria

Standing rules live in `/CLAUDE.md`. Facts and corrections live in `REPORT_ADDENDUM.md`. This file says what each stage delivers. `<engine>` is the container engine available on the machine (`docker`, `apptainer` or `singularity`).

A stage is **DONE** only when every acceptance line was run and its output shown. If no container engine is available, finish the static part, mark the stage **DONE (static) — runs pending**, and list the exact commands the user must run.

⚠️ = stop and wait for the user's answer.

---

## Stage 0 — Baseline (optional, decision N9)

Run legacy MUFFIN on the OSF subset (`https://osf.io/9xmh4/download`) with Nextflow 25.04 or `NXF_SYNTAX_PARSER=v1`, modes short / long / hybrid-metaspades / hybrid-metaflye, with the workarounds in report §8 and §11.

- Deliver `baseline/metrics.tsv` (metrics listed in report §8, Stage 0) and notes on which metrics are affected by the bugs in report §3.
- If the legacy pipeline cannot be made to run in about half a day, stop, write down what failed in `baseline/NOTES.md`, and record "no baseline" in `PROGRESS.md`.

## Stage 1 — Scaffold

⚠️ Confirm before generating: decisions N1, N2, N5, N10 and the template features. Proposed `skip_features`: `igenomes`, `fastqc`, `is_nfcore`, `email`, `seqera_platform`, `gpu`, `codespaces`, plus `license` if GPL-3.0 is kept.

1. `git mv` all legacy code (`main.nf`, `nextflow.config`, `modules/`, `readflows/`, `configs/`, `obsolete_process/`, `MUFFIN_params.yml`) into `legacy/`. Keep `bin/`, `.figure/`, `LICENSE`, `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, `README.md`, `CLAUDE.md`, `.claude/`, `docs/migration/` in place. Commit.
2. Generate the template into a **scratch directory outside the repo** with `nf-core pipelines create --template-yaml <file>` (name `muffin`, organisation `RVanDamme`, version `3.0.0dev`). Copy its files, not its `.git`, into the repo. Resolve collisions explicitly (`README.md`, `.gitignore`, `docs/`, `LICENSE`) and list each resolution in the plan. If N10 is yes, tell the user how to push the scratch repo's `TEMPLATE` branch.
3. Remove template leftovers that do not apply (`params.fasta` / `getGenomeAttribute` if igenomes is skipped). Set `manifest.nextflowVersion` per N1.
4. Samplesheet schemas per report §5.2: `assets/schema_input.json` (`sample, short_reads_1, short_reads_2, long_reads, rna_reads_1, rna_reads_2`; at least one of short pair / long; pairs complete), `assets/schema_coverage.json`, `assets/schema_bins.json`. In `PIPELINE_INITIALISATION` derive `meta.mode`, `meta.has_rna`, `meta.single_end`.
5. All parameters of report §5.3 in `nextflow.config` with defaults; rebuild `nextflow_schema.json` (`nf-core pipelines schema build --no-prompts`, then fill descriptions/enums by editing the JSON).
6. `conf/test.config` pointing at small hosted data. Propose the dataset (nf-core/test-datasets is reachable; the OSF subset is larger) and its size before wiring it.
7. Start `CHANGELOG.md` at `v3.0.0dev` with the breaking CLI change (D2).

**Accept**
- `nextflow lint .` → 0 errors; 0 warnings outside `modules/nf-core/**` and `legacy/` (exclude `legacy/` from lint or state its count separately).
- `nf-core pipelines lint` → no failures other than listed, expected template TODOs. Any lint test disabled in `.nf-core.yml` is listed with its reason.
- `nextflow run . -profile test,<engine> -stub --outdir results_stub` succeeds with an empty `MUFFIN` workflow.

## Stage 2 — nf-core modules and `ext.args`

1. `nf-core modules install` for: chopper, fastp, spades, flye, minimap2/align, racon, medaka, pilon, gunzip, bwa/index, bwa/mem, samtools/index, samtools/fastq, metabat2/jgisummarizebamcontigdepths, metabat2/metabat2, semibin/singleeasybin, comebin/runcomebin, seqkit/stats, cat/cat, checkm2/predict, checkm2/databasedownload, sourmash/sketch, sourmash/gather, eggnogmapper, trinity, salmon/index, salmon/quant, untar, multiqc. Add samtools/merge only if N4 is "no"; add sourmash/taxannotate or gtdbtk/classifywf only if N3 says so.
2. For each module, put in the stage report: real `input:` / `output:` copied from the installed `main.nf` → report §6 row → `ext.args` set. Flag every difference from `REPORT_ADDENDUM.md` §3.
3. Write `conf/modules.config`: `withName` selectors scoped by subworkflow where a module is reused; flags from report §6 as corrected by the addendum (SemiBin options in `ext.args2`; no `--compression none`; eggNOG search mode is an input; salmon `-l A`; sourmash `-p k=31,scaled=10000`, ONT 5000; `gather -k 31`; chopper `-q/-l/--headcrop 100`; flye `--meta`; spades `--meta`).
4. ⚠️ Verify in the pinned containers, report, and wait if something is missing: medaka model list and auto-selection; flye `--plasmids` (expected: drop); SemiBin2 `--environment` / `--training-mode`; COMEBin `-l`; CheckM2 extension flag with gzipped bins; `SALMON_QUANT` with `gtf = []`. Without a container engine, check upstream docs/source at the pinned version and mark each item "source-checked, not executed".
5. Regenerate container configs if the template requires it.

**Accept**: `modules.json` complete; `nf-core modules lint --all` clean for `modules/nf-core`; `nextflow lint .` → 0 errors.

## Stage 3 — Local modules

`METAQUAST`, `SOURMASH_TAXMETAGENOME`, `SOURMASH_TAXGENOME` (per N3), `PANKEGG`, and only if downloads are wanted `SOURMASH_DOWNLOADDB`, `EGGNOGMAPPER_DOWNLOADDB`. Each via `nf-core modules create`, with `meta.yml`, `environment.yml` (exact versions), pinned container, `stub:`, versions topic, and an nf-test with a real and a stub case.

- Before writing each one, run the `search-existing-modules` check once more; nf-core/modules moves fast.
- `PANKEGG`: standard library only, so use a pinned plain Python container. ⚠️ N6 before deciding how the KEGG REST call is handled. Header/footer stripping of eggNOG output happens inside the parser.

**Accept**: `nf-test test modules/local/` passes (real + stub); `nf-core modules lint --local` passes.

## Stage 4 — Subworkflows and main workflow

One subworkflow per directory under `subworkflows/local/<name>/` with `main.nf` and `meta.yml`, shapes as in report §7 as corrected by the addendum.

- Routing by `meta.mode` with `branch`: short → FASTP → SPADES `--meta`; long → CHOPPER → FLYE → LONGREAD_POLISH; hybrid → FASTP + CHOPPER → SPADES (hybrid) or FLYE → LONGREAD_POLISH according to `params.hybrid_assembler`.
- `READ_MAPPING`: own reads always; every `--coverage_input` row mapped to **every** assembly with `combine`; BAMs grouped per assembly with `groupTuple` on `meta.id`.
- `BINNING`: one binner per run (D7). MetaBAT2 gets all BAMs through the depth step. SemiBin2/COMEBin input per N4. COMEBin `-l` from N50 parsed out of `SEQKIT_STATS` into `meta.n50` with `map`.
- `BIN_QC`: CheckM2 DB channel exists in both the `--checkm2_db` and the `--download_checkm2_db` case, in the `[dbmeta, db]` shape. Filtering is channel operators on the CheckM2 TSV (`Completeness - 5*Contamination > params.bin_min_score`), matching bins by name with `.fa`, `.fa.gz`, `.fna`, `.fna.gz` stripped.
- `BIN_POLISH_PILON`: only `meta.mode == 'hybrid' && params.hybrid_assembler == 'metaflye'`; `GUNZIP` first; up to two aliased rounds gated by `params.pilon_rounds`.
- `UNBINNED_READ_CLASSIFICATION`, `BIN_TAXONOMY`, `BIN_ANNOTATION`, `TRANSCRIPTOMICS`, `REPORTING` per report §6/§7. `--input_bins` starts at `BIN_QC`; `--skip_*` flags per §5.3.
- Versions through the `versions` topic into MultiQC using the template utilities.

**Accept**
- `nextflow lint .` → 0 errors, 0 warnings in owned code.
- `-profile test,<engine> -stub` passes for short, long, hybrid-metaspades, hybrid-metaflye and `--input_bins`.
- Real `-profile test,<engine>` passes for at least short and long.
- The DAG (`-with-dag`) or `nextflow inspect` output is attached to the stage report.

## Stage 5 — Publishing

`publishDir` rules in `conf/modules.config`, reproducing the legacy layout in spirit (`<sample>/assemble|classify|annotate/…`). Write `docs/migration/OUTPUT_MAP.md` (old path → new path) and `docs/output.md`.

**Accept**: every path in `OUTPUT_MAP.md` exists after the real test run; nothing is written outside `--outdir` and the work directory.

## Stage 6 — Tests and CI

Pipeline-level nf-test per mode (`tests/*.nf.test`) with snapshots of stable outputs only; GitHub Actions running `nextflow lint`, `nf-core pipelines lint`, and the nf-test matrix. DB-heavy steps use stub or tiny DB subsets in `test`; full DBs only in `test_full`.

**Accept**: `nf-test test` passes locally; CI green on the branch.

## Stage 7 — Docs and release preparation

`README.md`, `docs/usage.md` (with the old → new parameter table of report §5.3), `docs/output.md`, `CITATIONS.md`, `CHANGELOG.md` for v3.0.0. ⚠️ Ask before deleting `legacy/`.

**Accept**: `nf-core pipelines lint` → 0 failures.

## Stage 8 — Comparison and (optional) Platform

Compare against the baseline within the tolerances of report §9, or against documented expectations if there is no baseline. Seqera Platform / Launchpad steps only if N8 says so.

---

## Stage report format

```
## Stage N report
Status: DONE | DONE (static) — runs pending | BLOCKED (reason)
Environment: nextflow <ver>, nf-core <ver>, nf-test <ver>, engine <name|none>
Commits: <hash> <message> (one per line)
Files changed: <list>
Verification (command → result, with the relevant output lines):
  - nextflow lint .            → <errors>/<warnings> (owned code: <n>)
  - nf-core pipelines lint     → <passed>/<warned>/<failed>
  - nf-core modules lint --all → …
  - runs / nf-test             → <command> → <result>
Not verified (and why): <list or "none">
Behaviour changes (→ CHANGELOG): <list or "none">
Deviations from report/addendum: <list with file:line evidence, or "none">
Open questions: <list or "none">
Next stage plan: <3–6 bullets>
```
