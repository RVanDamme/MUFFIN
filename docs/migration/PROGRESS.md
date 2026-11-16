# Migration progress

Single source of truth for where the rebuild stands. Updated at every checkpoint (see `CLAUDE.md`, "Traceability and resuming"). A new session reads the **Now** block first.

## Now

- **Stage:** 2 — nf-core modules + ext.args (in progress, started 2026-11-16)
- **Current task:** next after 2.1
- **Next action:** install modules of 2.2
- **In flight / half-decided:** nothing. HPC2N memory request behaviour (`--mem` vs per-core memory) to be checked on the first real run (`docs/hpc2n.md`).
- **Last checkpoint:** task 2.1 (2026-11-16)

## Stage status

| Stage | Status | Tag | Log | Notes |
|---|---|---|---|---|
| Project setup | DONE 2026-11-16 | — | — | No pipeline code changed |
| 0 Baseline | not started | | | optional, see N9 |
| 1 Scaffold | DONE 2026-11-16 | `stage-1-done` (a19dd17; later Stage 1 commits are docs + `hpc2n` profile) | `log/stage-1.md` | cloud session: no container engine, no nf-test |
| 2 nf-core modules + ext.args | in progress | | `log/stage-2.md` | cloud session: no container engine |
| 3 Local modules | not started | | | PANKEGG tax-genome test waits on N12 |
| 4 Subworkflows + workflow | not started | | | needs N4, N7, N11 |
| 5 Publishing | not started | | | |
| 6 Tests + CI | not started | | | |
| 7 Docs | not started | | | |
| 8 Comparison / Platform | not started | | | needs N8 |

## Task checklist — current stage (Stage 2)

_Format: `- [ ] 2.1 <task> — <commit hash when done>`._

- [x] 2.1 `nf-core modules update multiqc` (fixes the versions-topic lint failure) and adapt the call if its inputs changed
- [ ] 2.2 Install read QC / assembly / polishing modules: chopper, fastp, spades, flye, minimap2/align, racon, medaka, pilon, gunzip
- [ ] 2.3 Install mapping / binning modules: bwa/index, bwa/mem, samtools/index, samtools/fastq, metabat2/jgisummarizebamcontigdepths, metabat2/metabat2, semibin/singleeasybin, comebin/runcomebin, seqkit/stats, cat/cat
- [ ] 2.4 Install QC / taxonomy / annotation / RNA modules: checkm2/predict, checkm2/databasedownload, sourmash/sketch, sourmash/gather, gtdbtk/classifywf, eggnogmapper, trinity, salmon/index, salmon/quant, untar
- [ ] 2.5 Real signatures of every installed module vs addendum §3 (table in the stage report)
- [ ] 2.6 `conf/modules.config`: scoped `withName` selectors and `ext.args`
- [ ] 2.7 ⚠️ Tool-option checks (source-checked, no containers here) — report and wait if something is missing
- [ ] 2.8 Container configs, lints, stage-reviewer, report, README_v2, tag

`samtools/merge` is not installed: N4 is taken as "yes" (recommended) for the module list; it is formally confirmed in Stage 4.

### Stage 1 checklist (done)

- [x] 1.0 Record decisions and checklist; install `seqeralabs/nextflow-plugin` skills (`npx skills add`) — be238b2
- [x] 1.1 `git mv` legacy code into `legacy/`, keeping the same directory names; files still needed stay in place — 291707c
- [x] 1.2 Generate template (scratch dir, `2.0.0dev`), copy into repo, resolve collisions, push `TEMPLATE` branch, create `README_v2.md` — c7cef5c
- [x] 1.3 Remove template leftovers; `manifest.nextflowVersion`; exclude `legacy/` from lint/prettier — 539b978
- [x] 1.4 Samplesheet schemas (`schema_input`, `schema_coverage`, `schema_bins`) and `meta` derivation in `PIPELINE_INITIALISATION` — 9a9adce
- [x] 1.5 Parameters (report §5.3 + addendum §4a) and `nextflow_schema.json` — 11ea817
- [x] 1.6 `conf/test.config` + `assets/samplesheet_test.csv` (nf-core/test-datasets `mag` minigut, ~20 MB) — 276e96c
- [x] 1.7 `CHANGELOG.md` at `v2.0.0dev` — 5c81f12
- [x] 1.9 `hpc2n` profile (`conf/hpc2n.config`, `docs/hpc2n.md`)
- [x] 1.8 Acceptance checks, stage-reviewer pass, stage report, tag `stage-1-done` (tag withheld until the container stub run passes)

## Decisions

Report decisions (defaults stand until changed here):

| # | Decision | Status |
|---|---|---|
| D1 | Standalone nf-core-style pipeline `RVanDamme/muffin` | default |
| D2 | Samplesheet `--input` replaces `--illumina/--ont/--rna/--mode` (v2.0.0) | default |
| D3 | Metric-level equivalence | default |
| D4 | User-supplied DB paths, opt-in downloads | default |
| D5 | Local modules first, upstream later | default |
| D6 | sourmash `lca classify` → `gather` + `tax` | **replaced by N3** (three methods) |
| D7 | One binner per run | default |
| D8 | Publishing as the template generates = `publishDir` | default |
| — | Short mode uses `spades --meta` (legacy ran isolate SPAdes) | default — needs explicit sign-off (behaviour change) |

Addendum decisions (options in `REPORT_ADDENDUM.md` §6):

| # | Question | Answer | Date |
|---|---|---|---|
| N1 | Minimum Nextflow version | **`!>=26.04.0`** (confirmed by the user) | 2026-11-16 |
| N2 | Licence GPL-3.0 vs MIT | **Keep GPL-3.0** (`license` template feature skipped) | 2026-11-16 |
| N3 | Bin taxonomy method | **`--bin_taxonomy`: `sourmash` (gather + tax genome, default), `sourmash_lca` (legacy option), `gtdbtk` (higher accuracy, heavy resources, documented). `PANKEGG_MAKEDB` supports all three** (addendum §4a) | 2026-11-16 |
| N4 | All BAMs to COMEBin/SemiBin, drop merge | OPEN (recommended yes) | |
| N5 | Agents | **Claude does every stage and the review** | 2026-11-16 |
| N6 | PANKEGG | **Standalone PANKEGG, `pankegg_make_db` only; legacy `bin/` scripts not ported** | 2026-11-16 |
| N7 | GPU for SemiBin2/COMEBin | OPEN (recommended CPU only) | |
| N8 | Seqera Platform stage | OPEN (recommended skip) | |
| N9 | Legacy baseline run | OPEN | |
| N10 | Keep a `TEMPLATE` branch | **Yes** | 2026-11-16 |
| N11 | RNA results reporting without a PANKEGG RNA mode | **Publish Trinity / salmon / eggNOG outputs + MultiQC; no RNA in the PANKEGG DB** | 2026-11-16 |
| N12 | PANKEGG release reading `sourmash tax genome` output (tag + flag) | OPEN — user is updating PANKEGG; default method skips the PANKEGG step with a warning until then | |
| N13 | Version of the rebuild | **`2.0.0dev` → v2.0.0** (not v3; overrides the report and earlier docs) | 2026-11-16 |
| N14 | Template features | **skip `igenomes fastqc is_nfcore email seqera_platform gpu codespaces license`; keep the rest** | 2026-11-16 |
| N15 | Test data | **nf-core/test-datasets `mag` minigut (2 short pairs + 1 ONT file, ~20 MB); no RNA test data yet** | 2026-11-16 |
| N16 | Human-readable status | **`README_v2.md` at the repo root, updated at every stage; replaces `README.md` in Stage 7** | 2026-11-16 |
| N17 | Legacy layout | **Files not needed any more move to `legacy/<same dirname>/`; files still needed stay where they are** | 2026-11-16 |
| N18 | Skills | **`npx skills add seqeralabs/nextflow-plugin` run; skills committed under `.agents/skills/` (symlinked from `.claude/skills/`); keep them for the duration of the upgrade and remove `.agents/`, `.claude/skills/`, `skills-lock.json` before merging to the main branch** | 2026-11-16 |
| N19 | Run environment | **HPC2N Kebnekaise with Apptainer (Slurm); `-profile hpc2n`, account from `SBATCH_ACCOUNT`** | 2026-11-16 |
| — | Push policy | **Push to `origin Upgrade` after every checkpoint commit; never force-push** | 2026-11-16 |

## Verification still owed

**Carried to Stage 2:** `nf-core modules lint --all` fails for `multiqc` ("does not emit software versions as topic") → `nf-core modules update multiqc`; `tests/nextflow.config` still points at the template's `test-datasets/muffin/` path; decide `bin_min_score` integer vs number; handle `--input` together with `--input_bins` (Stage 4).

Items marked *(verify)* in `REPORT_ADDENDUM.md`: medaka model list / auto-selection; CheckM2 extension flag and `Name` column with gzipped bins; `SALMON_QUANT` with empty `gtf`; whether container configs must be regenerated after module installs; PANKEGG availability on Bioconda/PyPI; `sourmash lca classify` in sourmash 4.9.4; `tax genome` csv_summary columns; GTDB-Tk DB size / memory for the pinned release and its handling of gzipped bins; bin-name agreement between eggNOG prefix, CheckM2 `Name` and classification files.

## Session log

One line per session: date — what was done — last commit.

- 2026-11-16 — Project configured; report re-checked on Nextflow 26.10.6 / nf-core tools 4.2.2; corrections in `REPORT_ADDENDUM.md`. PANKEGG switched to the standalone tool; traceability protocol and reviewer added. Bin taxonomy set to three methods (N3).
- 2026-11-16 — Stage 1 started (cloud session, no container engine). Decisions N1, N2, N10, N11, N13–N18 recorded.
- 2026-11-16 — Stage 1 static work complete (be238b2..HEAD); reviewer pass done; container stub run and tag pending.
- 2026-11-16 — User confirmed N1 and the skills decision; pushed tag `stage-1-done`; container stub run still owed (command was run with a literal `<engine>`).
- 2026-11-16 — Container stub run passed on Kebnekaise (user); Stage 1 DONE. `hpc2n` profile and `docs/hpc2n.md` added.
