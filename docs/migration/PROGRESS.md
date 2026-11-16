# Migration progress

Single source of truth for where the rebuild stands. Updated at every checkpoint (see `CLAUDE.md`, "Traceability and resuming"). A new session reads the **Now** block first.

## Now

- **Stage:** 1 — Scaffold (not started)
- **Current task:** none yet
- **Next action:** ask the user decisions N1, N2, N10 and the template features (N3 is decided; N12 is waiting on the PANKEGG update); print the Stage 1 plan; write the Stage 1 task checklist below; wait for "go".
- **In flight / half-decided:** nothing
- **Last checkpoint:** project setup pushed to `origin/Upgrade` (2026-11-16)

## Stage status

| Stage | Status | Tag | Log | Notes |
|---|---|---|---|---|
| Project setup | DONE 2026-11-16 | — | — | No pipeline code changed |
| 0 Baseline | not started | | | optional, see N9 |
| 1 Scaffold | not started | | | |
| 2 nf-core modules + ext.args | not started | | | |
| 3 Local modules | not started | | | PANKEGG tax-genome test waits on N12 |
| 4 Subworkflows + workflow | not started | | | needs N4, N7, N11 |
| 5 Publishing | not started | | | |
| 6 Tests + CI | not started | | | |
| 7 Docs | not started | | | |
| 8 Comparison / Platform | not started | | | needs N8 |

## Task checklist — current stage

_Written at plan time, before any edit. Format: `- [ ] 1.1 <task> — <commit hash when done>`._

## Decisions

Report decisions (defaults stand until changed here):

| # | Decision | Status |
|---|---|---|
| D1 | Standalone nf-core-style pipeline `RVanDamme/muffin` | default |
| D2 | Samplesheet `--input` replaces `--illumina/--ont/--rna/--mode` (v3.0.0) | default |
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
| N1 | Minimum Nextflow version | OPEN (recommended `!>=26.04.0`) | |
| N2 | Licence GPL-3.0 vs MIT | OPEN (recommended keep GPL-3.0) | |
| N3 | Bin taxonomy method | **`--bin_taxonomy`: `sourmash` (gather + tax genome, default), `sourmash_lca` (legacy option), `gtdbtk` (higher accuracy, heavy resources, documented). `PANKEGG_MAKEDB` supports all three** (addendum §4a) | 2026-11-16 |
| N4 | All BAMs to COMEBin/SemiBin, drop merge | OPEN (recommended yes) | |
| N5 | Agents | **Claude does every stage and the review** | 2026-11-16 |
| N6 | PANKEGG | **Standalone PANKEGG, `pankegg_make_db` only; legacy `bin/` scripts not ported** | 2026-11-16 |
| N7 | GPU for SemiBin2/COMEBin | OPEN (recommended CPU only) | |
| N8 | Seqera Platform stage | OPEN (recommended skip) | |
| N9 | Legacy baseline run | OPEN | |
| N10 | Keep a `TEMPLATE` branch | OPEN (recommended yes) | |
| N11 | RNA results reporting without a PANKEGG RNA mode | OPEN (recommended publish + MultiQC) | |
| N12 | PANKEGG release reading `sourmash tax genome` output (tag + flag) | OPEN — user is updating PANKEGG; default method skips the PANKEGG step with a warning until then | |
| — | Push policy | **Push to `origin Upgrade` after every checkpoint commit; never force-push** | 2026-11-16 |

## Verification still owed

Items marked *(verify)* in `REPORT_ADDENDUM.md`: medaka model list / auto-selection; CheckM2 extension flag and `Name` column with gzipped bins; `SALMON_QUANT` with empty `gtf`; whether container configs must be regenerated after module installs; PANKEGG availability on Bioconda/PyPI; `sourmash lca classify` in sourmash 4.9.4; `tax genome` csv_summary columns; GTDB-Tk DB size / memory for the pinned release and its handling of gzipped bins; bin-name agreement between eggNOG prefix, CheckM2 `Name` and classification files.

## Session log

One line per session: date — what was done — last commit.

- 2026-11-16 — Project configured; report re-checked on Nextflow 26.10.6 / nf-core tools 4.2.2; corrections in `REPORT_ADDENDUM.md`. PANKEGG switched to the standalone tool; traceability protocol and reviewer added. Bin taxonomy set to three methods (N3).
