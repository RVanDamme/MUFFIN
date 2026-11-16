# Migration progress

Single source of truth for where the rebuild stands. Update it at the end of every stage and whenever a decision is taken.

## Stage status

| Stage | Status | Date | Notes |
|---|---|---|---|
| Project setup (CLAUDE.md, docs, skills) | DONE | 2026-11-16 | No pipeline code changed |
| 0 Baseline | not started | | optional, see N9 |
| 1 Scaffold | not started | | |
| 2 nf-core modules + ext.args | not started | | |
| 3 Local modules | not started | | |
| 4 Subworkflows + workflow | not started | | |
| 5 Publishing | not started | | |
| 6 Tests + CI | not started | | |
| 7 Docs | not started | | |
| 8 Comparison / Platform | not started | | |

## Decisions

Report decisions (defaults stand until changed here):

| # | Decision | Status |
|---|---|---|
| D1 | Standalone nf-core-style pipeline `RVanDamme/muffin` | default |
| D2 | Samplesheet `--input` replaces `--illumina/--ont/--rna/--mode` (v3.0.0) | default |
| D3 | Metric-level equivalence | default |
| D4 | User-supplied DB paths, opt-in downloads | default |
| D5 | Local modules first, upstream later | default |
| D6 | sourmash `lca classify` → `gather` + `tax` | default — needs explicit sign-off (method change) |
| D7 | One binner per run | default |
| D8 | Publishing as the template generates = `publishDir` | default |
| — | Short mode uses `spades --meta` (legacy ran isolate SPAdes) | default — needs explicit sign-off (behaviour change) |

Addendum decisions (see `REPORT_ADDENDUM.md` §6 for the options):

| # | Question | Answer |
|---|---|---|
| N1 | Minimum Nextflow version | OPEN (recommended `!>=26.04.0`) |
| N2 | Licence GPL-3.0 vs MIT | OPEN (recommended keep GPL-3.0) |
| N3 | Bin taxonomy method | OPEN (recommended local `SOURMASH_TAXGENOME`) |
| N4 | All BAMs to COMEBin/SemiBin, drop merge | OPEN (recommended yes) |
| N5 | One agent vs three-agent split | OPEN (recommended one agent) |
| N6 | PANKEGG network call / packaging | OPEN |
| N7 | GPU for SemiBin2/COMEBin | OPEN (recommended CPU only) |
| N8 | Seqera Platform stage | OPEN (recommended skip) |
| N9 | Legacy baseline run | OPEN |
| N10 | Keep a `TEMPLATE` branch | OPEN (recommended yes) |

## Verification still owed

Items marked *(verify)* in `REPORT_ADDENDUM.md`: medaka model list / auto-selection; CheckM2 extension flag and `Name` column with gzipped bins; `SALMON_QUANT` with empty `gtf`; whether container configs must be regenerated after module installs.

## Log

- 2026-11-16 — Project configured. Report claims re-checked on Nextflow 26.10.6 and nf-core/tools 4.2.2; corrections written to `REPORT_ADDENDUM.md`.
