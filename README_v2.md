# MUFFIN v2 (in development)

> **Status page for the MUFFIN v2 rebuild.** This file is kept up to date at the end of every stage and will replace `README.md` when v2.0.0 is ready (Stage 7). The current `README.md` still describes MUFFIN v1, whose code now lives in `legacy/`.

MUFFIN is a metagenomics workflow: hybrid assembly, differential-coverage binning, bin quality control, taxonomic classification, functional annotation and optional metatranscriptomics. Version 2 is a rebuild on the [nf-core](https://nf-co.re) template (Nextflow ≥ 26.04, strict syntax, nf-core modules). It is an nf-core-*style* pipeline, not an official nf-core pipeline.

## Where the rebuild stands

| Stage | What | Status |
|---|---|---|
| 1 | Pipeline skeleton, parameters, samplesheets, test profile | done |
| 2 | nf-core modules installed and configured | not started |
| 3 | MUFFIN-specific modules (metaQUAST, sourmash tax, PANKEGG DB) | not started |
| 4 | Subworkflows and main workflow | not started |
| 5 | Output layout | not started |
| 6 | Tests and CI | not started |
| 7 | Documentation, this file becomes `README.md` | not started |
| 8 | Comparison with v1 results | not started |

**What works today:** nothing analytical yet. The skeleton validates parameters and samplesheets and has been exercised as a stub run with Apptainer on HPC2N Kebnekaise (MultiQC stubbed); no assembly, binning or annotation step is wired.

Details for developers: `docs/migration/PROGRESS.md` (current task, decisions), `docs/migration/log/` (what was run), `CHANGELOG.md` (what changes for users).

## How v2 will be run

```bash
nextflow run RVanDamme/MUFFIN -r Upgrade -profile <docker|apptainer|singularity> \
    --input samplesheet.csv --outdir results
```

Requirements: Nextflow ≥ 26.04 and Docker, Apptainer or Singularity.

On HPC2N Kebnekaise use `-profile hpc2n` (Slurm + Apptainer); setup and the `sbatch` launcher are in `docs/hpc2n.md`.

### Samplesheet (`--input`)

One row per sample. The columns present decide how the sample is processed; one run can mix sample types.

```csv
sample,short_reads_1,short_reads_2,long_reads,rna_reads_1,rna_reads_2
S1,s1_R1.fq.gz,s1_R2.fq.gz,s1_ont.fq.gz,s1_rna_R1.fq.gz,s1_rna_R2.fq.gz
S2,s2_R1.fq.gz,s2_R2.fq.gz,,,
S3,,,s3_ont.fq.gz,,
```

| Columns filled | Mode | Assembly |
|---|---|---|
| short pair + long reads | hybrid | metaSPAdes (hybrid) or metaFlye + polishing (`--hybrid_assembler`) |
| short pair only | short | metaSPAdes |
| long reads only | long | metaFlye + polishing |

RNA columns are optional and add the metatranscriptomics branch for that sample.

Two more optional samplesheets:

- `--coverage_input` (`sample,short_reads_1,short_reads_2,long_reads`): extra read sets mapped to every assembly for differential-coverage binning.
- `--input_bins` (`sample,bins_dir`): start from existing bins (classification and annotation only).

## What changes compared with v1

| v1 | v2 |
|---|---|
| `--illumina`, `--ont`, `--rna`, `--mode` | `--input` samplesheet |
| `--extra_ill`, `--extra_ont` | `--coverage_input` samplesheet |
| `--bin_classify`, `--bin_annotate` | `--input_bins` samplesheet |
| `--modular` | `--skip_assembly`, `--skip_classification`, `--skip_annotation` |
| `--assembler` | `--hybrid_assembler` |
| `--bintool` | `--binner` |
| `--output` | `--outdir` |
| `--cpus`, `--memory` | per-process resources, capped with `resourceLimits` in a config |

The full old → new table is in `docs/migration/MUFFIN_nfcore_migration_report.md` §5.3 and will move to `docs/usage.md`.

Bin taxonomy has three methods (`--bin_taxonomy`): `sourmash` (gather + tax genome, default, light), `sourmash_lca` (v1 behaviour, deprecated database format) and `gtdbtk` (most accurate, needs on the order of 100 GB of reference data and tens of GB of RAM).

The pathway report database is built with the standalone [PANKEGG](https://github.com/RVanDamme/PANKEGG) tool; open it with `pankegg_app`.

## Repository layout during the rebuild

| Path | Content |
|---|---|
| `main.nf`, `workflows/`, `subworkflows/`, `modules/`, `conf/`, `assets/` | v2 pipeline |
| `legacy/` | v1 code, same directory names as before (`legacy/modules/`, `legacy/bin/`, …); kept until the v1/v2 comparison is done |
| `docs/migration/` | plan, progress, logs and reports of the rebuild |
| `README.md` | v1 documentation (replaced by this file at release) |

## Known limitations right now

- No analysis steps yet (see the table above).
- Only a stub run has been done with containers so far; no real analysis run.
- The default taxonomy method cannot feed PANKEGG until PANKEGG reads `sourmash tax genome` output.
