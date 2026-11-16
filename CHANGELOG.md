# RVanDamme/muffin: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

MUFFIN v2 is a rebuild of the pipeline on the [nf-core](https://nf-co.re/) template (nf-core/tools 4.2.2) for Nextflow >= 26.04 with strict syntax. The v1 code is kept under `legacy/` during development.

### `Added`

- nf-core template: parameter validation and `--help` (nf-schema), MultiQC report, per-process resource labels, institutional config profiles.
- `--coverage_input` samplesheet for extra read sets used in differential-coverage binning.
- `--input_bins` samplesheet to start from existing bins.
- `--bin_taxonomy` with three methods: `sourmash` (gather + tax genome, default), `sourmash_lca` (v1 method) and `gtdbtk`; `--sourmash_lca_db`, `--gtdbtk_db`.
- `--bin_min_score` (default 50) for the bin filter `Completeness - 5 x Contamination`, hard-coded in v1.

### `Changed`

- **Breaking:** reads are given with a samplesheet (`--input`: `sample,short_reads_1,short_reads_2,long_reads,rna_reads_1,rna_reads_2`) instead of `--illumina`, `--ont`, `--rna` and `--mode`. The mode (short, long, hybrid) is derived per sample, so one run can mix modes.
- **Breaking:** parameters renamed: `--extra_ill`/`--extra_ont` → `--coverage_input`; `--bin_classify`/`--bin_annotate` → `--input_bins`; `--modular` → `--skip_assembly`, `--skip_classification`, `--skip_annotation`; `--assembler` → `--hybrid_assembler`; `--bintool` → `--binner`; `--bining_model`/`--environment` → `--semibin_mode`/`--semibin_environment` (`pretrained` or `self`; the semi-supervised mode is deprecated in SemiBin2 and was dropped); `--skip_ill_qc`/`--skip_ont_qc` → `--skip_shortread_qc`/`--skip_longread_qc`; `--short_qc`/`--ont_min_qc` → `--longread_min_length`/`--longread_min_quality`; `--model` → `--medaka_model`; `--polish_iteration` → `--pilon_rounds`; `--checkm2db`/`--checkm2_low` → `--checkm2_db`/`--checkm2_lowmem`; `--sourmash_db_lineage` → `--sourmash_lineages` (`--sourmash_db_full` dropped, `--sourmash_db` is the GTDB zip); `--eggnog_db` → `--eggnog_data_dir`; `--skip_bin_sorting` → `--skip_bin_filter`; `--skip_bad_reads_recovery` → `--skip_unbinned_read_classification`; `--quast` → `--run_quast`; `--output` → `--outdir`.
- Short-read samples are assembled with metaSPAdes (`spades --meta`); v1 ran SPAdes in isolate mode for short-only data despite documenting metaSPAdes.
- Medaka selects its model automatically from the reads unless `--medaka_model` is given (v1 always used `r941_min_high_g303`).
- Databases are user-supplied paths; downloads are opt-in (`--download_checkm2_db`, `--download_sourmash_db`, `--download_eggnog_db`).
- Minimum Nextflow version is 26.04.0.

### `Removed`

- `--cpus`, `--memory`, `--cachedir`, `--runinfodir`, `--reassembly`, `--reference`, `--db_path`, `--db_file`, `--checkm2db_force_update`, `--profile` and the `<executor>,<engine>` profile pairs: resources come from process labels and `resourceLimits`, engines from the standard nf-core profiles (`-profile docker`, `apptainer`, `singularity`, …). `--help` is now generated from the parameter schema.

### `Fixed`

- Flye is no longer called with `--plasmids`, which Flye >= 2.9 ignores.

### `Dependencies`

| Tool | Version |
| --- | --- |
| Nextflow | >= 26.04.0 |
| nf-core/tools template | 4.2.2 |
| nf-schema | 2.7.3 |

### `Deprecated`

- `--bin_taxonomy sourmash_lca`: kept for continuity with v1 results; the LCA database format is deprecated upstream.
