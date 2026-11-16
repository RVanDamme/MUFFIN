# Migration report addendum — verified facts and corrections

Verified on 2026-11-16 against the live toolchain and `nf-core/modules@c87ca93` (2026-11-14).
Where this file and `MUFFIN_nfcore_migration_report.md` disagree, **this file wins**.
Items marked *(verify)* were read from source but not executed; check them in the pinned container before relying on them.

## 1. Toolchain (what "latest" means today)

| Tool | Verified version | Consequence |
|---|---|---|
| Nextflow | **26.10.6** is the current stable (the report audited 26.04.4) | Develop and test on 26.10.x. `nextflow lint` on 26.10 also type-checks. |
| nf-core/tools | **4.2.2** | Template sets `nextflowVersion = '!>=25.10.4'`; raise it (decision N1). Strict syntax is mandatory for template ≥ 4.0. |
| nf-schema | 2.7.3 (pinned by the template) | — |
| nf-test | template loads `nft-utils@1.2.1` | — |
| Agent skills | `nextflow-io/agent-skills` is **archived**; successor is `seqeralabs/nextflow-plugin` (v0.2.2, 28 skills) | Use the successor. See `CLAUDE.md` for which skills apply. |

Confirmed from the report: legacy code gives **150 errors / 29 warnings in 43 files** under 26.10.6 too; all 30 nf-core modules listed in §6 exist; all of them emit versions through `topic: versions`.

New under 26.10: nf-core modules still contain `when:` blocks, which now raise a deprecation **warning** (seen in the fresh template for `fastqc` and `multiqc`). "0 warnings" is therefore only achievable for code we own, never for `modules/nf-core/**`.

## 2. What `nf-core pipelines create` (4.2.2) actually generates

- It creates a **new git repository** (branches `master`, `dev`, `TEMPLATE`) in its output directory. It cannot be pointed at the existing repo root. Generate into a scratch directory with `--template-yaml`, then copy the files (not `.git`) into the `Upgrade` branch.
- Template feature keys for `skip_features`: `github github_badges changelog license ci test_config igenomes modules multiqc fastqc nf_schema nf_core_configs is_nfcore seqera_platform gpu codespaces vscode code_linters citations documentation rocrate email`.
- Publishing is **`publishDir` in `conf/modules.config`** (plus `outputDir = params.outdir` for compatibility). So D8 means `publishDir`, not the `output {}` block. The `migrate-nextflow-code` skill also says workflow outputs should wait for static typing, which is deferred.
- It generates `conf/containers_*.config` files. Regenerate with `nf-core pipelines containers create-configs` after adding modules *(verify whether lint requires it)*.
- The template ships an **MIT** `LICENSE`; MUFFIN is **GPL-3.0** (decision N2).
- Entry workflow is named `RVANDAMME_MUFFIN`; `main.nf` still references igenomes/`params.fasta` unless `igenomes` is skipped.
- `nf-core modules create` has no "local mode" flag. Run inside the pipeline repo it writes to `modules/local/`. Useful flags: `--label`, `--meta`, `--conda-name`, `--conda-package-version`, `--empty-template`, `--author`.

## 3. Real module signatures (corrections to report §6 and §7)

Read from `modules/nf-core/<tool>/main.nf`. Re-read them after installing; this table is a snapshot.

| Module (pinned tool) | Real inputs | Outputs that matter | Differs from the report |
|---|---|---|---|
| `CHOPPER` (0.14.1) | `[meta, fastq]`, `fasta` | `fastq` | as reported |
| `FASTP` (1.4.0) | `[meta, reads, adapter_fasta]`, `discard_trimmed_pass`, `save_trimmed_fail`, `save_merged` | `reads`, `json`, `html`, `log` | **adapter FASTA is inside the tuple**: `[meta, reads, []]`. Uses `meta.single_end`. |
| `SPADES` (4.1.0) | `[meta, illumina, pacbio, nanopore]`, `yml`, `hmm` | `contigs` (`.fa.gz`), `scaffolds` | as reported |
| `FLYE` (2.9.5) | `[meta, reads]`, `mode` | `fasta` (`.fasta.gz`) | as reported |
| `MINIMAP2_ALIGN` (2.30) | `[meta, reads]`, `[meta2, reference]`, `bam_format`, `bam_index_extension`, `cigar_paf_format`, `cigar_bam` | `paf`, `bam`, `index` | as reported |
| `RACON` (1.4.20) | `[meta, reads, assembly, paf]` | `improved_assembly` (`.fasta.gz`) | output is gzipped |
| `MEDAKA` (2.2.2) | `[meta, reads, assembly]` | `assembly` (`.fa.gz`) | output is gzipped; runs `medaka_consensus $args` |
| `PILON` (1.24) | `[meta, fasta]`, `[meta2, bam, bai]`, `pilon_mode` | `improved_assembly` (`.fasta`) | passes `--genome $fasta` unchanged, so it **needs an uncompressed FASTA** → add `nf-core/gunzip` |
| `BWA_MEM` (0.7.19) | `[meta, reads]`, `[meta2, index]`, `[meta3, fasta]`, `sort_bam` | `bam` | as reported |
| `SAMTOOLS_MERGE` (1.25) | `[meta, bams, indexes]`, `[meta2, fasta, fai, gzi]`, `index_format` | `bam`, `index` | three inputs, not one |
| `SAMTOOLS_FASTQ` (1.25) | `[meta, bam]`, `interleave` | `fastq` (`_1/_2`), `other`, `singleton` | long reads come out on **`other`**, not `fastq` |
| `METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS` (2.17) | `[meta, bams, bais]` | `depth` (`.txt.gz`) | as reported |
| `METABAT2_METABAT2` (2.17) | `[meta, fasta, depth]` | `fasta` (bins, **`.fa.gz`**), `unbinned`, `tooshort`, `lowdepth` | bins are gzipped |
| `SEMIBIN_SINGLEEASYBIN` (2.4.1) | `[meta, fasta, bam]` | `output_fasta` (**`.fa.gz`**) | command is `SemiBin2 ${args} single_easy_bin … ${args2}` → `--environment` / `--training-mode` go in **`ext.args2`**. `--compression none` must not be set (the output glob expects `.fa.gz`). `bam` may be a list. Has label `process_gpu`. |
| `COMEBIN_RUNCOMEBIN` (1.1.0) | `[meta, assembly, bams]` (staged as `bam/*`) | `bins` (**`.fa.gz`**) | **takes a directory of BAMs**, so no merged BAM is needed (decision N4). Has label `process_gpu`. |
| `SEQKIT_STATS` (2.13.0) | `[meta, reads]` | `stats` (tsv) | default args already `--all` (includes N50) |
| `CHECKM2_PREDICT` (1.1.0) | `[meta, bins]`, **`[dbmeta, db]`** | `checkm2_tsv` = **`${prefix}_checkm2_report.tsv`**, `checkm2_output` (dir) | DB is a tuple; the report file is not called `quality_report.tsv` at the channel level |
| `CHECKM2_DATABASEDOWNLOAD` | `val db_zenodo_id` (null → default record) | `database` = `[meta, .dmnd]` | emits the tuple shape `CHECKM2_PREDICT` wants |
| `SOURMASH_SKETCH` (4.9.4) | `[meta, library]`, `merge_sigs` | `signatures` | extra `val` input |
| `SOURMASH_GATHER` (4.9.4) | `[meta, signature]`, `database`, `save_unassigned`, `save_matches_sig`, `save_prefetch`, `save_prefetch_csv` | `result` (`.csv.gz`) | six inputs |
| `SOURMASH_TAXANNOTATE` (4.9.4) | `[meta, gather_csv]`, `taxonomy` | `result` (`.with-lineages.csv.gz`) | only **appends lineages to gather rows**; it does not produce one classification per bin (decision N3) |
| `EGGNOGMAPPER` (2.1.13) | `[meta, fasta]`, **`[search_mode, db]`**, `eggnog_data_dir` | `annotations`, `orthologs`, `hits` | `-m diamond` comes from the `search_mode` input, not `ext.args`; `--dbmem` is automatic above 40 GB. **Correction 2026-11-16 (Stage 2 review):** gz input is NOT handled: the module tests `fasta.extension == '.gz'` (`modules/nf-core/eggnogmapper/main.nf:31`), which is never true, so run `GUNZIP` before `EGGNOGMAPPER` |
| `TRINITY` (2.15.2) | `[meta, reads]` | `transcript_fasta` (`.fa.gz`) | labels `process_high` + `process_high_memory` |
| `SALMON_QUANT` (2.11.0) | `[meta, reads]`, `[meta2, index, gtf, transcript_fasta]` | `results` (dir) | no lib-type input → `-l A` goes in `ext.args`. Script always emits `--geneMap ${gtf}`; with `gtf = []` that leaves a bare `--geneMap` *(verify — likely needs a patch or a gene map)*. Uses `meta.single_end`. |

Not in nf-core/modules (confirmed): `metaquast`, `sourmash tax metagenome`, `sourmash tax genome`, `sourmash lca classify`, an eggNOG DB downloader, `pankegg_make_db`. `gtdbtk/classifywf` exists.

### Consequences for the wiring

1. **Every assembly and every bin is gzipped** in the nf-core modules. Legacy MUFFIN used plain `.fa`. Bin-name matching in `BIN_QC` must strip `.fa.gz`/`.fa`/`.fna(.gz)`; CheckM2's extension option must match *(verify the exact flag and that the `Name` column drops the full extension)*; PILON needs `GUNZIP` first.
2. **`meta` needs `single_end`** for every short-read and RNA channel (`false` for pairs), because `FASTP` and `SALMON_QUANT` read it.
3. Report §7 `READ_MAPPING` emits `merged_bam`. With COMEBin taking a BAM directory and SemiBin accepting a BAM list, `SAMTOOLS_MERGE` is probably unnecessary (decision N4).
4. Report §6 row 24 offers `sourmash/taxannotate` for bin taxonomy. The closest `tax` equivalent of `lca classify` is `sourmash tax genome`, which has no nf-core module and whose output PANKEGG cannot read yet (§4a, N12).
5. Medaka 2.x: `r941_min_high_g303` *(verify with `medaka tools list_models`)*. Medaka 2 can pick the model from the read headers, so `--medaka_model` defaulting to "auto" (no `-m`) is the robust choice *(verify)*.

## 4. PANKEGG (decided 2026-11-16: use the standalone tool, DB generation only)

The legacy `bin/pankegg_bin.py`, `pankegg_bin_RNA.py`, `main.py`, `db_utils.py`, `sql_commands.py`, `checkm_sourmash_parser.py`, `ko.txt` and `pathway.txt` are **not ported**. They move to `legacy/` and are replaced by `pankegg_make_db` from https://github.com/RVanDamme/PANKEGG. The web app (`pankegg_app`) is not run by the pipeline; the pipeline publishes the SQLite DB and documents how to open it.

Read from the PANKEGG source at tag `2.0.0` (`0ef3aee`; `master` was `3c6e652` on 2026-11-01):

- Install: `pip install` from the repo (`setup.py`, version 2.0.0, console script `pankegg_make_db`). It is **not on Bioconda or PyPI** *(verify)*, so the local module needs a Wave-built container from a pinned environment that pip-installs the tagged release. Licence MIT.
- CLI: `pankegg_make_db -i <samples.csv> -o <name> --output_dir <dir> [--gtdbtk]` → `<dir>/<name>.db`.
- `samples.csv` header: `Sample name,Annotation_dir,classification_dir,Checkm2_dir`; the three path columns are **glob patterns**, one row per sample. The module must stage the files and write this CSV itself, with paths relative to the task directory.
- Annotations: eggNOG `*.annotations` files. **Correction (run 2026-11-16):** PANKEGG 2.0.0 skips only one leading `##` line and rejects raw emapper output (several `##` header lines plus a footer). Raw-file support is part of the PANKEGG update (N12); the pipeline must not add a stripping process. The **bin name is the first two dot-separated parts of the file name** (`<a>.<b>.…`), and must equal the `Name` column of the CheckM2 report and the bin ID in the classification file. Set `ext.prefix` for `EGGNOGMAPPER` accordingly and check the CheckM2 `Name` values for gzipped bins *(verify)*.
- CheckM2: needs columns `Name`, `Completeness`, `Contamination`; the nf-core `${prefix}_checkm2_report.tsv` is a copy of `quality_report.tsv`, so it fits.
- Classification, default mode: expects the **`sourmash lca classify` CSV** (header starts `ID,status,…`, one data line per file, `.fa` stripped from the ID). With `--gtdbtk`: GTDB-Tk summary files (`user_genome`, `classification`).
- No network access at run time; KEGG tables ship inside the package (`data/`).
- No RNA mode. The legacy RNA report has no equivalent (decision N11).

**Consequence.** PANKEGG 2.0.0 cannot read `sourmash tax genome` output, which is the pipeline's default taxonomy method. The user will add that support to PANKEGG (N12); see §4a.

## 4a. Bin taxonomy: three methods (decided 2026-11-16, replaces D6)

Selected with `--bin_taxonomy`. One method per run. All three feed `PANKEGG_MAKEDB`.

| `--bin_taxonomy` | Role | Processes | Databases (params) | Cost |
|---|---|---|---|---|
| `sourmash` (**default**) | Fast, light, current sourmash method | `SOURMASH_SKETCH` → `SOURMASH_GATHER` → local `SOURMASH_TAXGENOME` (per bin) | `--sourmash_db` (GTDB k=31 zip), `--sourmash_lineages` (csv) | Minutes per bin, a few GB RAM, DB around 10 GB |
| `sourmash_lca` | Legacy MUFFIN behaviour, kept for continuity with old results | `SOURMASH_SKETCH` → local `SOURMASH_LCACLASSIFY` (per bin) | `--sourmash_lca_db` (`*.lca.json[.gz]`) | Light, but the LCA JSON format is deprecated upstream and newer GTDB releases may not be published in it. Print a deprecation warning at start-up |
| `gtdbtk` | Highest accuracy and the community standard (marker-gene placement + ANI) | `nf-core/gtdbtk/classifywf` (2.7.2), once per sample on all bins | `--gtdbtk_db` (directory or tarball → `UNTAR`) | Reference data on the order of 100 GB on disk, tens of GB of RAM (`process_high_memory`), typically hours per sample. State the exact figures for the pinned release in `docs/usage.md` and the schema help text *(verify against the GTDB-Tk docs for the release used)* |

Notes for the wiring:
- `GTDBTK_CLASSIFYWF` inputs: `[meta, bins]` (staged as `bins/*`), `[db_name, db]`, `val use_pplacer_scratch_dir`; output `summary` = `${prefix}/classify/*.summary.tsv`. Check which genome extension it expects for gzipped bins (`--extension` in `ext.args`) *(verify)*.
- `sourmash lca classify` still has to exist in the pinned sourmash 4.9.4 *(verify with `sourmash lca classify --help` in the container; if it is gone, stop and report)*.
- `SOURMASH_TAXGENOME` runs `sourmash tax genome --gather-csv … --taxonomy … --output-format csv_summary` and emits the per-bin classification CSV. Expected columns in sourmash 4.9: `query_name, status, rank, fraction, lineage, query_md5, query_filename, f_weighted_at_rank, bp_match_at_rank, query_ani_at_rank` *(verify in the container and save one real file as a fixture)*. **This file is the contract for the PANKEGG update**; any change of output format must be agreed with the user first.
- Validation: the DB params required by the chosen method must be set, the others ignored; fail early with a clear message.
- `UNBINNED_READ_CLASSIFICATION` (sourmash `tax metagenome` on unmapped reads) is independent of `--bin_taxonomy` and still needs `--sourmash_db` / `--sourmash_lineages` unless skipped.

### PANKEGG_MAKEDB and the three methods

| `--bin_taxonomy` | Classification files passed | `pankegg_make_db` invocation | Works with PANKEGG 2.0.0? |
|---|---|---|---|
| `sourmash` | `tax genome` CSVs, one per bin | flag or auto-detection **to be defined by the PANKEGG update** | **No.** Needs the new PANKEGG release (N12) |
| `sourmash_lca` | `lca classify` CSVs, one per bin | default mode | Yes |
| `gtdbtk` | `*.summary.tsv` | `--gtdbtk` | Yes |

The module takes the method as a `val` input and maps it to the right arguments in the script block (no `params` inside the module). Until N12 is resolved: build and test the module for `sourmash_lca` and `gtdbtk`; for `sourmash`, implement the wiring, mark the real test as pending, and make the pipeline skip `PANKEGG_MAKEDB` with a visible warning rather than fail. Do not add a format-conversion process to fake the old format.

## 5. Execution environment

- Claude Code cloud sessions have **no Docker daemon** and cannot download nf-test. There you can scaffold, edit, `nextflow lint`, `nf-core pipelines lint` and install modules, but **not** run `-profile test,docker`, `-stub` with containers, or nf-test.
- Real runs need a machine with Docker **or** Apptainer/Singularity. Acceptance criteria use `-profile test,<engine>`; `<engine>` is whatever that machine has.
- Reachable from the cloud sandbox: `get.nextflow.io`, PyPI, `raw.githubusercontent.com` (nf-core test-datasets), `osf.io`. Not reachable: GitHub release assets.

## 6. Open decisions (answer before the stage that needs them)

| # | Needed by | Question | Recommended default |
|---|---|---|---|
| N1 | Stage 1 | Minimum Nextflow in `manifest.nextflowVersion` | `!>=26.04.0`, develop and test on 26.10.6 |
| N2 | Stage 1 | Licence: keep GPL-3.0 or move to MIT like nf-core | Keep GPL-3.0; skip the `license` template feature and keep the existing `LICENSE` |
| N3 | — | Bin taxonomy method | **DECIDED 2026-11-16:** three selectable methods via `--bin_taxonomy`: `sourmash` (gather + tax genome, default), `sourmash_lca` (legacy), `gtdbtk` (highest accuracy, heavy). See §4a |
| N4 | Stage 4 | Feed COMEBin/SemiBin all BAMs instead of one merged BAM | Yes: keeps per-sample coverage, removes `SAMTOOLS_MERGE`. Behaviour change → CHANGELOG |
| N5 | — | One agent or the three-agent split of report §10 | **DECIDED 2026-11-16:** Claude does every stage **and the review**. Report §10 and all file-ownership rules are void |
| N6 | — | PANKEGG source | **DECIDED 2026-11-16:** standalone PANKEGG (`pankegg_make_db` only), see §4. Legacy `bin/` scripts are not ported |
| N7 | Stage 4 | GPU for SemiBin2/COMEBin | CPU only (do not request `accelerator`); document runtime |
| N8 | Stage 8 | Seqera Platform validation and Launchpad update | Optional; skip unless a Platform workspace is available |
| N9 | Stage 0 | Baseline run of legacy MUFFIN (report Stage 0) — where, and is it required? | Run locally if the legacy pipeline can be made to run; otherwise record "no baseline" and compare against expected ranges on the test data |
| N10 | Stage 1 | Keep a `TEMPLATE` branch so `nf-core pipelines sync` works later | Yes, push the scratch repo's `TEMPLATE` branch once |
| N11 | Stage 4 | PANKEGG has no RNA mode: how are transcriptome annotations and quantification reported? | Publish Trinity / salmon / eggNOG outputs and summarise in MultiQC; no RNA in the PANKEGG DB. Behaviour change → CHANGELOG |
| N12 | Stage 3/4 | PANKEGG release that reads `sourmash tax genome` CSVs: tag to pin, and the flag (or auto-detection) that selects it. The same release must also accept raw eggNOG output, `.fa.gz` bin names, and install correctly without a checkout (2.0.0's `pip install` omits its data files) | **OPEN — waiting for the user's PANKEGG update.** Until then the default method runs but skips `PANKEGG_MAKEDB` with a warning |

Report decisions D1–D8 stand as written, except D6, which is replaced by §4a, with D8 read as "`publishDir` via `conf/modules.config`".
