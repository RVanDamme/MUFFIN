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
| `EGGNOGMAPPER` (2.1.13) | `[meta, fasta]`, **`[search_mode, db]`**, `eggnog_data_dir` | `annotations`, `orthologs`, `hits` | `-m diamond` comes from the `search_mode` input, not `ext.args`; `--dbmem` is automatic above 40 GB; gz input handled |
| `TRINITY` (2.15.2) | `[meta, reads]` | `transcript_fasta` (`.fa.gz`) | labels `process_high` + `process_high_memory` |
| `SALMON_QUANT` (2.11.0) | `[meta, reads]`, `[meta2, index, gtf, transcript_fasta]` | `results` (dir) | no lib-type input → `-l A` goes in `ext.args`. Script always emits `--geneMap ${gtf}`; with `gtf = []` that leaves a bare `--geneMap` *(verify — likely needs a patch or a gene map)*. Uses `meta.single_end`. |

Not in nf-core/modules (confirmed): `metaquast`, `sourmash tax metagenome`, `sourmash tax genome`, an eggNOG DB downloader, pankegg.

### Consequences for the wiring

1. **Every assembly and every bin is gzipped** in the nf-core modules. Legacy MUFFIN used plain `.fa`. Bin-name matching in `BIN_QC` must strip `.fa.gz`/`.fa`/`.fna(.gz)`; CheckM2's extension option must match *(verify the exact flag and that the `Name` column drops the full extension)*; PILON needs `GUNZIP` first.
2. **`meta` needs `single_end`** for every short-read and RNA channel (`false` for pairs), because `FASTP` and `SALMON_QUANT` read it.
3. Report §7 `READ_MAPPING` emits `merged_bam`. With COMEBin taking a BAM directory and SemiBin accepting a BAM list, `SAMTOOLS_MERGE` is probably unnecessary (decision N4).
4. Report §6 row 24 offers `sourmash/taxannotate` for bin taxonomy. To replace `lca classify` like-for-like the right command is `sourmash tax genome`, which needs a local module (decision N3).
5. Medaka 2.x: `r941_min_high_g303` *(verify with `medaka tools list_models`)*. Medaka 2 can pick the model from the read headers, so `--medaka_model` defaulting to "auto" (no `-m`) is the robust choice *(verify)*.

## 4. PANKEGG

`bin/pankegg_bin.py`, `pankegg_bin_RNA.py`, `main.py`, `db_utils.py`, `sql_commands.py` import **only the Python standard library** (`sqlite3`, `csv`, `argparse`, `urllib`, …). No custom container build is needed: any pinned `python` biocontainer works.

`pankegg_bin.py` calls `http://rest.kegg.jp/get/<pathway>` at run time. A task that needs internet fails on offline HPC nodes and is not reproducible. Options: ship the lookups in `bin/pathway.txt`/`ko.txt`, or keep the call behind a flag (decision N6).

## 5. Execution environment

- Claude Code cloud sessions have **no Docker daemon** and cannot download nf-test. There you can scaffold, edit, `nextflow lint`, `nf-core pipelines lint` and install modules, but **not** run `-profile test,docker`, `-stub` with containers, or nf-test.
- Real runs need a machine with Docker **or** Apptainer/Singularity. Acceptance criteria use `-profile test,<engine>`; `<engine>` is whatever that machine has.
- Reachable from the cloud sandbox: `get.nextflow.io`, PyPI, `raw.githubusercontent.com` (nf-core test-datasets), `osf.io`. Not reachable: GitHub release assets.

## 6. Open decisions (answer before the stage that needs them)

| # | Needed by | Question | Recommended default |
|---|---|---|---|
| N1 | Stage 1 | Minimum Nextflow in `manifest.nextflowVersion` | `!>=26.04.0`, develop and test on 26.10.6 |
| N2 | Stage 1 | Licence: keep GPL-3.0 or move to MIT like nf-core | Keep GPL-3.0; skip the `license` template feature and keep the existing `LICENSE` |
| N3 | Stage 2/3 | Bin taxonomy method | Local `SOURMASH_TAXGENOME` after `SKETCH` → `GATHER` (closest to legacy). Alternative: GTDB-Tk (`nf-core/gtdbtk/classifywf`), the current standard but ~100 GB DB |
| N4 | Stage 4 | Feed COMEBin/SemiBin all BAMs instead of one merged BAM | Yes: keeps per-sample coverage, removes `SAMTOOLS_MERGE`. Behaviour change → CHANGELOG |
| N5 | Stage 1 | One agent or the three-agent split of report §10 | One agent (Claude Code) does all stages; Co-Scientist reviews. File-ownership rules dropped |
| N6 | Stage 3 | PANKEGG: keep the KEGG REST call at run time? Use the standalone PANKEGG package instead of `bin/` scripts? | Keep `bin/` scripts, make the network call optional |
| N7 | Stage 4 | GPU for SemiBin2/COMEBin | CPU only (do not request `accelerator`); document runtime |
| N8 | Stage 8 | Seqera Platform validation and Launchpad update | Optional; skip unless a Platform workspace is available |
| N9 | Stage 0 | Baseline run of legacy MUFFIN (report Stage 0) — where, and is it required? | Run locally if the legacy pipeline can be made to run; otherwise record "no baseline" and compare against expected ranges on the test data |
| N10 | Stage 1 | Keep a `TEMPLATE` branch so `nf-core pipelines sync` works later | Yes, push the scratch repo's `TEMPLATE` branch once |

Report decisions D1–D8 stand as written, with D8 read as "`publishDir` via `conf/modules.config`".
