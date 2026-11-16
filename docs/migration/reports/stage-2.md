## Stage 2 report

Status: DRAFT — review in progress
Environment: nextflow 26.10.6, nf-core 4.2.2, nf-test not installed, engine none (Claude Code cloud session)
Commits: 25c3fbd..HEAD of this report (range), tag `stage-2-done` after the review; log: docs/migration/log/stage-2.md
Files changed: `modules/nf-core/**` (30 modules via `nf-core modules install/update/remove`), `modules.json`, `conf/modules.config`, `nextflow_schema.json`, `CHANGELOG.md`, `README_v2.md`, `CLAUDE.md`, `docs/migration/*`.

Verification (command → result):
  - `nextflow lint . -exclude legacy -exclude .agents -exclude .claude` → 0 errors / 30 warnings (owned code: 0; all are `when:` deprecations in `modules/nf-core/**`).
  - `nf-core pipelines lint` → 325 passed / 26 warned / 0 failed (29 ignored).
  - `nf-core modules lint --all` → 1852 passed / 32 warned / 1 failed: multiqc `main_nf_version_topic`, by design upstream (MultiQC must not emit to the `versions` topic). Accepted exception to "clean".
  - `nf-core pipelines schema lint` → valid, 53 params.
  - runs → `nextflow run . -profile test -stub-run` without engine → `[SUCCESS] completed=1 failed=0` (config still loads; modules are not called yet).

Not verified (and why): every tool-option item of log 2.7 is source-checked, not executed (no container engine). To confirm in containers: COMEBin `-l` in 1.1.0; `sourmash lca classify` and `tax genome` csv_summary header in 4.9.4; GTDB-Tk `--extension gz`; CheckM2 `Name` for `.fa.gz` bins.

Behaviour changes (→ CHANGELOG): `--semibin_mode` is `pretrained|self` (semi-supervised dropped); Flye without `--plasmids`.

Deviations from report/addendum:
  - `cat/cat` → `find/concatenate` (deprecated upstream, `assert false` in its script) — N22.
  - chopper options in `ext.args2`, not `ext.args` (`modules/nf-core/chopper/main.nf`).
  - No patch of `salmon/quant`; transcript-to-gene map built in Stage 4 — N21.
  - `samtools/merge` not installed (N4 assumed yes; formally decided in Stage 4).
  - CheckM2 `Name` keeps `.fa` for gzipped bins → name matching must strip it (Stage 4).

### Module signatures (installed `main.nf`)

| Module | Real inputs (installed `main.nf`) | Emits | `ext.args` set | Difference from addendum §3 |
|---|---|---|---|---|
| `bwa/index` | `tuple val(meta), path(fasta)`; `val index_dir_name` | index | — | extra `val index_dir_name` (not in addendum) |
| `bwa/mem` | `tuple val(meta) , path(reads)`; `tuple val(meta2), path(index)`; `tuple val(meta3), path(fasta)`; `val sort_bam` | bam, cram, sam, csi, crai | — | as addendum |
| `checkm2/databasedownload` | `val(db_zenodo_id)` | database | — | as addendum |
| `checkm2/predict` | `tuple val(meta), path(fasta, stageAs: "input_bins/*")`; `tuple val(dbmeta), path(db)` | checkm2_output, checkm2_tsv | `--lowmem` if set | `Name` keeps `.fa` for `.fa.gz` bins (source-checked) |
| `chopper` | `tuple val(meta), path(fastq, stageAs: 'input/*')`; `path fasta` | fastq | args2: `-q/-l/--headcrop 100` | chopper options are `ext.args2` (args/args3 = pigz) |
| `comebin/runcomebin` | `tuple val(meta), path(assembly), path(bam, stageAs: "bam/*")` | bins, tsv, log, embeddings, covembeddings, embedding_ids | `-l 0.07|0.15` from `meta.n50` | as addendum |
| `eggnogmapper` | `tuple val(meta), path(fasta)`; `tuple val(search_mode), path(db)`; `path(eggnog_data_dir)` | annotations, orthologs, hits | `--itype metagenome|CDS --translate --go_evidence non-electronic --target_orthologs all` | as addendum |
| `fastp` | `tuple val(meta), path(reads), path(adapter_fasta)`; `val discard_trimmed_pass`; `val save_trimmed_fail`; `val save_merged` | reads, json, html, log, reads_fail, reads_merged | prefix `<id>_rna` in TRANSCRIPTOMICS | as addendum |
| `find/concatenate` | `tuple val(meta), path(files_in, stageAs: 'to_concatenate/*', arity: '1..*')` | file_out | — | replaces deprecated `cat/cat` |
| `flye` | `tuple val(meta), path(reads)`; `val mode` | fasta, gfa, gv, txt, log, json | `--meta` | as addendum |
| `gtdbtk/classifywf` | `tuple val(meta) , path("bins/*")`; `tuple val(db_name), path(db)`; `val use_pplacer_scratch_dir` | gtdb_outdir, summary, tree, markers, msa, user_msa, filtered, failed, log, warnings | `--extension gz` | `--extension gz` for gzipped bins (source-checked) |
| `gunzip` | `tuple val(meta), path(archive)` | gunzip | — | as addendum |
| `medaka` | `tuple val(meta), path(reads), path(assembly)` | assembly | `-m <model>` only if set | as addendum |
| `metabat2/jgisummarizebamcontigdepths` | `tuple val(meta), path(bam), path(bai)` | depth | — | as addendum |
| `metabat2/metabat2` | `tuple val(meta), path(fasta, stageAs: 'input_bins/'), path(depth)` | tooshort, lowdepth, unbinned, membership, fasta | — | extra `membership` output |
| `minimap2/align` | `tuple val(meta), path(reads)`; `tuple val(meta2), path(reference)`; `val bam_format`; `val bam_index_extension`; `val cigar_paf_format`; `val cigar_bam` | paf, bam, index | `-x map-ont --split-prefix tmp_idx` | as addendum |
| `multiqc` | `tuple val(meta), path(multiqc_files, stageAs: "?/*"), path(multiqc_config, stageAs: "?/*"), path(multiqc_logo), path(replace_names), path(sample_names)` | report, data, plots | — | updated; versions not on topic by design |
| `pilon` | `tuple val(meta), path(fasta)`; `tuple val(meta2), path(bam), path(bai)`; `val pilon_mode` | improved_assembly, vcf, change_record, tracks_bed, tracks_wig | — | as addendum |
| `racon` | `tuple val(meta), path(reads), path(assembly), path(paf)` | improved_assembly | — | as addendum |
| `salmon/index` | `tuple val(meta), path(transcript_fasta), path(genome_fasta)` | index | — | `[meta, transcript_fasta, genome_fasta]`; genome `[]` |
| `salmon/quant` | `tuple val(meta), path(reads)`; `tuple val(meta2), path(index), path(gtf), path(transcript_fasta)` | results, json_info, lib_format_counts | `--libType A` | always `--geneMap ${gtf}` → gene map built in Stage 4 (N21) |
| `samtools/fastq` | `tuple val(meta), path(input)`; `val interleave` | fastq, interleaved, singleton, other | `-f 4` (unbinned reads) | as addendum |
| `samtools/index` | `tuple val(meta), path(input)` | index | — | as addendum |
| `semibin/singleeasybin` | `tuple val(meta), path(fasta), path(bam)` | csv, model, tsv, output_fasta | args2: `--environment <env>` or `--self-supervised` | subcommand options in `ext.args2` |
| `seqkit/stats` | `tuple val(meta), path(reads)` | stats | — | as addendum |
| `sourmash/gather` | `tuple val(meta), path(signature)`; `path database`; `val save_unassigned`; `val save_matches_sig`; `val save_prefetch`; `val save_prefetch_csv` | result, unassigned, matches, prefetch, prefetchcsv | `-k 31` | as addendum |
| `sourmash/sketch` | `tuple val(meta), path(library, stageAs: 'library/*')`; `val merge_sigs` | signatures | `dna -p k=31,scaled=10000` (ONT 5000, reads with `abund`) | default args must be overridden |
| `spades` | `tuple val(meta), path(illumina), path(pacbio), path(nanopore)`; `path yml`; `path hmm` | scaffolds, contigs, transcripts, gene_clusters, gfa, warnings, log | `--meta` | as addendum |
| `trinity` | `tuple val(meta), path(reads, stageAs: "input*/*", arity: '1..*')` | transcript_fasta, log | — | as addendum |
| `untar` | `tuple val(meta), path(archive)` | untar | — | as addendum |
