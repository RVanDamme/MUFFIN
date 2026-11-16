# MUFFIN → nf-core-style pipeline on Nextflow 26.04+ — Migration Report

| | |
|---|---|
| Repository | https://github.com/RVanDamme/MUFFIN |
| Audited commit | `7e9447e7de85a8e11a09e70ebed824e06ef35457` (2025-02-09) |
| Nextflow used for audit | 26.04.4 (strict syntax parser on by default) |
| nf-core/modules snapshot | `master`, cloned 2026-11-16 |
| Report date | 2026-11-16 |
| Status | Analysis only — no code has been changed |

---

## 0. Executive summary

MUFFIN in its current form **does not run on Nextflow ≥ 26.04** (150 lint errors / 29 warnings over 43 files; 32 of 35 active `.nf` files fail). Rather than patching the legacy code line-by-line, the agreed direction is to **rebuild MUFFIN on the nf-core template**, replacing ~23 of the ~30 custom processes with nf-core modules, and writing the remaining MUFFIN-specific logic as nf-core-compliant local modules.

Estimated effort: **14–20 focused developer-days** (AI-assisted). Static types / records are deferred until nf-core modules adopt them.

During the audit, several **latent functional bugs** in the current pipeline were found (section 3). Some of them mean the *current* pipeline does not do what the README says (e.g. differential-coverage binning with extra reads never actually happens). The rebuild must implement the **intended** behaviour, and the baseline comparison must take these bugs into account.

---

## 1. Decisions (defaults used by this plan — confirm or override before starting)

| # | Decision | Default assumed in this plan |
|---|---|---|
| D1 | Standalone nf-core-*style* MUFFIN vs. contributing MUFFIN-specific parts to nf-core/mag | **Standalone nf-core-style pipeline** (`RVanDamme/muffin`), built with the nf-core template but not submitted as an official nf-core pipeline |
| D2 | Input interface | **Samplesheet (`--input`) replaces `--illumina/--ont/--rna` directories and `--mode`** (breaking CLI change, documented in CHANGELOG as v3.0.0) |
| D3 | Result equivalence | **Metric-level equivalence** vs. baseline (tool versions change, byte-identical outputs are not expected) |
| D4 | Databases | **User-supplied paths preferred**, with opt-in download (`--download_checkm2_db`, `--download_sourmash_db`, `--download_eggnog_db`) |
| D5 | Upstream contributions | **Local modules first**, written so they can later be contributed to nf-core/modules (metaquast, sourmash tax metagenome, pankegg) |
| D6 | Bin classification method | sourmash `lca classify` (LCA JSON DB, deprecated format) → **`sourmash gather` + `sourmash tax`** with the GTDB zip + lineages DB (method change — must be accepted by the user) |
| D7 | Binner selection | Keep **one binner per run** (`--binner metabat2|semibin2|comebin`) for parity; multi-binner + DAS Tool/Binette refinement is a later feature |
| D8 | Output publishing | **Follow whatever the current nf-core template generates** (do not hand-roll publishing). |

---

## 2. Current-state audit

### 2.1 Lint summary (Nextflow 26.04.4, `nextflow lint .`)

| Category | Count | Where |
|---|---|---|
| `Invalid process directive` (`errorStrategy = {...}`, `maxRetries = 5`) | 116 | every module |
| Conditional `include` / `include ... params(...)` | 5 files | `main.nf`, `readflows/*.nf`, `modules/modules_inclusion.nf` |
| Top-level statements mixed with declarations | 1 file | `main.nf` (log.info, version check, help, profile regex, `params.x = ...`) |
| `shell:` block deprecated | 19+ | `seqtk_retrieve_reads`, `sourmash`, `eggnog`, `polish`, `flye`, `ont_qc`, `trinity_and_salmon`, ... |
| `switch` statements | 4 | readflows + `main.nf` |
| Empty `input:` section → invalid process | 2 | `readme_output.nf`, `test_data_dll.nf` |
| Undefined / duplicate variables | 9 | see 3.1 |
| `workflow.onComplete` at top level | 3 | readflows |
| `nextflow.enable.dsl=2`, `Channel.` | 4 / 13 | everywhere |
| Dead code | — | `obsolete_process/` (~100 errors, delete) |

Readflow files stop parsing at line 6, so **their real error count is unknown** — irrelevant for the rebuild, but relevant if anyone tries a minimal patch.

### 2.2 Config / platform issues

| Issue | Location | Impact |
|---|---|---|
| `google-lifesciences` executor removed in Nextflow 25.10 | `nextflow.config` profile `gcloud` | Profile is unusable |
| Profile regex in `main.nf` only accepts `<executor>,<engine>` | `main.nf` | Launchpad profile `test, local` exits immediately on any NF version |
| `test` profile is empty | `nextflow.config` | No CI / smoke test possible |
| Duplicate params (`modular`, `ont`, `illumina`, `checkm2db`) | `nextflow.config` | Confusing; later value wins |
| `channels = [...]` instead of `conda.channels` | conda profile | Silently ignored |
| **`configs/container.config` defines every label twice in one `process {}` block — the second (old `nanozoo/*`) definitions override the Wave ones** | `configs/container.config` | With `-profile docker` the pipeline actually uses e.g. sourmash **2.0.1**, medaka **1.0.3**, flye **2.8**, eggnog-mapper **2.0.1**, checkm2 `latest`. sourmash 2.0.1 has no `sketch`/`tax` subcommands → classification likely fails under docker |
| Resources: every label = `params.cpus` / `params.memory` | `configs/*.config` | No per-tool sizing; replaced by nf-core `process_*` labels |
| Processes write into `params.output` with `cp` | `bins_tools.nf:separateBins` (unused) | Breaks on object storage |
| `storeDir`/`publishDir` to hard-coded `gs://gcloud_storage/...` | `sourmashgetdatabase.nf`, `eggnog_get_databases.nf` | Breaks outside that GCP project |
| `params.db_path = '/home/user/databases'` default | `nextflow.config`, `main.nf` | CheckM2 DB written outside work dir — not portable |

---

## 3. Latent functional bugs found (must be fixed *by design* in the rebuild)

### 3.1 Lint-detected
| Bug | Location | Effect |
|---|---|---|
| `${assmebly}` typo | `modules/quast.nf:18` | `--quast` always fails |
| `N50` undefined in `comebin` (env input) / n50 module | `modules/comebin.nf:106,126` | Fragile; NF ≥ 26 rejects |
| `${bins}` used as command | `modules/metaphlan.nf:16` | MetaPhlAn process can never work (it is also never called) |
| `name` undefined | `modules/ext_check.nf:8` | — |
| `name` declared twice in input tuples | `parser.nf`, `polish.nf` (pilon, pilon2), `sourmash.nf` (sourmash_ill) | Second tuple silently overrides first |
| `metaquast` also publishes to undefined `params.outdir` | `modules/quast.nf` | — |

### 3.2 Logic bugs (found by reading the dataflow)
| Bug | Location | Effect |
|---|---|---|
| **Extra reads never mapped**: `extra_ont_ch` / `extra_ill_ch` emit bare paths, then `assembly_ch.join(extra_ch)` joins on the *path* as key → never matches | all 3 readflows | `--extra_ont` / `--extra_ill` (differential-coverage binning) silently does nothing; `metabat2_extra` never runs |
| **CheckM2 depends on `checkm_download_db.out` even when the DB is provided** | all 3 readflows (`checkm2(classify_ch, checkm_download_db.out)`) | If `--checkm2db` is set or DB exists, `checkm_download_db` never runs → `checkm2` never runs (empty channel) |
| `get_fasta_path` is an empty process used to glob `*.fa` | `bins_tools.nf` | Should be a channel operation |
| `get_wrong_bin` awk paths assume `.fa` extension; comebin/semibin outputs may be `.fa.gz`/`.fna` depending on version | `bins_tools.nf` | Bin filtering may select nothing |
| `sourmash_ill` takes R1 and R2 as two *separate* tuples, both keyed `name` | `sourmash.nf` | Classifies R1 and R2 independently (2× work, 2 kreports) |
| `minimap2_bin` sorts without `-@ task.cpus` | `minimap2.nf` | Performance only |
| pankegg report (`parser_bin`, `parser_bin_RNA`) is **commented out** in all readflows | readflows | No final HTML report is produced today |
| `sourmash_checkm_parser` parses **CheckM v1** text output, never called | `checkm_sourmash_parser.nf` | Dead code |
| `cat_all_bins`, `separateBins`, `pilon`, `pilong`, `metaphlan`, `test`, `ont/illumina_reads_retrieval`, `sourmash_genome_size`, `pankegg_db` | various | Never called — **do not port** |
| `--short_qc` (default `2000`) is passed to chopper as **minimum ONT read length** `-l` | `ont_qc.nf` | Misleading name → rename `--longread_min_length` |
| chopper `--headcrop 100` hard-coded | `ont_qc.nf` | Move to `ext.args` |
| flye `--plasmids` | `flye.nf` | Option is deprecated/ignored in Flye ≥ 2.9 → verify and drop |
| medaka model `r941_min_high_g303` | `polish.nf` | **Verify** it still exists in medaka 2.x shipped by nf-core; otherwise expose `--medaka_model` with a current default |
| eggnog `--dbmem` always on | `eggnog.nf` | nf-core module enables `--dbmem` only when `task.memory > 40 GB` — acceptable, document |

---

## 4. Current dataflow (what the pipeline actually does today)

```
SHORT  : fastp → metaSPAdes(short) → bwa → [metaquast] → binner(asm, ill.bam)
LONG   : chopper → metaFlye → minimap2(paf) → racon → medaka → minimap2(bam) → binner(asm, ont.bam)
HYBRID : chopper + fastp →
           metaspades:  metaSPAdes(ill + ont)
           metaflye:    metaFlye → racon → medaka
         → minimap2(ont) + bwa(ill) → binner(asm, ont.bam + ill.bam; semibin/comebin get samtools-merged bam)

CLASSIFY (all modes):
  bins (from binner, or --bin_classify CSV)
  → checkm2 → bin_filter  (keep bin if Completeness − 5×Contamination > 50)
             get_wrong_bin (the rest, published only)
  → [unless --skip_bad_reads_recovery and no --bin_classify]:
       cat good bins → re-map reads (bwa / minimap2) → samtools view -f4 → seqtk subseq
       → sourmash sketch/gather/tax metagenome (kreport) on UNMAPPED reads (GTDB rs214 zip + lineages)
  → HYBRID+metaflye only: pilon × polish_iteration on each good bin (bwa + samtools + pilon loop)
  → per-bin: sourmash sketch + `lca classify` (GTDB rs214 LCA JSON)

ANNOTATE:
  bins (classified bins, or --bin_annotate / --bin_classify CSV)
  → per-bin eggNOG-mapper (diamond, --itype metagenome, --translate, non-electronic GO, all orthologs)
  → [--rna]: fastp_rna → Trinity → salmon (align_and_estimate_abundance.pl) → eggNOG (--itype CDS)
  → (pankegg HTML report — currently disabled)
```

---

## 5. Target architecture

### 5.1 Repository layout (nf-core template)

```
muffin/
├── main.nf                         # entry: PIPELINE_INITIALISATION → MUFFIN → PIPELINE_COMPLETION
├── workflows/muffin.nf
├── subworkflows/local/
│   ├── utils_nfcore_muffin_pipeline/   # template-generated
│   ├── preprocess_reads/               # FASTP (+ CHOPPER) per read type
│   ├── assembly/                       # metaSPAdes / metaFlye + LONGREAD_POLISH
│   ├── longread_polish/                # MINIMAP2_ALIGN(paf) → RACON → MEDAKA
│   ├── read_mapping/                   # BWA_INDEX/BWA_MEM, MINIMAP2_ALIGN, SAMTOOLS_MERGE
│   ├── binning/                        # METABAT2 / SEMIBIN / COMEBIN (+ SEQKIT_STATS for N50)
│   ├── bin_qc/                         # CHECKM2 + score-based filtering (channel ops)
│   ├── bin_polish_pilon/               # BWA + PILON rounds (hybrid metaflye only)
│   ├── unbinned_read_classification/   # re-map to good bins → unmapped reads → sourmash
│   ├── bin_taxonomy/                   # SOURMASH_SKETCH → GATHER → TAX
│   ├── bin_annotation/                 # EGGNOGMAPPER per bin
│   ├── transcriptomics/                # FASTP → TRINITY → SALMON_INDEX/QUANT → EGGNOGMAPPER
│   └── reporting/                      # PANKEGG + MULTIQC
├── modules/nf-core/…                   # installed with `nf-core modules install`, never hand-edited
├── modules/local/<tool>/{main.nf,meta.yml,environment.yml,tests/}
├── conf/{base.config,modules.config,test.config,test_full.config}
├── assets/{schema_input.json,schema_bins.json,samplesheet*.csv,multiqc_config.yml}
├── bin/                                 # pankegg + parser scripts (from current bin/)
├── docs/{usage.md,output.md}
├── tests/                               # pipeline-level nf-test
├── nextflow.config, nextflow_schema.json, modules.json, .nf-core.yml
└── CHANGELOG.md, CITATIONS.md, README.md
```

### 5.2 Samplesheet (`--input`, validated by nf-schema)

```csv
sample,short_reads_1,short_reads_2,long_reads,rna_reads_1,rna_reads_2
S1,s1_R1.fq.gz,s1_R2.fq.gz,s1_ont.fq.gz,s1_rna_R1.fq.gz,s1_rna_R2.fq.gz
S2,s2_R1.fq.gz,s2_R2.fq.gz,,,
S3,,,s3_ont.fq.gz,,
```

- Read type is **derived per sample** (`meta.mode = hybrid|short|long`), replacing `--mode`. A single run may mix modes.
- `meta` fields: `id`, `mode`, `has_rna`.
- Optional `--coverage_input` samplesheet (`sample,short_reads_1,short_reads_2,long_reads`) for *extra* coverage samples (replaces `--extra_ill/--extra_ont`); every row is mapped against **every** assembly (`combine`, not `join`) and the BAMs are added to the binner's depth input → fixes bug 3.2-1.
- Optional `--input_bins` samplesheet (`sample,bins_dir`) to start at classification/annotation (replaces `--bin_classify` / `--bin_annotate`).

### 5.3 Parameter mapping (old → new)

| Old | New | Notes |
|---|---|---|
| `--illumina`, `--ont`, `--rna`, `--mode` | `--input` | samplesheet |
| `--extra_ill`, `--extra_ont` | `--coverage_input` | samplesheet |
| `--bin_classify`, `--bin_annotate` | `--input_bins` | samplesheet |
| `--modular full/assemble/classify/annotate/…` | `--skip_assembly` (implied by `--input_bins`), `--skip_classification`, `--skip_annotation` | nf-core convention |
| `--assembler metaspades|metaflye` | `--hybrid_assembler metaspades|metaflye` | only used for hybrid samples; short → metaSPAdes, long → metaFlye |
| `--bintool` | `--binner metabat2|semibin2|comebin` | |
| `--bining_model`, `--environment` | `--semibin_mode pretrained|self|semi`, `--semibin_environment` | |
| `--skip_ill_qc`, `--skip_ont_qc` | `--skip_shortread_qc`, `--skip_longread_qc` | |
| `--short_qc` (2000), `--ont_min_qc` (10) | `--longread_min_length` (2000), `--longread_min_quality` (10) | chopper |
| `--model` | `--medaka_model` | verify default against medaka 2.x |
| `--polish_iteration`, `--skip_pilon` | `--pilon_rounds` (0–2, default 2), `--skip_pilon` | |
| `--checkm2db`, `--checkm2_low`, `--db_path`, `--db_file`, `--checkm2db_force_update` | `--checkm2_db`, `--checkm2_lowmem`, `--download_checkm2_db` | |
| `--sourmash_db`, `--sourmash_db_full`, `--sourmash_db_lineage` | `--sourmash_db` (zip), `--sourmash_lineages` (csv), `--download_sourmash_db` | LCA JSON dropped (D6) |
| `--eggnog_db` | `--eggnog_data_dir`, `--download_eggnog_db` | |
| *(hard-coded 50)* | `--bin_min_score` (default 50) | `Completeness − 5×Contamination > bin_min_score` |
| `--skip_bin_sorting` | `--skip_bin_filter` | |
| `--skip_bad_reads_recovery` | `--skip_unbinned_read_classification` | |
| `--quast` | `--run_quast` | |
| `--output` | `--outdir` | |
| `--cpus`, `--memory`, `--cachedir`, `--profile`, `--runinfodir`, `--reassembly`, `--reference`, `--help` | removed | nf-core labels / `resourceLimits` / nf-schema help |

---

## 6. Process-by-process mapping

`ext.args` shown are the **exact flags MUFFIN uses today**; they must be placed in `conf/modules.config` with `withName:` selectors scoped by subworkflow (e.g. `'.*:TRANSCRIPTOMICS:FASTP'`).

| # | Current process (file) | Current command (core flags) | Replacement | `ext.args` / wiring notes |
|---|---|---|---|---|
| 1 | `chopper` (ont_qc) | `chopper -q ${ont_min_qc} -l ${short_qc} --headcrop 100` | `nf-core/chopper` | `"-q ${params.longread_min_quality} -l ${params.longread_min_length} --headcrop 100"`; `fasta` input = `[]` |
| 2 | `fastp` (fastp) | `fastp -i R1 -I R2 -o -O` (defaults) | `nf-core/fastp` | no extra args; `adapter_fasta=[]`, `false,false,false` |
| 3 | `fastp_rna` (fastp) | same | `nf-core/fastp` (alias in TRANSCRIPTOMICS) | separate `withName` |
| 4 | `spades` (spades) | `spades.py --meta -1 -2 --nanopore` | `nf-core/spades` | `ext.args = '--meta'`; input `[meta, [r1,r2], [], ont]`, `yml=[]`, `hmm=[]`; output `contigs` (gz) |
| 5 | `spades_short` (spades) | `spades.py -1 -2` (**no `--meta`!**) | `nf-core/spades` | **Decision:** current short mode runs *isolate* SPAdes, not metaSPAdes, despite README. Recommend `--meta` (metagenomic), flag as behaviour change |
| 6 | `flye` (flye) | `flye --nano-raw --plasmids --meta` | `nf-core/flye` | `mode = '--nano-raw'`, `ext.args = '--meta'` (drop `--plasmids` after verifying) |
| 7 | `minimap_polish` (minimap2) | `minimap2 -x map-ont --split-prefix temp_index` → paf | `nf-core/minimap2/align` | `ext.args = '-x map-ont --split-prefix tmp_idx'`; `bam_format=false` |
| 8 | `racon` (polish) | `racon -t reads paf asm` | `nf-core/racon` | input `[meta, reads, assembly, paf]` |
| 9 | `medaka` (polish) | `medaka_consensus -m ${model}` | `nf-core/medaka` | `ext.args = "-m ${params.medaka_model}"` (verify model) |
| 10 | `minimap2`, `extra_minimap2`, `minimap2_bin` | `minimap2 -ax map-ont --split-prefix` → sort | `nf-core/minimap2/align` (3 aliases) | `ext.args = '-ax map-ont --split-prefix tmp_idx'`; `bam_format=true`, `bam_index_extension='bai'` |
| 11 | `bwa`, `extra_bwa`, `bwa_bin` | `bwa index -a bwtsw` + `bwa mem` → sort | `nf-core/bwa/index` + `nf-core/bwa/mem` (or `fastq_align_bwa` subworkflow) | `sort_bam=true`; then `samtools/index` |
| 12 | `bam_merger`, `bam_merger_extra` | `samtools merge` | `nf-core/samtools/merge` | used for semibin2/comebin (single BAM input); **not needed for metabat2** |
| 13 | `metabat2`, `metabat2_single`, `metabat2_extra` | `jgi_summarize_bam_contig_depths` + `metabat2 -i -a` | `nf-core/metabat2/jgisummarizebamcontigdepths` + `nf-core/metabat2/metabat2` | all BAMs (sample + coverage) in one depth call; `maxForks 1` dropped |
| 14 | `semibin2` | `SemiBin2 single_easy_bin --compression none` + `--environment X` or `--training-mode=self/semi` | `nf-core/semibin/singleeasybin` | `ext.args` built from `params.semibin_mode` / `params.semibin_environment` |
| 15 | `n50` + `comebin` | `n50 asm` → `run_comebin.sh -l (N50>10000 ? 0.07 : 0.15)` | `nf-core/seqkit/stats` (`-a`) + `nf-core/comebin/runcomebin` | parse N50 from seqkit TSV into `meta.n50` with `.map`; `ext.args = { "-l ${meta.n50 > 10000 ? 0.07 : 0.15}" }` |
| 16 | `metaquast` (quast) | `metaquast -t -o` | **local `METAQUAST`** (nf-core `quast` runs `quast.py` only) | container: quast (same package ships `metaquast.py`) |
| 17 | `checkm_download_db` | `checkm2 database --download --path` | `nf-core/checkm2/databasedownload` | only if `--download_checkm2_db`; fixes bug 3.2-2 |
| 18 | `checkm2` | `checkm2 predict --input bins/*.fa [--lowmem]` | `nf-core/checkm2/predict` | `ext.args = params.checkm2_lowmem ? '--lowmem' : ''`; add `-x fa` handling for bin extension |
| 19 | `bin_filter`, `get_wrong_bin`, `get_fasta_path`, `separateBins` | awk on `quality_report.tsv` | **channel operators** (no process) | `splitCsv(sep:'\t', header:true)` on CheckM2 report → `branch { good: Completeness − 5×Contamination > params.bin_min_score }` → join with per-bin files keyed by bin name |
| 20 | `bin_merger` | `cat bins > merged_bin_assembly.fa` | `nf-core/cat/cat` | |
| 21 | `unmapped_illumina_retrieve`, `unmapped_ont_retrieve` | `samtools view -f4` → ids → `seqtk subseq` | `nf-core/samtools/fastq` with `ext.args='-f 4'` (direct unmapped extraction) **or** `samtools/view` + `seqtk/subseq` | Recommend `SAMTOOLS_FASTQ` (simpler; near-equivalent — current code subsets the original FASTQ by unmapped read *names*, so a pair is kept if either mate is unmapped; reproduce with `-f 4` + mate handling or keep `seqtk/subseq` if exact parity is required) |
| 22 | `sourmash_download_db*` | `wget` GTDB rs214 | **local `SOURMASH_DOWNLOADDB`** or direct `file(url)` staging | only if `--download_sourmash_db`; update to current GTDB release |
| 23 | `sourmash_ill`, `sourmash_ont` | `sketch dna -p k=31,scaled=10000|5000` → `gather -k 31` → `tax metagenome --output-format kreport` | `nf-core/sourmash/sketch` + `nf-core/sourmash/gather` + **local `SOURMASH_TAXMETAGENOME`** | R1+R2 sketched together (fix 3.2-5); ONT `scaled=5000` via separate `withName` |
| 24 | `sourmash_bins` | `sketch dna k=31,scaled=10000` → `lca classify --db *.lca.json` | `nf-core/sourmash/sketch` + `nf-core/sourmash/gather` + `nf-core/sourmash/taxannotate` (or local `sourmash tax genome`) | **Method change D6** |
| 25 | `pilon2` | loop: bwa index/mem → sort → index → `pilon -Xmx(40% mem)` × N | `BWA_INDEX` + `BWA_MEM` + `SAMTOOLS_INDEX` + `nf-core/pilon` — **two explicit aliased rounds** (`PILON_ROUND1`, `PILON_ROUND2`) gated by `params.pilon_rounds` | strict syntax has no loops over processes; per-bin `transpose()` |
| 26 | `eggnog_download_db` | `download_eggnog_data.py -y` | **local `EGGNOGMAPPER_DOWNLOADDB`** | only if `--download_eggnog_db` |
| 27 | `eggnog_bin` | `emapper.py -m diamond --itype metagenome --translate --go_evidence non-electronic --target_orthologs all --dbmem` + strip header/footer | `nf-core/eggnogmapper` | `ext.args = '--itype metagenome --translate --go_evidence non-electronic --target_orthologs all'`; header stripping → done in the pankegg parser, not a process |
| 28 | `de_novo_transcript_and_quant` | Trinity → `align_and_estimate_abundance.pl --est_method salmon` | `nf-core/trinity` + `nf-core/salmon/index` + `nf-core/salmon/quant` | Trinity's salmon wrapper ≈ salmon index + quant; verify `-l A` |
| 29 | `eggnog_rna` | `emapper.py --itype CDS ...` | `nf-core/eggnogmapper` (alias) | `ext.args = '--itype CDS --translate --go_evidence non-electronic --target_orthologs all'` |
| 30 | `parser_bin`, `parser_bin_RNA` (bin/pankegg_*.py) | currently **disabled** | **local `PANKEGG`** module | not on Bioconda → build container with Wave from a pip/conda env; scripts stay in `bin/` or become a package |
| 31 | `readme_output`, `test`, `ext_check`, `help` | wget README / test data, file-ext check, help text | removed | `docs/output.md`, `conf/test.config`, nf-schema validation, `--help` |
| 32 | — | — | `nf-core/multiqc` | fastp, quast, samtools/flagstat (optional), salmon |
| — | `metaphlan`, `cat_all_bins`, `separateBins`, `pilon`, `pilong`, `ont/illumina_reads_retrieval`, `sourmash_genome_size`, `pankegg_db`, `sourmash_checkm_parser`, all of `obsolete_process/` | — | **not ported** | dead code |

All 20 nf-core modules above that were linted pass `nextflow lint` on 26.04.4 with **0 errors**, and already emit versions via the `versions` topic channel.

### 6.1 Local modules to write (nf-core style)

Each in `modules/local/<tool>/` with `main.nf`, `meta.yml`, `environment.yml`, `tests/main.nf.test` (+ snapshot), `tag "$meta.id"`, `process_*` label, `ext.args`/`ext.prefix`, `stub:` block, versions via `topic: versions`, Wave/biocontainers container pinned to an exact version.

| Module | Purpose | Inputs → outputs |
|---|---|---|
| `METAQUAST` | assembly QC (metagenome) | `[meta, contigs]` → `[meta, dir]`, `report.tsv` |
| `SOURMASH_TAXMETAGENOME` | kreport of unbinned reads | `[meta, gather_csv]`, `lineages` → `[meta, *.kreport.txt]` |
| `SOURMASH_DOWNLOADDB` *(optional)* | fetch GTDB zip + lineages | `val url` → `path db`, `path lineages` |
| `EGGNOGMAPPER_DOWNLOADDB` *(optional)* | `download_eggnog_data.py` | → `path eggnog_data_dir` |
| `PANKEGG` | KEGG pathway HTML/SQLite report from bins (± RNA) | `[meta, eggnog_annotations, checkm2_report, sourmash_results, (rna_annotations, quant)]` → html, db |

---

## 7. Subworkflow contracts (channel shapes)

| Subworkflow | take | emit |
|---|---|---|
| `PREPROCESS_READS` | `ch_short [meta, [r1,r2]]`, `ch_long [meta, fq]` | `short`, `long` (same shapes) |
| `ASSEMBLY` | `ch_short`, `ch_long` (both keyed by `meta.id`) | `contigs [meta, fasta.gz]` |
| `LONGREAD_POLISH` | `[meta, draft]`, `[meta, long]` | `polished [meta, fasta]` |
| `READ_MAPPING` | `contigs`, `ch_short`, `ch_long`, `ch_coverage` | `bams [meta, [bam...], [bai...]]`, `merged_bam [meta, bam]` |
| `BINNING` | `contigs`, `bams`, `merged_bam` | `bins [meta, [fa...]]` |
| `BIN_QC` | `bins`, `checkm2_db` | `good_bins [meta, [fa...]]`, `bad_bins`, `checkm2_report [meta, tsv]` |
| `BIN_POLISH_PILON` | `good_bins`, `ch_short` | `bins [meta, [fa...]]` (only `meta.mode=='hybrid' && hybrid_assembler=='metaflye'`) |
| `UNBINNED_READ_CLASSIFICATION` | `good_bins`, `ch_short`, `ch_long`, `sourmash_db`, `lineages` | `kreports [meta, txt]` |
| `BIN_TAXONOMY` | `bins`, `sourmash_db`, `lineages` | `taxonomy [meta, csv]` |
| `BIN_ANNOTATION` | `bins`, `eggnog_data_dir` | `annotations [meta, [tsv...]]` |
| `TRANSCRIPTOMICS` | `ch_rna`, `eggnog_data_dir` | `transcripts`, `quant`, `annotations` |
| `REPORTING` | all of the above, `ch_versions` (topic) | `multiqc_report`, `pankegg` |

Rules: `meta.id` is the join key everywhere; bins are always `[meta, List<Path>]` and are split per bin with `transpose()` only where a tool is per-bin (sourmash, eggnog, pilon); per-bin `meta` gets `bin_id` added. Optional DB/side inputs use `[]` not dummy files. No `params.*` reads inside modules — only inside workflows/`modules.config`.

---

## 8. Staged plan with acceptance criteria

| Stage | Work | Acceptance criteria | Est. |
|---|---|---|---|
| **0 Baseline** | Run current MUFFIN (Nextflow 25.04 or `NXF_SYNTAX_PARSER=v1`, `-profile local,docker` with the **Wave containers only** — fix the duplicate label override locally for the baseline) on the OSF subset (`https://osf.io/9xmh4/download`) in short, long, hybrid(metaspades), hybrid(metaflye). | `baseline/metrics.tsv` with: assembly total length, #contigs, N50; #bins; per-bin CheckM2 completeness/contamination; #good bins; sourmash top taxon per bin; #eggNOG annotated genes per bin; Trinity #transcripts. Known bugs from section 3 noted next to affected metrics. | 1–1.5 d |
| **1 Template** | `nf-core pipelines create` (name `muffin`, org `RVanDamme`, keep: nf-schema, nf-test, MultiQC, CI, igenomes **off**, fastqc **off**); samplesheet schemas; `conf/test.config` using hosted test data. | `nextflow lint .` 0 errors; `nf-core pipelines lint` passes template checks; `-profile test,docker -stub` runs | 1–1.5 d |
| **2 nf-core modules** | `nf-core modules install` for every module in section 6; `conf/modules.config` with the `ext.args` above. | `modules.json` complete; no file under `modules/nf-core/` hand-edited (`nf-core modules lint --all` clean) | 1.5–2 d |
| **3 Local modules** | Section 6.1 modules + Wave containers. | each has passing `nf-test` (real + stub), `nf-core modules lint` style checks pass | 2–3 d |
| **4 Subworkflows & workflow** | Section 7 wiring; bug fixes from section 3 by design. | `-profile test,docker` passes all 4 modes; channel shapes as section 7 | 4–6 d |
| **5 Publishing** | As generated by template (D8); output layout documented in `docs/output.md`. | every file listed in `docs/output.md` exists after test run | 0.5–1 d |
| **6 Tests & CI** | pipeline nf-test per mode (`tests/*.nf.test`), GitHub Actions: `nextflow lint`, `nf-core pipelines lint`, nf-test matrix. | CI green | 2–3 d |
| **7 Docs** | `README.md`, `docs/usage.md` (with old→new param table §5.3), `docs/output.md`, `CITATIONS.md`, `CHANGELOG.md` (v3.0.0). | `nf-core pipelines lint` 0 failures | 1–1.5 d |
| **8 Platform** | Run on Seqera Platform compute env `renaud-van-damme-compute-env-5993`; compare with baseline; update Launchpad pipeline *Muffin* (ID 198241367522843): new revision tag, profile `test,docker` (or no profile + samplesheet), schema-driven form. | metrics within tolerance (§9); Launchpad launch succeeds | 1 d |
| *Later* | Static types / records; multi-binner + refinement; official nf-core/mag contribution. | — | — |

**Total ≈ 14–20 days.**

---

## 9. Equivalence tolerances (baseline vs. new)

| Metric | Tolerance |
|---|---|
| Assembly total length / N50 | ± 5 % (same assembler; version differences) |
| # bins (MetaBAT2) | ± 10 % |
| # good bins (score > 50) | ± 1 bin on test data |
| CheckM2 completeness / contamination of matched bins | ± 5 points |
| sourmash genus of matched bins | identical genus for ≥ 90 % of bins (method change D6) |
| eggNOG annotated genes per bin | ± 5 % |
| Behaviour changes (expected differences) | extra-reads now actually used; CheckM2 runs with pre-supplied DB; short mode uses `--meta` (if accepted); pankegg report now produced |

---

## 10. Work split (Co-Scientist / Claude Code / Codex)

File ownership prevents merge conflicts. One branch per stage; merge before the next stage starts.

| Agent | Owns | Tasks |
|---|---|---|
| **Seqera Co-Scientist** | `assets/schema_*.json`, `nextflow_schema.json`, `conf/test*.config`, `baseline/` | Stage 0 baseline on Platform; samplesheet + schema; verify tool option compatibility (medaka model, flye flags); Wave container for PANKEGG/METAQUAST; per-stage Platform validation; Launchpad update (Stage 8) |
| **Claude Code** | `main.nf`, `workflows/`, `subworkflows/local/`, `conf/modules.config`, `nextflow.config` | Stage 1 scaffold; Stage 4 wiring; `ext.args` translation; Stage 5 publishing |
| **Codex** | `modules/local/**`, `modules/nf-core/**` (install only), `tests/`, `docs/`, `.github/` | Stage 2 bulk install; Stage 3 one PR per local module; Stage 6 tests/CI; Stage 7 docs |

Sequencing: (0 ‖ 1) → (2+3 by Codex ‖ 4 by Claude Code against `meta.yml` contracts) → Co-Scientist validation → (5 ‖ 6 ‖ 7) → 8.

---

## 11. Risks

| Risk | Mitigation |
|---|---|
| Baseline itself fails (old containers, sourmash 2.0.1, CheckM2 DB bug) | Use Wave container set only; pre-supply DB **and** work around 3.2-2 locally for the baseline; document what could not be baselined |
| Large DBs (CheckM2 ~3 GB, eggNOG ~50 GB, GTDB sourmash ~10 GB) slow down tests | Test profile uses stub runs for DB-heavy steps + tiny DB subsets; full DB only in `test_full` |
| medaka 2.x dropped `r941_min_high_g303` | expose `--medaka_model`, pick current default, note in CHANGELOG |
| COMEBin / SemiBin2 GPU containers | nf-core modules use CPU builds — confirm runtime on test data |
| Method change sourmash LCA → gather+tax | D6 sign-off; document in CHANGELOG |
| pankegg not packaged | Wave-built container from `bin/` scripts + `environment.yml` |
| nf-core template / tools version drift | pin `nf-core/tools` version in `.nf-core.yml`; record it in CHANGELOG |
