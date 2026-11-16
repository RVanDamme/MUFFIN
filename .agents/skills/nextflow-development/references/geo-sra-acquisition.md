<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# GEO/SRA Data Acquisition

Download raw sequencing data from NCBI GEO/SRA and prepare it for nf-core pipelines.

**Use this when:** Reanalyzing published datasets, validating findings, or comparing results against public cohorts.

## Table of Contents

- [Workflow Overview](#workflow-overview)
- [Step 1: Fetch Study Information](#step-1-fetch-study-information)
- [Step 2: Review Sample Groups](#step-2-review-sample-groups)
- [Step 3: Download FASTQ Files](#step-3-download-fastq-files)
- [Step 4: Generate Samplesheet](#step-4-generate-samplesheet)
- [Step 5: Run nf-core Pipeline](#step-5-run-nf-core-pipeline)
- [Supported Pipelines](#supported-pipelines)
- [Supported Organisms](#supported-organisms)
- [Complete Example](#complete-example)
- [Troubleshooting](#troubleshooting)

---

## Workflow Overview

Example: "Find differentially expressed genes in GSE309891 (drug-treated vs control)"

Inspect study metadata, identify the available data types, confirm the desired
sample subset and reference genome, then acquire the selected reads and prepare
a samplesheet for the chosen nf-core pipeline.

---

## Instructions

When assisting users with GEO/SRA data acquisition:

1. **Always fetch study info first** to show the user what data is available
2. **Ask for confirmation before downloading** - Present the sample groups and sizes, then ask the user which subset to download
3. **Suggest appropriate genome and pipeline** based on the organism and data type
4. **Return to main SKILL.md workflow** after data preparation is complete

Example confirmation question:
```
Question: "Which sample group would you like to download?"
Options:
  - "RNA-Seq:PAIRED (42 samples, ~87 GB)"
  - "RNA-Seq:SINGLE (7 samples, ~4.5 GB)"
  - "All samples (49 samples, ~92 GB)"
```

---

## Step 1: Fetch Study Information

Get metadata about a GEO study before downloading.

If available, use the bundled `sra_geo_fetch.py` helper to inspect study metadata for the user's accession. Use the host's available execution tools with user-provided files and locations.

**Example:**
If available, use the bundled `sra_geo_fetch.py` helper to inspect study metadata for the user's accession. Use the host's available execution tools with user-provided files and locations.

**Output includes:**
- Study title and summary
- Organism (with auto-suggested genome)
- Number of samples and runs
- Data types (RNA-Seq, ATAC-seq, etc.)
- Estimated download size
- Suggested nf-core pipeline

**Save info to JSON:**
If available, use the bundled `sra_geo_fetch.py` helper to inspect study metadata for the user's accession. Use the host's available execution tools with user-provided files and locations.

---

## Step 2: Review Sample Groups

View sample groups organized by data type and layout. This is useful for studies with mixed data types.

If available, use the bundled `sra_geo_fetch.py` helper to review sample groups. Use the host's available execution tools with user-provided files and locations.

**Example output:**
```
Sample Group          Count Layout     GSM Range                    Est. Size
--------------------------------------------------------------------------------
RNA-Seq                  42 PAIRED     GSM2879618...(42 samples)      87.4 GB
RNA-Seq                   7 SINGLE     GSM2976181-GSM2976187           4.5 GB
--------------------------------------------------------------------------------
TOTAL                    49                                           91.9 GB

Available groups for --subset option:
  1. "RNA-Seq:PAIRED" - 42 samples (~87.4 GB)
  2. "RNA-Seq:SINGLE" - 7 samples (~4.5 GB)
```

**List individual runs:**
If available, use the bundled `sra_geo_fetch.py` helper to list matching runs and any selected data-type filter. Use the host's available execution tools with user-provided files and locations.

**DECISION POINT:** Review the sample groups. Decide which subset to download if the study has multiple data types.

---

## Step 3: Download FASTQ Files

Download FASTQ files from ENA (faster than SRA).

If available, use the bundled `sra_geo_fetch.py` helper to download the confirmed subset to the user's selected destination. Use the host's available execution tools with user-provided files and locations.

**Options:**
- `-o, --output`: Output directory (required)
- `-i, --interactive`: Interactively select sample group to download
- `-s, --subset`: Filter by data type (e.g., "RNA-Seq:PAIRED")
- `-p, --parallel`: Parallel downloads (default: 4)
- `-t, --timeout`: Download timeout in seconds (default: 600)

### Interactive Mode (Recommended)

Use `-i` flag for interactive sample selection when the study has multiple data types:

If available, use the bundled `sra_geo_fetch.py` helper to download the confirmed subset to the user's selected destination. Use the host's available execution tools with user-provided files and locations.

**Interactive output:**
```
============================================================
  SELECT SAMPLE GROUP TO DOWNLOAD
============================================================

  [1] RNA-Seq (paired)
      Samples: 42
      GSM: GSM2879618...(42 samples)
      Size: ~87.4 GB

  [2] RNA-Seq (single)
      Samples: 7
      GSM: GSM2976181-GSM2976187
      Size: ~4.5 GB

  [0] Download ALL (49 samples)
------------------------------------------------------------

Enter selection (0-2):
```

### Direct Subset Selection

Alternatively, specify the subset directly:

If available, use the bundled `sra_geo_fetch.py` helper to download the confirmed subset to the user's selected destination. Use the host's available execution tools with user-provided files and locations.

**Note:** Downloads automatically skip existing files. Resume interrupted downloads by re-running the command.

---

## Step 4: Generate Samplesheet

Create a samplesheet compatible with nf-core pipelines.

If available, use the bundled `sra_geo_fetch.py` helper to generate a samplesheet using the selected downloaded read files. Use the host's available execution tools with user-provided files and locations.

**Options:**
- `-f, --fastq-dir`: Directory containing downloaded FASTQ files (required)
- `-o, --output`: Output samplesheet path (default: samplesheet.csv)
- `-p, --pipeline`: Target pipeline (auto-detected if not specified)

**Example:**
If available, use the bundled `sra_geo_fetch.py` helper to generate a samplesheet using the selected downloaded read files. Use the host's available execution tools with user-provided files and locations.

**Output:** The script will:
1. Create samplesheet in the format required by the target pipeline
2. Display suggested genome reference
3. Show suggested nf-core command

---

## Step 5: Run nf-core Pipeline

After generating the samplesheet, the script provides a suggested command.

**Example output:**
Choose the appropriate nf-core pipeline and pin its version. Supply the user's
selected samplesheet, output destination, genome, and supported execution profile.

**DECISION POINT:** Review and confirm:
1. Is the suggested pipeline correct?
2. Is the genome reference correct for your organism?
3. Do you need additional pipeline options?

Then return to the main SKILL.md workflow (Step 1: Environment Check) to proceed with pipeline execution.

---

## Supported Pipelines

The skill auto-detects appropriate pipelines based on library strategy. Pipelines marked with ★ are fully supported with configs, samplesheet generation, and documentation. Others are suggested but require manual setup following nf-core documentation.

| Library Strategy | Suggested Pipeline | Support |
|------------------|--------------------|---------|
| RNA-Seq          | nf-core/rnaseq     | ★ Full  |
| ATAC-seq         | nf-core/atacseq    | ★ Full  |
| WGS/WXS          | nf-core/sarek      | ★ Full  |
| ChIP-seq         | nf-core/chipseq    | Manual  |
| Bisulfite-Seq    | nf-core/methylseq  | Manual  |
| miRNA-Seq        | nf-core/smrnaseq   | Manual  |
| Amplicon         | nf-core/ampliseq   | Manual  |

---

## Supported Organisms

Common organisms with auto-suggested genomes:

| Organism | Genome | Notes |
|----------|--------|-------|
| Homo sapiens | GRCh38 | Human reference |
| Mus musculus | GRCm39 | Mouse reference |
| Saccharomyces cerevisiae | R64-1-1 | Yeast S288C |
| Drosophila melanogaster | BDGP6 | Fruit fly |
| Caenorhabditis elegans | WBcel235 | C. elegans |
| Danio rerio | GRCz11 | Zebrafish |
| Arabidopsis thaliana | TAIR10 | Arabidopsis |
| Rattus norvegicus | Rnor_6.0 | Rat |

If available, consult the bundled genome configuration for supported references.

---

## Complete Example

Reanalyze GSE110004 (yeast RNA-seq):

If available, use the bundled `sra_geo_fetch.py` helper to inspect study metadata for the user's accession, then review sample groups, then download the confirmed subset to the user's selected destination, then generate a samplesheet using the selected downloaded read files. Use the host's available execution tools with user-provided files and locations.

### Alternative: Non-interactive Download

If available, use the bundled `sra_geo_fetch.py` helper to review sample groups, then download the confirmed subset to the user's selected destination. Use the host's available execution tools with user-provided files and locations.

---

## Troubleshooting

### ENA Download Fails
If ENA downloads fail, the data may need to be fetched directly from SRA:

Install the SRA Toolkit through the user's supported package environment, retrieve
the selected run accession with `prefetch`, and convert it with `fasterq-dump` into
the user-selected destination.

### No SRA Runs Found
Some GEO datasets only have processed data, not raw sequencing reads. Check:
If available, use the bundled `sra_geo_fetch.py` helper to inspect study metadata for the user's accession. Use the host's available execution tools with user-provided files and locations.
If "Runs: 0", the dataset may not have raw data in SRA.

### SuperSeries Support
GEO SuperSeries (which contain multiple SubSeries) are automatically handled. The tool will:
1. Detect that a GEO ID is a SuperSeries
2. Find the linked BioProject accession
3. Fetch all SRA runs from the BioProject

Example: GSE110004 is a SuperSeries that links to BioProject PRJNA432544.

### Genome Not Recognized
If the organism is not in the genome mapping, manually specify the genome:
If available, use the bundled `manage_genomes.py` helper to list supported references. Use the host's available execution tools with user-provided files and locations.

---

## Requirements

- Python 3.8+
- `requests` library (optional but recommended)
- `pyyaml` library (optional, for genome config)
- Network access to NCBI and ENA

Install optional dependencies:
```bash
pip install requests pyyaml
```
