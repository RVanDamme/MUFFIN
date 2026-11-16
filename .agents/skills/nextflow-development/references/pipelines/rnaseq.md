<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nf-core/rnaseq

**Version:** 3.22.2

**Official Documentation:** https://nf-co.re/rnaseq/3.22.2/
**GitHub:** https://github.com/nf-core/rnaseq

> **Note:** When updating to a new version, check the [releases page](https://github.com/nf-core/rnaseq/releases) for breaking changes and update the version in commands below.

## Contents
- [Test command](#test-command)
- [Samplesheet format](#samplesheet-format)
- [Parameters](#parameters)
- [Output files](#output-files)
- [Downstream analysis](#downstream-analysis)

## Test command

```bash
nextflow run nf-core/rnaseq -r 3.22.2 -profile test,docker
```

Expected: ~15 min, produces a MultiQC report.

## Samplesheet format

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,USER_CTRL1_R1_FQ_GZ,USER_CTRL1_R2_FQ_GZ,auto
CONTROL_REP2,USER_CTRL2_R1_FQ_GZ,USER_CTRL2_R2_FQ_GZ,auto
TREATMENT_REP1,USER_TREAT1_R1_FQ_GZ,USER_TREAT1_R2_FQ_GZ,auto
```

| Column | Required | Values |
|--------|----------|--------|
| sample | Yes | Alphanumeric, underscores allowed |
| fastq_1 | Yes | User-provided R1 read file reference |
| fastq_2 | No | User-provided R2 read file reference (empty for single-end) |
| strandedness | Yes | `auto`, `forward`, `reverse`, `unstranded` |

**Strandedness guide:**
- `auto`: Inferred from data (recommended)
- `forward`: TruSeq Stranded, dUTP protocols
- `reverse`: Ligation-based protocols
- `unstranded`: Non-stranded protocols

## Parameters

For each run, supply the user's selected samplesheet with `--input` and output
destination with `--outdir`. The fragments below show pipeline options without
prescribing input/output locations.

### Minimal run
```bash
nextflow run nf-core/rnaseq -r 3.22.2 -profile docker \
    --genome GRCh38
```

### Common parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--aligner` | `star_salmon` | Options: `star_salmon`, `star_rsem`, `hisat2` |
| `--genome` | - | `GRCh38`, `GRCh37`, `mm10`, `BDGP6` |
| `--pseudo_aligner` | - | Set to `salmon` for pseudo-alignment only |
| `--skip_trimming` | false | Skip adapter trimming |
| `--skip_alignment` | false | Pseudo-alignment only |

### Custom reference
For a custom reference, supply the user's genome FASTA with `--fasta` and
annotation with `--gtf`. Supply a compatible prebuilt STAR index with
`--star_index` if available; otherwise the pipeline can build one.

## Output files

Locate the count matrix, TPM matrix, alignment files, MultiQC report, and
execution reports from the selected run's actual output inventory.

**Key outputs:**
- `salmon.merged.gene_counts.tsv`: Input for DESeq2/edgeR
- `salmon.merged.gene_tpm.tsv`: Normalized expression

## Downstream analysis

```r
library(DESeq2)
counts <- read.delim("salmon.merged.gene_counts.tsv", row.names=1)
coldata <- data.frame(
    condition = factor(c("control", "control", "treatment", "treatment"))
)
dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = coldata,
    design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treatment", "control"))
```

## Troubleshooting

**STAR index fails**: Increase memory with `--max_memory '64.GB'` or provide pre-built `--star_index`.

**Low alignment rate**: Verify genome matches species; check FastQC for adapter contamination.

**Strandedness detection fails**: Specify explicitly with `--strandedness reverse`.

## More Information

- **Full parameter list:** https://nf-co.re/rnaseq/3.22.2/parameters/
- **Output documentation:** https://nf-co.re/rnaseq/3.22.2/docs/output/
- **Usage documentation:** https://nf-co.re/rnaseq/3.22.2/docs/usage/
