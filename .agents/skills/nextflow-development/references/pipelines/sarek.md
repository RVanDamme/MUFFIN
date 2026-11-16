<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nf-core/sarek

**Version:** 3.7.1

**Official Documentation:** https://nf-co.re/sarek/3.7.1/
**GitHub:** https://github.com/nf-core/sarek

> **Note:** When updating to a new version, check the [releases page](https://github.com/nf-core/sarek/releases) for breaking changes and update the version in commands below.

## Contents
- [Test command](#test-command)
- [Samplesheet format](#samplesheet-format)
- [Variant calling modes](#variant-calling-modes)
- [Parameters](#parameters)
- [Output files](#output-files)

## Test command

```bash
nextflow run nf-core/sarek -r 3.7.1 -profile test,docker
```

Expected: ~20 min, creates aligned BAMs and variant calls.

## Samplesheet format

### From FASTQ
```csv
patient,sample,lane,fastq_1,fastq_2
patient1,tumor,L001,USER_TUMOR_L001_R1_FQ_GZ,USER_TUMOR_L001_R2_FQ_GZ
patient1,tumor,L002,USER_TUMOR_L002_R1_FQ_GZ,USER_TUMOR_L002_R2_FQ_GZ
patient1,normal,L001,USER_NORMAL_R1_FQ_GZ,USER_NORMAL_R2_FQ_GZ
```

### From BAM/CRAM
```csv
patient,sample,bam,bai
patient1,tumor,USER_TUMOR_BAM,USER_TUMOR_BAM_BAI
patient1,normal,USER_NORMAL_BAM,USER_NORMAL_BAM_BAI
```

### With tumor/normal status
```csv
patient,sample,lane,fastq_1,fastq_2,status
patient1,tumor,L001,tumor_R1.fq.gz,tumor_R2.fq.gz,1
patient1,normal,L001,normal_R1.fq.gz,normal_R2.fq.gz,0
```

`status`: 0 = normal, 1 = tumor

## Variant calling modes

For each run, supply the user's selected samplesheet with `--input` and output
destination with `--outdir`. The fragments below show pipeline options without
prescribing input/output locations.

### Germline (single sample)
```bash
nextflow run nf-core/sarek -r 3.7.1 -profile docker \
    --genome GRCh38 \
    --tools haplotypecaller,snpeff
```

### Somatic (tumor-normal pair)
```bash
nextflow run nf-core/sarek -r 3.7.1 -profile docker \
    --genome GRCh38 \
    --tools mutect2,strelka,snpeff
```

### WES (exome)
```bash
nextflow run nf-core/sarek -r 3.7.1 -profile docker \
    --genome GRCh38 \
    --wes \
    --tools haplotypecaller,snpeff
```

### Joint germline (cohort)
```bash
--tools haplotypecaller --joint_germline
```

## Parameters

### Available tools

**Germline callers:**
- `haplotypecaller`: GATK HaplotypeCaller
- `freebayes`: FreeBayes
- `deepvariant`: DeepVariant (GPU optional)
- `strelka`: Strelka2 germline

**Somatic callers:**
- `mutect2`: GATK Mutect2
- `strelka`: Strelka2 somatic
- `manta`: Structural variants

**CNV callers:**
- `ascat`: Copy number
- `controlfreec`: CNV detection
- `tiddit`: SV calling

**Annotation:**
- `snpeff`: Functional annotation
- `vep`: Variant Effect Predictor

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tools` | - | Comma-separated list of tools |
| `--genome` | - | `GRCh38`, `GRCh37` |
| `--wes` | false | Exome mode (requires `--intervals`) |
| `--intervals` | - | BED file for targeted regions |
| `--joint_germline` | false | Joint calling for cohorts |
| `--skip_bqsr` | false | Skip base quality recalibration |

## Output files

Locate analysis-ready alignments, germline or somatic variant calls,
annotations, and the MultiQC report from the selected run's actual output inventory.

## Troubleshooting

**BQSR fails**: Check known sites available for genome. Skip with `--skip_bqsr` for non-standard references.

**Mutect2 no variants**: Verify tumor/normal pairing in samplesheet (check `status` column).

**Out of memory**: `--max_memory '128.GB'` for WGS.

**DeepVariant GPU**: Ensure NVIDIA Docker runtime configured.

## More Information

- **Full parameter list:** https://nf-co.re/sarek/3.7.1/parameters/
- **Output documentation:** https://nf-co.re/sarek/3.7.1/docs/output/
- **Usage documentation:** https://nf-co.re/sarek/3.7.1/docs/usage/
