<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nf-core Test Data Reference

## Using nf-core Test Data

nf-core provides standardized test data via `params.test_data`. Available in nf-core pipelines/modules.

### Configuration

Ensure test profile loads test data config:

```groovy
// User-provided Nextflow test configuration
params {
    test_data_base = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data'
}

// Load test data paths
includeConfig 'https://raw.githubusercontent.com/nf-core/modules/master/tests/config/test_data.config'
```

### Common Test Data Paths

#### SARS-CoV-2 (Small, Fast Tests)

```groovy
// Illumina reads
params.test_data['sarscov2']['illumina']['test_1_fastq_gz']
params.test_data['sarscov2']['illumina']['test_2_fastq_gz']
params.test_data['sarscov2']['illumina']['test_single_end_bam']
params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam']
params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai']

// Genome
params.test_data['sarscov2']['genome']['genome_fasta']
params.test_data['sarscov2']['genome']['genome_fasta_fai']
params.test_data['sarscov2']['genome']['genome_dict']
params.test_data['sarscov2']['genome']['genome_gtf']
params.test_data['sarscov2']['genome']['genome_gff3']
params.test_data['sarscov2']['genome']['genome_sizes']

// VCF
params.test_data['sarscov2']['illumina']['test_vcf']
params.test_data['sarscov2']['illumina']['test_vcf_gz']
params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi']
```

#### Homo sapiens (Larger, Comprehensive Tests)

```groovy
// Illumina reads
params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz']
params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz']
params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam']
params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai']

// Genome (chr22 subset)
params.test_data['homo_sapiens']['genome']['genome_fasta']
params.test_data['homo_sapiens']['genome']['genome_fasta_fai']
params.test_data['homo_sapiens']['genome']['genome_dict']
params.test_data['homo_sapiens']['genome']['genome_gtf']

// 10x Genomics
params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_bam']
params.test_data['homo_sapiens']['10xgenomics']['cellranger']['test_10x_barcodes']
```

#### Other Organisms

```groovy
// E. coli
params.test_data['bacteroides_fragilis']['genome']['genome_fasta']
params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz']

// Yeast
params.test_data['saccharomyces_cerevisiae']['genome']['genome_fasta']
```

### Usage in Tests

```groovy
nextflow_process {
    name "Test FASTQC"
    script "main.nf"
    process "FASTQC"

    test("sarscov2 - single_end") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:true ],
                    file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
                ]
                """
            }
        }
        then {
            assert process.success
        }
    }

    test("sarscov2 - paired_end") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ],
                    [
                        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
                    ]
                ]
                """
            }
        }
        then {
            assert process.success
        }
    }
}
```

### User-Provided Test Data

Use the user's verified test input references in the module's actual input
contract. Confirm availability in the intended execution environment and use
checkIfExists when appropriate. Do not prescribe a fixture folder.

To derive smaller inputs, use format-aware tools on the user's selected data:
seqtk for deterministic FASTQ sampling, samtools for a suitable BAM interval,
and bcftools for a suitable VCF interval. Preserve pairing, headers, indexing,
and the biological context required by the test. Let the user choose where
derived inputs are stored.

### Test Data Guidelines

1. **Use sarscov2 for fast tests** - Small genome, quick execution
2. **Use homo_sapiens for comprehensive tests** - More realistic but slower
3. **Always use `checkIfExists: true`** - Fail fast on missing data
4. **Cache test data** - Configure caching in CI to avoid repeated downloads
5. **Version test data** - Use specific commits/tags for reproducibility
