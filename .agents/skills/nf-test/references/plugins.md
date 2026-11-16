<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nf-test Plugin Ecosystem

**ALWAYS prefer nf-test plugins over custom file parsing logic.** Plugins provide type-safe, validated assertions for common bioinformatics file formats. Writing manual parsing code (reading lines, checking magic bytes, splitting on tabs) is an anti-pattern when a plugin exists.

## Anti-Patterns — Do NOT Do This

Avoid these examples of custom assertion code:

```groovy
// ❌ WRONG — manual VCF line counting
assert path(process.out.vcf[0][1]).linesGzip.findAll { !it.startsWith("#") }.size() == 42

// ❌ WRONG — manual VCF header parsing
def vcfLines = path(process.out.vcf[0][1]).linesGzip
assert vcfLines[0].startsWith("##fileformat=VCF")
def headerLine = vcfLines.find { it.startsWith("#CHROM") }
def samples = headerLine.split("\t").drop(9)
assert samples.size() == 51

// ❌ WRONG — manual VCF variant parsing
def variants = path(process.out.vcf[0][1]).linesGzip.findAll { !it.startsWith("#") }
def chroms = variants.collect { it.split("\t")[0] }.unique() as Set
assert chroms == ["chr20"] as Set

// ❌ WRONG — manual BAM magic byte check
assert path(process.out.bam[0][1]).bytes[0..3] == [0x1f, 0x8b, 0x08, 0x04] as byte[]
```

**Instead, use the plugin equivalents below.**

## Plugin Setup

Add plugins to `nf-test.config`:

```groovy
config {
    plugins {
        load "nft-utils@0.0.9"
        load "nft-vcf@1.0.7"
        load "nft-bam@0.6.1"
    }
}
```

nf-test manages plugin downloading and caching. Use the installed version's configuration and reported cache details.

### CLI Override

Plugins can also be activated via the `--plugins` flag (colon-separated):

```bash
nf-test test --plugins nft-vcf@1.0.7:nft-bam@0.6.1
```

### CI / Parallel Execution

The plugin download mechanism is not safe for parallel execution. Pre-download plugins before running tests in parallel:

```bash
nf-test update-plugins   # download all plugins before parallel test runs
```

If plugins become corrupted or you need to force a fresh download:

```bash
nf-test clean             # clear cache and force redownload on next run
```

### Custom Repositories

```groovy
config {
    plugins {
        repository "https://github.com/my-org/nf-test-plugins/blob/main/plugins.json"
        load "my-custom-plugin@1.0.0"
        // or load a local JAR directly
        // loadFromFile userProvidedPluginJar
    }
}
```

Full plugin catalog: https://plugins.nf-test.com/
Official docs: https://www.nf-test.com/docs/plugins/using-plugins/

## nft-vcf (v1.0.7) — VCF File Assertions

Extends `path` with a `.vcf` property. Uses [HTSJDK](https://github.com/samtools/htsjdk) under the hood.

**Source:** https://github.com/seppinho/nft-vcf

### Basic Properties

```groovy
def vcfFile = path(process.out.vcf[0][1]).vcf

assert vcfFile.chromosomes == ['20'] as Set   // Set of chromosome names
assert vcfFile.chromosome == "20"             // First chromosome as String
assert vcfFile.sampleCount == 51              // Number of samples
assert vcfFile.phased                         // Whether VCF is phased
assert vcfFile.variantCount == 7824           // Total variant count
```

### with() Block Pattern (Preferred)

```groovy
with(path(process.out.vcf[0][1]).vcf) {
    assert chromosomes == ['20'] as Set
    assert sampleCount == 51
    assert phased
    assert variantCount == 7824
    assert chromosome == "20"
}
```

### Header (HTSJDK VCFHeader)

Returns the [VCFHeader](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html) instance.

```groovy
assert path(process.out.vcf[0][1]).vcf.header.getColumnCount() == 4
```

### Summary

Returns VCF summary attributes (chromosomes, variantCount, sampleCount, phasing status).

```groovy
path(process.out.vcf[0][1]).vcf.summary
```

### Variant Queries

```groovy
// Single variant by position (returns HTSJDK VariantContext)
def variant = path(process.out.vcf[0][1]).vcf.getVariant("chr20", 123)
assert variant.getContig() == "chr20"
assert variant.getHetCount() > 0
assert variant.getAttribute("AF") != null

// All variants as array of VariantContext
assert path(process.out.vcf[0][1]).vcf.variants.size() > 0

// First n variants
assert path(process.out.vcf[0][1]).vcf.getVariants(100).size() == 100

// Variants as strings (first n lines)
assert path(process.out.vcf[0][1]).vcf.getVariantsAsStrings(100).size() == 100

// MD5 of all variants (ideal for regression testing)
assert path(process.out.vcf[0][1]).vcf.variantsMD5 == "expected_md5"

// Variants by genomic range
assert path(process.out.vcf[0][1]).vcf.getVariantsByRange("chr20", 1, 10000).size() > 0
```

### INFO Field Queries

```groovy
// Get R2 INFO value at position
path(process.out.vcf[0][1]).vcf.getInfoR2("chr20", 1)

// Get any INFO tag value at position
path(process.out.vcf[0][1]).vcf.getInfoTag("R2", "chr20", 1)
```

### Indexing

```groovy
// Create tabix index (needed for range queries on unindexed VCFs)
path(process.out.vcf[0][1]).vcf.createIndex()
```

### Side-by-Side: Wrong vs Right

| Task | ❌ Custom Parsing | ✅ nft-vcf Plugin |
|------|------------------|-------------------|
| Count variants | `linesGzip.findAll { !it.startsWith("#") }.size()` | `vcf.variantCount` |
| Get samples | `headerLine.split("\t").drop(9).size()` | `vcf.sampleCount` |
| List chromosomes | `variants.collect { it.split("\t")[0] }.unique()` | `vcf.chromosomes` |
| Check phasing | manual GT field parsing | `vcf.phased` |
| Regression test | `md5` on full file (includes headers) | `vcf.variantsMD5` (variants only) |
| Query variant | line-by-line search + split | `vcf.getVariant("chr20", 123)` |
| Range query | filter + parse positions | `vcf.getVariantsByRange("chr20", 1, 10000)` |
| INFO field | split line, find field index | `vcf.getInfoTag("AF", "chr20", 123)` |

## nft-bam (v0.6.1) — BAM/CRAM File Assertions

Extends `path` with a `.bam` property.

**Documentation:** https://nvnieuwk.github.io/nft-bam

### Setup

```groovy
config {
    plugins {
        load "nft-bam@0.6.1"
    }
}
```

### Usage

```groovy
def bamFile = path(process.out.bam[0][1]).bam
```

See the [full nft-bam documentation](https://nvnieuwk.github.io/nft-bam) for all available properties and methods.

## nft-utils (v0.0.9) — Utility Functions & Path Extensions

Core utility functions for nf-test: path extensions (md5, json, yaml, lines), snapshot helpers, file collection, output sanitization, and nf-core dependency management.

**Source:** https://github.com/nf-core/nft-utils
**Docs:** https://nf-co.re/nft-utils

For the full API reference, see nft-utils.md.

## nft-fasta (v1.0.0) — FASTA File Assertions

```groovy
config {
    plugins {
        load "nft-fasta@1.0.0"
    }
}
```

## nft-fastq (v0.1.0) — FASTQ File Assertions

```groovy
config {
    plugins {
        load "nft-fastq@0.1.0"
    }
}
```

## Other Available Plugins

| Plugin | Description |
|--------|-------------|
| nft-compress | Compressed file assertions |
| nft-csv | CSV/TSV file assertions |
| nft-anndata | AnnData (h5ad) file assertions |
| nft-tiff | TIFF image file assertions |
| nft-parquet | Parquet file assertions |

Full catalog: https://plugins.nf-test.com/

## Best Practices

1. **Always use plugins** for format-specific assertions instead of manual parsing
2. **Load only needed plugins** in nf-test.config to keep test setup lean
3. **Use `with()` blocks** for multiple assertions on the same file
4. **Prefer `.variantsMD5`** for regression testing VCF content
5. **Use `.getVariant()`** for targeted variant validation instead of iterating all variants
6. **Never parse VCF/BAM/FASTA manually** when a plugin method exists — the plugin handles edge cases (compressed files, indices, malformed records) that custom parsing misses
