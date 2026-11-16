<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nf-test Assertions Reference

## Basic Assertions

```groovy
// Workflow/process success
assert workflow.success
assert process.success
assert function.success

// Failure testing
assert workflow.failed
assert workflow.exitStatus == 1
assert workflow.errorReport.contains("Missing required parameter")
```

## Path Assertions (nft-utils)

```groovy
def outputFile = path(process.out.txt[0][1])

// Existence
assert outputFile.exists()
assert !outputFile.isEmpty()

// Content
assert outputFile.text.contains("SUCCESS")
assert outputFile.readLines().size() == 10
assert outputFile.bytes.length > 0

// MD5 checksum
assert outputFile.md5 == "d41d8cd98f00b204e9800998ecf8427e"

// JSON/YAML
def json = outputFile.json
assert json.version == "1.0"

def yaml = outputFile.yaml
assert yaml.params != null

// Lines (for large files)
assert outputFile.lines.any { it.contains("PASS") }
assert outputFile.lines.count { it.startsWith("#") } == 5

// Gzipped files (for generic text — NOT for VCF/BAM/FASTA, use plugins instead)
assert path(process.out.fastq[0][1]).linesGzip.size() > 100
```

## Channel Assertions

```groovy
// Size
assert process.out.reads.size() == 1
assert workflow.out.bam.size() == 2

// Structure
assert process.out.reads[0][0] instanceof Map  // meta map
assert process.out.reads[0][0].id == "test"
assert process.out.reads[0][0].single_end == true

// File in channel
assert process.out.reads[0][1].name.endsWith(".fastq.gz")
assert path(process.out.reads[0][1]).exists()

// Iterate all items
process.out.reads.each { meta, reads ->
    assert meta.id != null
    assert reads != null
}

// Collect specific values
def ids = process.out.reads.collect { it[0].id }
assert ids.contains("sample1")
```

## Snapshot Assertions

```groovy
// Basic snapshot
assert snapshot(process.out).match()

// Named snapshot (multiple per test)
assert snapshot(process.out.bam).match("bam_output")
assert snapshot(process.out.bai).match("bai_output")

// Combined snapshot
assert snapshot(
    process.out.bam,
    process.out.bai,
    process.out.stats
).match("all_outputs")

// Selective snapshot (exclude unstable)
assert snapshot(
    process.out.versions,
    path(process.out.log[0][1]).readLines().findAll { !it.contains("timestamp") }
).match()

// File content snapshot
assert snapshot(
    path(process.out.txt[0][1]).text
).match("text_content")

// Sorted snapshot (for order-independent)
assert snapshot(
    process.out.ids.sort()
).match("sorted_ids")
```

## Version Assertions

nf-core modules emit versions:

```groovy
assert process.out.versions
assert snapshot(process.out.versions).match("versions")

// Parse versions.yml
def versions = path(process.out.versions[0]).yaml
assert versions.TOOL_NAME.tool == "1.2.3"
```

## assertAll Helper

Group assertions - all run even if some fail:

```groovy
then {
    assertAll(
        { assert process.success },
        { assert process.out.bam.size() == 1 },
        { assert process.out.bai.size() == 1 },
        { assert snapshot(process.out).match() }
    )
}
```

## Workflow-specific Assertions

```groovy
// Trace
assert workflow.trace.succeeded().size() == 5
assert workflow.trace.failed().size() == 0
assert workflow.trace.tasks().any { it.name == "FASTQC" }

// Stdout/Stderr
assert workflow.stdout.contains("Pipeline completed successfully")
assert !workflow.stderr.contains("ERROR")
```

## Process-specific Assertions

```groovy
// Exit status
assert process.exitStatus == 0

// Command executed
assert process.out.command.contains("--threads 4")

// Output directory
assert process.out.outdir.list().size() > 0
```

## File Type Assertions — Use Plugins

**ALWAYS use nf-test plugins for format-specific file assertions.** See plugins.md for full API docs.

```groovy
// VCF — use nft-vcf plugin (NEVER parse lines manually)
with(path(process.out.vcf[0][1]).vcf) {
    assert variantCount > 0
    assert sampleCount == 51
    assert chromosomes == ['20'] as Set
    assert variantsMD5 == "expected_md5"
}

// BAM — use nft-bam plugin (NEVER check magic bytes manually)
def bamFile = path(process.out.bam[0][1]).bam

// FASTQ — basic check (or use nft-fastq plugin)
def fastq = path(process.out.fastq[0][1])
assert fastq.linesGzip[0].startsWith("@")
```

## Custom Matchers

```groovy
// Define in test file or shared library
def matchesFastqFormat(path) {
    def lines = path.linesGzip.take(4)
    return lines[0].startsWith("@") &&
           lines[2] == "+" &&
           lines[1].length() == lines[3].length()
}

// Use in test
assert matchesFastqFormat(path(process.out.fastq[0][1]))
```
