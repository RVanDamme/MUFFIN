<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Nextflow DSL2 Reference for R Script Conversion

## Process syntax

```nextflow
process PROCESS_NAME {
    tag "$meta.id"
    label 'process_medium'

    conda 'r-base=4.3 r-tidyverse=2.0'
    container 'rocker/tidyverse:4.3'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.csv"), emit: results
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript - <<'RSCRIPT'
    library(tidyverse)

    data <- read_csv("${reads}")
    write_csv(data, "${prefix}.results.csv")

    writeLines(
        c("'${task.process}':", paste("  r-base:", as.character(getRversion()))),
        "versions.yml"
    )
    RSCRIPT
    """
}
```

## Workflow syntax

```nextflow
workflow {
    // Create channel from input
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> [ [id: row.sample], file(row.fastq) ] }

    // Chain processes
    PREPROCESS(ch_input)
    ANALYZE(PREPROCESS.out.results)
    PLOT(ANALYZE.out.results)
}
```

## Channel factories

```nextflow
// From file paths
Channel.fromPath(params.input)
Channel.fromFilePairs(params.reads)

// From values
Channel.of(1, 2, 3)
Channel.value("reference.fa")

// From CSV/TSV samplesheet
Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> [ [id: row.sample], file(row.file) ] }
```

## Channel operators

```nextflow
// Transform
.map { it -> [it.id, it.file] }
.transpose()                  // Unpack list elements (prefer over flatMap)
.filter { it -> it.size > 0 }

// Combine
ch_a.join(ch_b)              // Join by key (first element)
ch_a.combine(ch_b)           // Cross product
ch_a.mix(ch_b)               // Merge
ch_a.concat(ch_b)            // Concatenate in order

// Split/collect
.splitCsv(header: true)
.splitText(by: 1000)
.collect()                    // Gather all items
.groupTuple()                // Group by key
.toSortedList()

// Branch
ch.branch {
    pass: it.quality > 30
    fail: true
}
```

> **Tip:** When a process emits `[meta, [file1, file2, file3]]`, use `.transpose()`
> to unpack into `[meta, file1], [meta, file2], [meta, file3]`. It's cleaner than
> `.flatMap { meta, files -> files.collect { f -> [meta, f] } }` and automatically
> preserves metadata. Only use `.flatMap()` when you need to filter or transform
> elements during unpacking.

## Input qualifiers

```nextflow
input:
val(x)           // Simple value
path(file)       // File (staged into work dir)
tuple val(meta), path(reads)  // Tuple of mixed types
env(VAR_NAME)    // Environment variable
stdin             // Pipe to stdin
```

## Output qualifiers

```nextflow
output:
path "*.csv"                    // Glob pattern
path "output.txt", emit: result // Named emit
tuple val(meta), path("*.bam")  // Tuple
env(MY_VAR)                     // Capture env var
stdout                          // Capture stdout
```

## Process directives

```nextflow
process FOO {
    // Resources
    cpus 4
    memory '8 GB'
    time '2h'
    label 'process_high'

    // Execution
    conda 'r-base=4.3 bioconductor-deseq2=1.42'
    container 'quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_0'
    publishDir params.outdir, mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    // Caching
    cache 'lenient'
    storeDir params.storeDir
}
```

## Config (nextflow.config)

```nextflow
params {
    input     = null
    outdir    = null
    genome    = 'GRCh38'
}

process {
    withLabel: 'process_low' {
        cpus   = 2
        memory = '4 GB'
    }
    withLabel: 'process_medium' {
        cpus   = 4
        memory = '8 GB'
    }
    withLabel: 'process_high' {
        cpus   = 8
        memory = '32 GB'
    }
}

profiles {
    docker {
        docker.enabled = true
    }
    conda {
        conda.enabled = true
    }
}
```

## nf-core conventions

- Process names: `UPPERCASE_SNAKE_CASE`
- Use `meta` map for sample metadata: `[id: 'sample1', single_end: false]`
- Emit `versions.yml` from every process
- Use `task.ext.args` for optional CLI arguments
- Use `task.ext.prefix` for output file naming
