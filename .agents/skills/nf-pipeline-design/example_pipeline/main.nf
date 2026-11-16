/*
 * Test pipeline for nflint — fast FASTA composition analysis.
 *
 * Follows the conventions in SKILL.md strictly:
 *   - main.nf stays thin and calls subworkflows only
 *   - every subworkflow executes one module per execution path
 *   - modules publish versions via typed topic sections (no versions.yml files)
 *
 * Running it end-to-end should take only a few seconds on a laptop.
 */

nextflow.enable.types = true

include { validateParameters } from 'plugin/nf-schema'
include { SPLIT_INPUT_FASTA       } from './subworkflows/local/split_input_fasta'
include { VALIDATE_SEQUENCES      } from './subworkflows/local/validate_sequences'
include { ANALYZE_GC              } from './subworkflows/local/analyze_gc'
include { ANALYZE_KMER            } from './subworkflows/local/analyze_kmer'
include { GENERATE_REPORT_SUMMARY } from './subworkflows/local/generate_report_summary'
include { RENDER_REPORT_DOC       } from './subworkflows/local/render_report_doc'

params {
    input: Path = file("${projectDir}/data/sample.fasta")
    outdir: String = 'results'
    kmer_size: Integer = 3
    summary_title: String = 'FASTA composition report'
}

workflow {

    // ----------------------------
    // Parameter validation and setup
    // ----------------------------
    validateParameters()
    ch_input_fasta = channel.of(params.input)
    kmer_size      = params.kmer_size
    report_title   = params.summary_title

    // ----------------------------
    // Pipeline run
    // ----------------------------

    /*
    Split the input FASTA into per-record files and attach a per-record meta map.
    */
    ch_records = SPLIT_INPUT_FASTA(ch_input_fasta)

    /*
    Filter each record to a clean ACGT alphabet before any composition maths.
    */
    ch_validated = VALIDATE_SEQUENCES(ch_records)

    /*
    Compute GC fraction per validated sequence.
    */
    ch_gc = ANALYZE_GC(ch_validated)

    /*
    Compute the k-mer frequency profile per validated sequence.
    */
    ch_kmer = ANALYZE_KMER(ch_validated, kmer_size)

    /*
    Collate per-sequence GC and k-mer TSVs into a single summary table.
    */
    ch_summary = GENERATE_REPORT_SUMMARY(
        ch_gc,
        ch_kmer,
    )

    /*
    Render the summary table into a human-readable markdown report.
    */
    RENDER_REPORT_DOC(
        ch_summary,
        report_title,
    )

    // ----------------------------
    // Software versions (topic channel)
    // ----------------------------
    // All modules publish [tool, version] tuples to the 'versions' topic
    // via their `topic:` section. Subscribe once here to collect them into
    // a single YAML report.
    channel.topic('versions')
        .unique()
        .map { tool, version ->
            "${tool}: ${version}"
        }
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            newLine: true,
            sort: true,
        )
}
