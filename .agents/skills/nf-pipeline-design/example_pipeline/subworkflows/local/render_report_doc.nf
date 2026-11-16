//
// RENDER_REPORT_DOC: render the summary TSV into a human-readable markdown report.
//
include { RENDER_REPORT } from '../../modules/local/render_report'

workflow RENDER_REPORT_DOC {
    take:
    ch_summary       // path summary_tsv
    report_title     // val

    main:
    report_result = RENDER_REPORT(ch_summary, report_title)

    emit:
    report = report_result.report
}
