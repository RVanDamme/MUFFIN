nextflow.enable.types = true

process RENDER_REPORT {
    tag "composition_report"
    label 'process_single'

    container 'quay.io/nf-core/ubuntu:20.04'

    input:
    summary_tsv: Path
    report_title: String

    output:
    report: Path = file("composition_report.md")

    topic:
    tuple('awk', eval('awk --version | head -1')) >> 'versions'

    script:
    """
    {
        echo "# ${report_title}"
        echo
        echo "Generated: \$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo
        echo '## Composition table'
        echo
        awk 'NF==0 { print ""; next }
             /^sample\\t/ {
                 line = "|"; n = NF
                 for (i = 0; i < n; i++) line = line " --- |"
                 print "| " \$0 " |"
                 print line
                 next
             }
             { print "| " \$0 " |" }' ${summary_tsv} \
            | sed 's/\\t/ | /g'
    } > composition_report.md
    """
}
