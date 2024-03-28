process metaquast {
    maxForks 1
    label 'QUAST'

    conda 'bioconda::QUAST=5.2.0'

    publishDir "\${params.output}/\${name}/assemble/metaquast", mode: 'copy', pattern: "metaquast"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.outdir}/metaquast_results", mode: 'copy'

    input:
    tuple val(name), path(assembly), path(reference_files)

    output:
    tuple val(name), path("metaquast")

    script:
    """
    metaquast \\
        -o ${name}_metaquast_report \\
        -t \${task.cpus} \\
        -m \${params.min_contig} \\
        --max-ref-number \${params.max_ref_number} \\
        \${reference_files} \\
        \$assembly
    """
}
