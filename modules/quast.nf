process metaquast {
    maxForks 1
    label 'metaquast'
    publishDir "\${params.output}/\${name}/assemble/metaquast", mode: 'copy', pattern: "quality_dir"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.outdir}/metaquast_results", mode: 'copy'

    input:
    tuple val(name), path(assembly), path(reference_files)

    output:
    tuple val(name), path("quality_dir")

    script:
    """
    metaquast.py \\
        -o ${name}_metaquast_report \\
        -t \${task.cpus} \\
        -m \${params.min_contig} \\
        --max-ref-number \${params.max_ref_number} \\
        \${reference_files} \\
        \$assembly
    """
}
