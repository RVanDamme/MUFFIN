process metaquast {
    maxForks 1
    label 'QUAST'

    publishDir "${params.output}/${name}/assemble/metaquast", mode: 'copy', pattern: "*metaquast"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.outdir}/metaquast_results", mode: 'copy'

    input:
    tuple val(name), path(assembly)

    output:
    tuple val(name), path("*metaquast")

    script:
    """
    metaquast -t ${task.cpus} -o ${name}_metaquast ${assmebly}
    """
}
