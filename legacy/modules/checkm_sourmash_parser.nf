process sourmash_checkm_parser {
    //label 'python38'
    label 'ubuntu'
    publishDir "${params.output}/${name}/classify/", mode: 'copy', pattern: "classify_step_summary.csv"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(checkm)
    path(sourmash)
    output:
    path("classify_step_summary.csv")
    shell:
    """
    grep -v "] INFO: " !{checkm} | grep -v "\\-\\-\\-\\-\\-\\-\\-" | grep -v "Bin Id" | sed -e 's/^[ \\t]*//'|sed 's/[ \\t]*\$//' |sed -r 's/ +/,/g'|sed '/^\$/d' >checkm.csv
    for file in \$(ls bin*.txt); do tail -n 1 \$file | sed -e 's/.fa//' >>sourmash.csv; done
    checkm_sourmash_parser.py -c checkm.csv -s sourmash.csv
    """
}
