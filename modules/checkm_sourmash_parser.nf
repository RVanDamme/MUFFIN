process sourmash_checkm_parser {
    label 'python38'
    publishDir "${params.output}/${name}/step2_final_result/", mode: 'copy', pattern: "step2_summary.csv"
    input:
    set val(name), file(checkm)
    file(sourmash)
    output:
    file("step2_summary.csv")
    shell:
    """
    grep -v "] INFO: " !{checkm} | grep -v "\\-\\-\\-\\-\\-\\-\\-" | grep -v "Bin Id" | sed -e 's/^[ \\t]*//'|sed 's/[ \\t]*\$//' |sed -r 's/ +/,/g'|sed '/^\$/d' >checkm.csv
    for file in !{sourmash}; do tail -n 1 \$file | sed -e 's/.fa//' >>sourmash.csv; done
    checkm_sourmash_parser.py -c checkm.csv -s sourmash.csv
    """
}