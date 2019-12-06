process parser {
    label 'parser'
    publishDir "${params.output}/${name}/final_result", mode: 'copy', pattern: "*.html"
    publishDir "${params.output}/${name}/final_result", mode: 'copy', pattern: "*.csv"
    input:
        set val(name), file(rna_annot), file(quant)
        set val(name), val(bin_id), file(bins_annot)
    output:
        file("results/*.html") 
        file("results/*.csv") 
    shell:
        """
        python -m ${annotation_parser} -r ${rna_annot} -l ${quant} -b ${bins_annot} -o results
        """
    }