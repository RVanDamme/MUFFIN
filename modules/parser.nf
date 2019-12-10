process parser {
    label 'pankegg'
    publishDir "${params.output}/${name}/final_result", mode: 'copy', pattern: "*.html"
    publishDir "${params.output}/${name}/final_result", mode: 'copy', pattern: "*.csv"
    input:
        set val(name), file(rna_annot), file(quant)
        set val(name), val(bin_id), file(bins_annot)
    output:
        file("results/*.html") 
        file("results/*.csv") 
    script:
        """
        parser.py -b ${bin_annot} -l ${quant} -o result -r ${rna_annot}
        """
    }

    // """
    //     #!/usr/bin/python

    //     import PANKEGG
    //     import PANKEGG.parser
    //     from PANKEGG.parser import *
    //     import sys
    //     sys.argv = [sys.argv[0], '-b', '!{bin_annot}' , '-l' ,'!{quant}', '-o', 'result', '-r', '!{rna_annot}']
    //     main()

    //     """