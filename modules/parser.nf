process parser {
    label 'pankegg'
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "parser_result/*"
    input:
        set val(name), file(rna_annot), file(quant)
        set val(name), val(bin_id), file(bins_annot)
    output:
        file("parser_result/*") 
    script:
        """
        pankegg.py -b ${bins_annot} -r ${rna_annot} -l ${quant} -o parser_result 
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