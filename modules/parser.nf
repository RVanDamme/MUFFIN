process parser_bin_RNA {
    label 'ubuntu'
    publishDir "${params.output}/${name}/annotate/", mode: 'copy', pattern: "parser_result/*"
    errorStrategy { task.exitStatus in 14..14 ? 'retry' : 'finish'}
    maxRetries 3 
    input:
        tuple val(name), path(rna_annot), path(quant)
        tuple val(name), path(bins_annot)
    output:
        path("parser_result/*") 
    script:
        """
        pankegg_bin_RNA.py -b ${bins_annot} -r ${rna_annot} -l ${quant} -o parser_result 
        """
    }
process parser_bin {
    label 'ubuntu'
    publishDir "${params.output}/${name}/annotate/", mode: 'copy', pattern: "parser_result/*"
    errorStrategy { task.exitStatus in 14..14 ? 'retry' : 'finish'}
    maxRetries 3 
    input:
        tuple val(name), path(bins_annot)
    output:
        path("parser_result/*") 
    script:
        """
        pankegg_bin.py -b ${bins_annot} -o parser_result 
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