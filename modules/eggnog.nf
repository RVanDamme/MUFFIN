process eggnog_bin { 
  label 'eggnog' 

  conda 'bioconda::diamond anaconda::biopython bioconda::eggnog-mapper=2.0.1 '

  publishDir "${params.output}/${name}/annotate/bin_annotation/", mode: 'copy', pattern: "*.tsv"
  errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
  maxRetries = 5
  input:
    tuple val(name), path(bin), path(db)
  output:
    tuple val(name), path("*.annotations.tsv")
    path("*.seed_orthologs.tsv")
  shell:
    """
    emapper.py --data_dir ${db} -d bact -o ${name}  -m diamond --data_dir ${bin} --cpu ${task.cpus} --go_evidence non-electronic  --target_orthologs all --translate
    tac ${name}.emapper.annotations | sed "1,3d" | tac |sed "1,3d" > ${name}.annotations.tsv
    cp ${name}.emapper.seed_orthologs ${name}.seed_orthologs.tsv
    """
}
// process eggnog_bin { 
//   label 'eggnog' 

//   conda 'bioconda::diamond anaconda::biopython bioconda::eggnog-mapper=2.0.1 '

//   publishDir "${params.output}/${name}/annotate/bin_annotation/", mode: 'copy', pattern: "*.tsv"
//   errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//   maxRetries = 5
//   input:
//     tuple val(name), path(bin), path(db)
//   output:
//     tuple val(name), path("*.annotations.tsv")
//     path("*.seed_orthologs.tsv")
//   shell:
//     """
//     bin_id=\$(basename !{bin} | sed -r "s/\\.\\w+//2")
//     emapper.py --data_dir ${db} -d bact -o \$bin_id  -m diamond -i ${bin} --cpu ${task.cpus} --go_evidence non-electronic  --target_orthologs all --translate
//     tac \$bin_id.emapper.annotations | sed "1,3d" | tac |sed "1,3d" > \$bin_id.annotations.tsv
//     cp \$bin_id.emapper.seed_orthologs \$bin_id.seed_orthologs.tsv
//     """
// }

process eggnog_rna { 
  label 'eggnog' 

  conda 'bioconda::diamond anaconda::biopython bioconda::eggnog-mapper=2.0.1 '

  publishDir "${params.output}/${name}/annotate/rna_annotation/", mode: 'copy', pattern: "*.tsv"
  errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
  maxRetries = 5
  input:
    tuple val(name), path(transcript), path(quant), path(db)
  output:
    tuple val(name), path("*.annotations.tsv"), path(quant)
    path("*.seed_orthologs.tsv")
  shell:
    """
    emapper.py --data_dir ${db} -d bact -o ${name}_transcript  -m diamond -i ${transcript} --cpu ${task.cpus} --go_evidence non-electronic  --target_orthologs all --translate
    tac ${name}_transcript.emapper.annotations | sed "1,3d" | tac |sed "1,3d" > ${name}_transcript.annotations.tsv
    cp ${name}_transcript.emapper.seed_orthologs ${name}_transcript.seed_orthologs.tsv
    """
}