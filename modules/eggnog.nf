process eggnog_bin { 
        label 'eggnog' 
        publishDir "${params.output}/${name}/bin_annotated/${bin_id}/", mode: 'copy', pattern: "*.tsv"
      input:
        set val(name), val(bin_id), file(bin), file(db)
      output:
        set val(name), val(bin_id), file("*.annotations.tsv")
        file("*.seed_orthologs.tsv")
      shell:
        """
        emapper.py --data_dir ${db} -d bact -o ${bin_id}  -m diamond -i ${bin} --cpu ${task.cpus} --go_evidence non-electronic  --target_orthologs all --translate
        tac ${bin_id}.emapper.annotations | sed "1,3d" | tac |sed "1,3d" > ${bin_id}.annotations.tsv
        cp ${bin_id}.emapper.seed_orthologs ${bin_id}.seed_orthologs.tsv
        """
    }

process eggnog_rna { 
        label 'eggnog' 
        publishDir "${params.output}/${name}/rna_annotated/", mode: 'copy', pattern: "*.tsv"
      input:
        set val(name), val(transcript), file(quant), file(db)
      output:
        set val(name), file("*.annotations.tsv"), file(quant)
        file("*.seed_orthologs.tsv")
      shell:
        """
        emapper.py --data_dir ${db} -d bact -o ${name}_transcript  -m diamond -i ${transcript} --cpu ${task.cpus} --go_evidence non-electronic  --target_orthologs all --translate
        tac ${name}_transcript.emapper.annotations | sed "1,3d" | tac |sed "1,3d" > ${name}_transcript.annotations.tsv
        cp ${name}_transcript.emapper.seed_orthologs ${name}_transcript.seed_orthologs.tsv
        """
    }