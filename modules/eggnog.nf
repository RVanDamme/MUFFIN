process eggnog { 
        label 'eggnog' 
        publishDir "${params.output}/${name}/checkm_bins/${bin_id}/", mode: 'copy', pattern: "*.tsv"
      input:
        set val(name), val(bin_id), file(bin)
        file(db)
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