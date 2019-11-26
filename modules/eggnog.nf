process eggnog { 
        label 'eggnog' 
      input:
        set val(sample), val(bin_id), file(bin)
        file(db)
      output:
        set val(sample), val(bin_id), file("*.annotations.tsv")
      shell:
        """
        emapper.py --data_dir ${db} -d bact -o ${bin_id}  -m diamond -i ${bin} --cpu ${task.cpus} --go_evidence non-electronic  --target_orthologs all --translate
        tac ${bin_id}.emapper.annotations | sed "1,3d" | tac |sed "1,3d" > ${bin_id}.annotations.tsv
        cp ${bin_id}.emapper.seed_orthologs ${bin_id}.seed_orthologs.tsv
        """
    }