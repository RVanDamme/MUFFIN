process eggnog { 
        label 'eggnog' 
      input:
        val(sample), val(bin_id), file(bin)
        file(db)
      output:
        val(sample), val(bin_id), file("*.annotations.tsv")
      shell:
        """
        emapper.py --data-dir ${db} -d bact -o ${bin_id} --output_dir ${bin_id}_out -m diamond -i ${bin} --cpu ${task.cpus} --go_evidence non-electronic  --target_orthologs all --translate
        tac ${bin_id}.annotations | sed "1,3d" | tac |sed "1,3d" > ${bin_id}.annotations.tsv
        cp ${bin_id}.seed_orthologs ${bin_id}.seed_orthologs.tsv
        """
    }