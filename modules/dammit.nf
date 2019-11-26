process dammit { 
        label 'dammit' 
        publishDir "${params.output}/${name}/dammit/", mode: 'copy', pattern: "*.dammit.*"
      input:
        set val(name), file(transcript)
        file(db)
      output:
        set val(name), file("*.dammit.fasta"), file("*.dammit.gff3")
      script:
        if (params.dammit_user_db) {"""
        dammit annotate ${transcript} --user-databases ${params.dammit_user_db} --busco-group ${params.busco_db}  -o ${name}_out --database-dir ${db}
        cp ${name}_out/${transcript}.dammit.fasta ${name}.dammit.fasta
        cp ${name}_out/${transcript}.dammit.gff3 ${name}.dammit.gff3
        """}
        else {"""
        dammit annotate ${transcript}  --busco-group ${params.busco_db}  -o ${name}_out --n_threads ${task.cpus} --database-dir ${db}
        cp ${name}_out/${transcript}.dammit.fasta ${name}.dammit.fasta
        cp ${name}_out/${transcript}.dammit.gff3 ${name}.dammit.gff3
        """}
    }
            //--n newname for the transcrip?
            // --force?
            // --quick as a param?

            // SHOULD BE
        // dammit annotate ${transcript} --user-databases ${params.dammit_user_db} --busco-group ${params.busco_db} --n_threads ${task.cpus} -o ${name}_out --database-dir ${db} 
        // cp ${name}_out/${transcript}.dammit.fasta ${name}.dammit.fasta
        // cp ${name}_out/${transcript}.dammit.gff3 ${name}.dammit.gff3
