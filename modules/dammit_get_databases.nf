process dammit_download_db {
        // if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/dammit', mode: 'copy', pattern: "" }
        // else { storeDir......}
        storeDir 'nextflow-autodownload-databases/dammit' 
        label 'dammit' 
      output:
        file("dammit-db")
      script:
        """
        mkdir dammit-db
        dammit databases --install --databases-dir dammit-db --n_threads ${task.cpus} --busco-group ${params.busco_db}
        """
    }