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
        dammit databases --install --n_threads 1 --busco-group ${params.busco_db} --database-dir dammit-db
        """
    }