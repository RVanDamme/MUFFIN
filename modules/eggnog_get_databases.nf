process eggnog_download_db {
        // if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/dammit', mode: 'copy', pattern: "" }
        // else { storeDir......}
        storeDir 'nextflow-autodownload-databases/eggnog' 
        label 'eggnog' 
      output:
        file("eggnogdb")
      script:
        """
        mkdir eggnog-db
        download_eggnog_data.py --data-dir eggnog-db -y
        """
    }