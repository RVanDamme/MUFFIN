process eggnog_download_db {
        
        if (workflow.profile == 'conda') { storeDir 'nextflow-autodownload-databases/eggnog' }
        else if (workflow.profile == 'gcloud') {publishDir 'gs://nf-muffin20/databases-nextflow/eggnog', mode: 'copy', pattern: "eggnog-db"}
        else { publishDir 'nextflow-autodownload-databases/eggnog', mode: 'copy', pattern: "eggnog-db" }
        label 'eggnog' 
      output:
        path("eggnog-db")
      script:
        """
        mkdir eggnog-db
        download_eggnog_data.py --data_dir eggnog-db -y
        """
    }