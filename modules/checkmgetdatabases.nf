process checkm_download_db {
        // if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/checkm', mode: 'copy', pattern: "genbank-k31.lca.json" }
        if (workflow.profile == 'gcloud') {publishDir 'gs://nf-muffin20/databases-nextflow/checkm', mode: 'copy', pattern: "checkm_data_2015_01_16.tar.gz"}
        else { storeDir 'nextflow-autodownload-databases/checkm' }
        label 'ubuntu' 
      output:
        file("checkm_data_2015_01_16.tar.gz")
      script:
        """
        wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
        """
    }