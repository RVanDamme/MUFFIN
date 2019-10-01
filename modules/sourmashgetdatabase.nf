process sourmash_download_db {
        //if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/sourmash', mode: 'copy', pattern: "genbank-k31.lca.json" }
        //else 
        { storeDir 'nextflow-autodownload-databases/sourmash' }  
        label 'ubuntu' 
      output:
        file("genbank-k31.lca.json")
      script:
        """
        wget https://osf.io/4f8n3/download -O genbank-k31.lca.json.gz
        gunzip genbank-k31.lca.json.gz
        """
    }