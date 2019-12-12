process sourmash_download_db {
        if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/sourmash', mode: 'copy', pattern: "gtdb.lca.json" }
        else { storeDir 'nextflow-autodownload-databases/sourmash' }  
        //this condition is here only for gcloud usage 
        label 'ubuntu' 
      output:
        file("gtdb.lca.json")
      script:
        """wget https://ndownloader.figshare.com/files/18809423?private_link=ed98a281ef089c033352 -O gtdb.lca.json
        """
    }