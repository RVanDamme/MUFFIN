process sourmash_download_db {
        if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/sourmash', mode: 'copy', pattern: "genbank-k31.lca.json" }
        else { storeDir 'nextflow-autodownload-databases/sourmash' }  
        //this condition is here only for gcloud usage 
        label 'ubuntu' 
      output:
        file("gtdb_r89_k31_scaled10k.sbt.json"), file(".sbt.gtdb_r89_k31_scaled10k")
      script:
        """wget https://ndownloader.figshare.com/files/18666761?private_link=f7f668a136b56e04a48d -O gtdb_r89_k31_scaled10k.tar.gz
        tar -xzf gtdb_r89_k31_scaled10k.tar.gz
        """
    }