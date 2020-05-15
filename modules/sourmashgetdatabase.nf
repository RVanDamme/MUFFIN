process sourmash_download_db {
        if (workflow.profile == 'gcloud') { publishDir 'gs://databases-nextflow/databases/sourmash', mode: 'copy', pattern: "genbank-k31.lca.json.gz" }
        else if (workflow.profile == 'gcloud') {publishDir 'gs://nf-muffin20/databases-nextflow/sourmash', mode: 'copy', pattern: "genbank-k31.lca.json.gz"}
        else { storeDir 'nextflow-autodownload-databases/sourmash' }  
        //this condition is here only for gcloud usage 
        label 'ubuntu' 
      output:
        file("genbank-k31.lca.json.gz")
      script:
        """
        #wget https://ndownloader.figshare.com/files/18809423?private_link=ed98a281ef089c033352 -O gtdb.lca.json
        wget https://osf.io/4f8n3/download -O genbank-k31.lca.json.gz 
        """
    }

    The link use from osf.io require sourmash V3 (sourmash V3.3.0 is working))