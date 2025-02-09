process checkm_download_db {
        
  //if (workflow.profile.contains('gcloud')) {publishDir 'gs://gcloud_storage/databases-nextflow/checkm2', mode: 'copy', pattern: "checkm_data_2015_01_16.tar.gz"}
  //else { storeDir 'nextflow-autodownload-databases/checkm2' }
  label 'checkm2' 

  conda 'bioconda::checkm2=1.0.1'

  errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
  maxRetries = 5
  output:
    val(true)
  script:
  """
  checkm2 database --download ${params.db_path}
  """
}