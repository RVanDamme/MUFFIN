process checkm_download_db {
        
  if (workflow.profile.contains('gcloud')) {publishDir 'gs://gcloud_storage/databases-nextflow/checkm', mode: 'copy', pattern: "checkm_data_2015_01_16.tar.gz"}
  else { storeDir 'nextflow-autodownload-databases/checkm' }
  label 'ubuntu' 
  errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
  maxRetries = 5
  output:
    path("checkm_data_2015_01_16.tar.gz")
  script:
  """
  wget --no-check-certificate https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
  """
}