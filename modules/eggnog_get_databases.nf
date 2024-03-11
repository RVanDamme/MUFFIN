process eggnog_download_db {
        
  if (workflow.profile.contains('conda')) { storeDir 'nextflow-autodownload-databases/eggnog' }
  else if (workflow.profile.contains('gcloud')) {publishDir 'gs://gcloud_storage/databases-nextflow/eggnog', mode: 'copy', pattern: "eggnog-db"}
  else { publishDir 'nextflow-autodownload-databases/eggnog', mode: 'copy', pattern: "eggnog-db" }
  label 'eggnog' 

  conda 'bioconda::diamond anaconda::biopython bioconda::eggnog-mapper=2.0.1 '
  
  errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
  maxRetries = 5
  output:
    path("eggnog-db")
  script:
  """
  mkdir eggnog-db
  download_eggnog_data.py --data_dir eggnog-db -y
  """
}