process sourmash_download_db {
  if (workflow.profile.contains('gcloud')) {publishDir 'gs://gcloud_storage/databases-nextflow/sourmash', mode: 'copy', pattern: "gtdb-rs202-k31.lca.json.gz" }
  else { storeDir 'nextflow-autodownload-databases/sourmash' }
  //this condition is here only for gcloud usage
  errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
  maxRetries = 5
  label 'ubuntu'
  output:
  path("gtdb-rs202-k31.lca.json.gz")
  script:
  """
  wget --no-check-certificate https://osf.io/9xdg2/download -O gtdb-rs202-k31.lca.json.gz
  """
}

    // The link use from osf.io require sourmash V3 (sourmash V3.3.0 is working))
