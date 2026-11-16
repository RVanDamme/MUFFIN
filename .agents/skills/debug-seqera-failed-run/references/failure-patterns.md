<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Failed Pipeline Diagnosis Patterns

## Diagnosis rules

- **Separate root cause from symptoms.** If the workflow log shows launch-time
  validation, permissions, or staging errors, do not blame a downstream task that
  merely failed after the missing artifact/permission problem.
- **Prefer explicit log evidence over exit-code guesses.** Exit status 137 often
  means OOM, but Platform may report OOM as `OutOfMemoryError` with another exit
  code (for example 247). Empty output plus `Host EC2 ... terminated` points to
  host loss, not a tool error.
- **Classify the failure type** in the answer: infrastructure/resource,
  storage/auth, transient compute, workflow script/configuration, or
  bioinformatics/input-data/tool error.
- **Ignore noisy non-fatal warnings** when a later traceback or platform error
  gives a concrete root cause (for example pandas warnings before `Errno 28`).

## Category-first triage

Triage in this order, stopping at the first category with strong log evidence:

1. **Launch / workflow-level failure:** malformed params, sample sheet,
   workflow config, credentials, plugin download, or Nextflow script compilation
   before a task really starts.
2. **Storage / staging / publish failure:** missing object, bad path, data link
   mismatch, IAM denial, or local/cloud path mismatch.
3. **Transient compute / scheduler failure:** EC2 host termination, spot
   interruption, external scheduler kill, queue/walltime eviction.
4. **Resource exhaustion:** OOM, disk/scratch full, quota, timeout.
5. **Container / registry failure:** image missing or pull denied.
6. **Bioinformatics/input/tool failure:** tool rejected FASTQ/BAM/VCF/reference,
   sample metadata, or command-line options.
7. **Unknown:** give the best evidence-backed hypothesis and say what extra log
   would disambiguate it.

## Common Platform failure patterns

| Pattern | Cause | Fix |
|---------|-------|-----|
| `OutOfMemoryError`, `Out of memory`, `Killed`, `Exit status 137`, `oom-kill` | OOM killed | Increase process `memory`; for Java tools also increase `-Xmx`/heap consistently; resume the run |
| `No space left on device`, `Errno 28`, disk quota, `quota exceeded`, `not enough space` | Disk or scratch full | Increase task disk/scratch, point `TMPDIR` at a larger volume, or clean work dirs |
| `Host EC2 ... terminated`, missing exit status, empty output | Transient host/preemption | Resume/retry; check AWS Batch/spot interruption and CE stability if repeated |
| `terminated for an unknown reason -- Likely it has been terminated by the external system`, `exit status 143`, `SIGTERM`, `preempt`, `evict` | External scheduler termination/preemption/time limit | Retry/resume; inspect scheduler events, spot interruptions, and walltime limits |
| `AccessDenied`, `Access Denied`, `not authorized`, `403`, `Permission denied`, `Unable to locate credentials`, `ExpiredToken` | IAM / storage / registry permissions | Fix compute environment role, bucket policy, data link credentials, temporary credentials, or registry login |
| `No such file or directory`, `FileNotFoundException`, `NoSuchKey`, `Cannot find any reads matching`, `Unable to access path`, `MissingFileException` | Missing input/path/staging issue | Verify sample sheet paths, data links, object keys, and staging/publish configuration |
| `Malformed row`, `Invalid number of columns`, `samplesheet`, `schema validation`, `invalid accession`, incompatible reference | Bioinformatics input validation | Fix sample sheet/reference/input format according to pipeline docs |
| `No such variable`, `Script compilation error`, `Unknown process`, `Channel not found`, `Unknown parameter`, `Cannot invoke method`, `No such property` | Nextflow script/configuration error | Fix pipeline/module code or parameter names/types before retrying |
| `pull access denied`, `manifest unknown`, `image not found`, `repository does not exist`, `unauthorized: authentication required` | Container pull/auth problem | Verify container name/tag and registry credentials |
| `Exit status 139`, `segmentation fault`, `core dumped` | Tool/container crash | Check tool inputs and container version; consider updating/pinning image |
| `Exit status 1` + stderr | Tool/script error | Read task log for the actual exception/error before classifying |

## Key log patterns

- **Script/config errors:** `Script compilation error`, `No such variable`,
  `Missing process`, `Channel not found`, `Unknown parameter`, `No such
  property`, `Cannot invoke method`
- **Staging/storage errors:** `FileNotFoundException`, `NoSuchKey`,
  `Cannot find any reads matching`, `MissingFileException`, `AccessDenied`,
  `s3:`, `gs://`, `az://`
- **Resource/infra errors:** `OutOfMemoryError`, `killed`, `No space left`,
  `quota`, `timeout`, `Host EC2`, `terminated by the external system`,
  `SIGTERM`, `preempt`, `evict`
- **Bio/input errors:** `Malformed row`, `invalid accession`, malformed
  FASTQ/FASTA/VCF/BAM/reference messages from the tool, invalid contig/sample
  names, reference/index mismatch, schema validation failure
- **Container errors:** `pull access denied`, `manifest unknown`, `not found`,
  `unauthorized: authentication required`

## Resume / retry guidance

Recommend resume only after considering the category:

- **Transient compute/preemption:** resume/retry is usually the first action,
  because work may be reusable and the root cause may be host loss or scheduler
  eviction.
- **OOM or disk/scratch:** resume after increasing `memory`, heap, `disk`, or
  scratch/TMPDIR. A plain resume with unchanged resources is likely to fail in
  the same place.
- **Missing input/path/staging:** fix the sample sheet, object key, data link, or
  staging/publish path first; then resume.
- **Permissions/auth:** fix IAM, bucket policy, data link credentials, temporary
  credentials, or registry auth first; then resume.
- **Script/config/params:** fix the pipeline code/config/parameter names or
  values first. Resume may reuse completed work, but it is not the fix.
- **Bioinformatics/tool/input:** fix the rejected FASTQ/BAM/VCF/reference,
  sample metadata, or tool option before resuming.

If a resumed run fails with a different category or first-error fingerprint,
explain that the first issue may have been unblocked and the new failure should
be debugged as the current root cause.

## Handling ambiguous failures

When logs do not match a clear platform/resource pattern:

1. Extract the **failing process** from `Error executing process > '...'`.
2. Look at the **first real error line** in stderr, not only the final Nextflow
   wrapper line.
3. Mention the **repository/pipeline family** if available, because some errors
   are tool- or pipeline-specific.
4. Normalize repeated suffixes like `(1)`, sample IDs, and task tags before
   comparing multiple failures.
5. If several tasks failed, identify whether they share the same fingerprint
   (same process + first error line) or are cascade symptoms.
6. If classification is uncertain, say so explicitly and ask for the task stderr
   / `.command.err` / workflow log section needed to decide.

Useful answer language for ambiguous tool failures:

- "This appears to be a tool-level failure inside `<process>`, not an infra
  failure, because the container started and produced a domain-specific error."
- "I don't see evidence of OOM, disk, host termination, or permissions; the next
  thing to inspect is the first exception in `.command.err`."
- "Multiple downstream failures share the same missing input, so the root cause
  is the upstream staging/path problem rather than each failed task."
