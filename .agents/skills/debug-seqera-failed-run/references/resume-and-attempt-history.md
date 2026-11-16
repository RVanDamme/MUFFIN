<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Resume and Attempt History

Use this reference when the user says they retried/resumed, or multiple attempts
appear to share the same Nextflow session, run name, launch lineage, or workflow
family.

## Why this matters

A resumed failed attempt does **not** prove resume fixed anything; it often
means the same issue was retried and failed again.

If a later resumed attempt succeeded, compare concrete deltas between the failed
and successful attempts. Changed CPU, memory, command line, workflow params, or
work directory can be useful clues, but a changed field is not automatically the
causal fix unless it matches the failure evidence.

## Attempt-chain workflow

1. Identify the current failed attempt from the workflow details returned by
   Platform.
2. Find earlier and later attempts in the same chain when available:
   - prefer the same Nextflow session when returned by Platform
   - otherwise use same workspace + run name / launch lineage / repository + user
     + nearby timestamps
   - if a later resumed attempt succeeded, compare the failed attempt with the
     successful one before giving generic advice
3. Compare attempts:
   - failure category
   - failing process
   - first real error line
   - exit status
   - resource requests (`memory`, heap, `disk`, `cpus`)
   - workflow params / command line, preferably at individual parameter-key
     level
   - work directory / launch directory / input path params
   - workflow revision / commit
   - profile, compute environment, Wave/Fusion flags, or container image if
     available
4. Explain whether the attempted resume changed the failure or appears to have
   unblocked it.

## Failure fingerprint

Build a compact fingerprint for each failed attempt:

```text
normalized_process + first_real_error_line + failure_category
```

Normalize away sample IDs, task tags like `(1)`, and long object paths when
comparing retries.

## Interpretation

| Observation | Interpretation | Guidance |
| --- | --- | --- |
| Same category and same fingerprint after resume | The previous attempted fix likely did not address root cause | Do not recommend another plain resume; change the underlying resource/path/config/input first |
| Same category but different process | Possible cascade or next task reached after partial progress | Identify the earliest/root process and shared error |
| Different category after resume | The first issue may have been fixed or bypassed | Debug the new category as the current root cause |
| OOM/disk repeats after resume | Resources/storage still insufficient | Increase memory/heap/disk/scratch/TMPDIR, then resume |
| Transient host loss disappears after resume | Resume likely worked for the transient issue | Continue debugging only if a new deterministic failure appears |
| Missing path/permission/config repeats | Resume cannot fix deterministic setup problems | Fix object path, data link, IAM, config, or params before resuming |
| Later resumed attempt succeeded and params/command line changed | The likely fix may have been corrected inputs, paths, or launch parameters | Name the changed fields and tie them to the failed log evidence |
| Later resumed attempt succeeded and resource requests changed | Resource tuning may have unblocked the run, especially for OOM/disk/transient compute | Mention the before/after resource values and whether they match the failure evidence |
| Later resumed attempt succeeded and revision/commit changed | A pipeline/config/code fix may have been applied | Point to the changed revision/commit, but avoid claiming details absent from logs |
| Later resumed attempt succeeded but only incidental metadata changed | Success is useful context but does not identify the fix | Report the successful resume and ask for config/params/log details if needed |

## Failed-to-successful resume comparisons

When a successful resumed attempt exists, add an **Attempt history** section:

1. State that a later resumed attempt in the same chain succeeded.
2. List concrete deltas: memory/CPU, params, command line, work directory,
   profile/compute settings, Wave/Fusion, container, revision/commit.
3. Rank likely fixes by evidence quality:
   - exact failed-log symptom plus matching successful delta = strongest
   - successful delta without matching log evidence = suggestive, not proven
   - generic retry advice = weakest
4. Use the failure category to focus the comparison:
   - **OOM/disk/resource:** structured memory/CPU/disk fields, explicit
     `max_memory`/`max_cpus`, heap, scratch/TMPDIR, disk, profile
   - **Missing input/path:** `input`, `samplesheet`/`sample_sheet`, `bam`,
     `bed`, `vcf`, `fastq`, bucket/path params, command line, work/launch
     directory, reference params when the missing file is a reference/index
   - **Storage/auth:** `outdir`, work/publish paths, bucket params, data links,
     IAM/credentials-related settings, Fusion/Wave/profile
   - **Script/config:** changed config/options params, config files,
     revision/commit, command line, pipeline-specific params named in the error
   - **Transient compute:** retry/resume success may be enough, but changed
     resources/compute settings are still useful context
5. Deprioritize noisy param deltas unless they directly match the error:
   run IDs, dates, trace suffixes, sample IDs, department, sequencing instrument
   metadata such as `machine`, and unrelated output naming changes.

## Answer language

- "This run was already resumed, and it failed with the same fingerprint. Resume
  alone is unlikely to fix it."
- "The resumed attempt failed differently, so the original issue may have been
  unblocked. The current root cause is now ..."
- "Increase resources/change the input/config first, then resume so completed
  tasks can be reused."
- "I cannot prove the final fix from failed-only history; this conclusion is
  based on comparing failed attempts."
- "A later resumed attempt in the same session succeeded. The successful attempt
  changed `<fields>`, so the likely fix was `<hypothesis>`; confidence is higher
  because this matches the failed log line `<evidence>`."
- "The successful resumed attempt changed resources, but the failed log points to
  a missing input path. I would treat the resource change as context, not the
  primary fix, unless more logs show resource pressure."
- "The key-level param diff is more relevant than the whole params blob: `<param>`
  changed from `<old category>` to `<new category>`, matching the failed log's
  missing input or reference evidence."
- "I see changed metadata params like run ID/trace suffix/machine, but those do
  not explain this error, so I would not treat them as the fix."
