<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Failure Fingerprints and Cascades

Use this reference when many tasks failed, logs are noisy, or failures appear to
move between processes after retry/resume.

## Build a fingerprint

For each failed task or attempt, extract:

```text
normalized_process + first_real_error_line + failure_category
```

Where:

- `normalized_process` removes sample tags, retry suffixes like `(1)`, and long
  IDs when comparing failures.
- `first_real_error_line` is the first deterministic exception/tool/platform
  error, not the final Nextflow wrapper line.
- `failure_category` is the category from `failure-patterns.md`.

## Find the root fingerprint

When multiple tasks failed:

1. Group tasks by fingerprint.
2. Prefer the earliest deterministic platform/tool error over later wrapper
   failures.
3. Look for shared artifacts, paths, buckets, references, sample IDs, or outputs.
4. If downstream tasks complain about missing files produced by an earlier
   failed task, the upstream producer is the root cause.
5. If many tasks fail with the same permission/path error, report one shared
   storage/staging root cause.
6. Summarize secondary failures as cascade symptoms.

## Cascade patterns

| Pattern | Root cause | Do not blame |
| --- | --- | --- |
| Upstream process fails, downstream tasks report missing output | Upstream process/tool/resource/staging issue | Each downstream missing-output task |
| Many tasks show the same `AccessDenied` path | Shared IAM/data-link/bucket policy issue | Individual tools |
| Many tasks show the same `NoSuchKey` or `No such file` object | Shared path/sample sheet/object-key issue | Scheduler or compute env |
| Workflow-level validation fails before tasks | Config/params/sample sheet issue | Any task process |
| First task OOMs, later tasks cannot find its output | OOM in first task | Downstream missing files |

## Answer language

- "These failures share the same fingerprint, so I would treat them as one root
  cause rather than independent task errors."
- "The downstream missing-file errors are cascade symptoms; the upstream process
  that should have produced that file failed first."
- "I would fix/rerun the upstream producer, then resume the workflow so completed
  tasks can be reused."
- "If the logs are truncated, the next evidence needed is the first exception in
  `.command.err` for the earliest failed process."
