<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Fix Confidence and Evidence

Use this reference when the logs are incomplete, the failure lands in an
ambiguous/tool-level category, or the user asks what fix is likely to work.

## Confidence levels

| Confidence | Evidence | How to answer |
| --- | --- | --- |
| High | Explicit decisive line: `OutOfMemoryError`, `NoSuchKey`, `AccessDenied`, `Errno 28`, malformed sample sheet, script compilation error | State the root cause directly and give the concrete fix |
| Medium | Repeated fingerprint, process-level tool error, or partial traceback but missing decisive context | Give the likely diagnosis, mention uncertainty, and name the next log to confirm |
| Low | Only Nextflow wrapper text, generic `exit status 1`, truncated logs, or no stderr | Do not invent a fix; ask for `.command.err`, task stderr/stdout, and workflow log |

## Evidence-first response

For every diagnosis, separate:

1. **Evidence:** short decisive log line(s).
2. **Inference:** what category/root cause that evidence supports.
3. **Confidence:** high/medium/low.
4. **Fix:** concrete action, with resume guidance if appropriate.
5. **Missing evidence:** if confidence is not high.

## Avoid overclaiming

- Generic `exit status 1` is not enough to call infra, OOM, permissions, or a
  pipeline bug.
- A container that starts and emits domain-specific stderr is usually a
  tool/input-level failure unless there is explicit platform/resource evidence.
- Warnings before a later exception are usually noise; quote the exception.
- If logs only show Nextflow wrapper lines, ask for `.command.err` or task stderr.

## Answer language

- "Confidence: high — the log explicitly says ..."
- "Confidence: medium — this looks tool/input-level because the container
  started and failed inside `<process>`, but I need `.command.err` to identify
  the exact bad input/option."
- "Confidence: low — I only have the Nextflow wrapper error. I would not change
  memory/IAM/config yet without task stderr."
- "The next evidence to fetch is the first exception line from `.command.err` for
  `<process>`."
