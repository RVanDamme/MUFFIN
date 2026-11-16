<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Recurring Failures, Customized Pipelines, and Fix Handoff

Read this when the user presents several failed runs of the same pipeline,
asks which failures keep coming back, or wants a code fix proposed (issue,
branch, or pull request) rather than only a diagnosis.

## 1. Decide whether a code change is warranted at all

Classify each failure before proposing any edit:

| Class | Typical evidence | Code change? |
|---|---|---|
| **Pipeline defect** | same process + same first error line across runs and inputs; parse/logic error; wrong channel shape; tool called with bad options | yes |
| **Environment / config** | IAM or credential denial, missing container, compute-env limits, wrong profile, quota | no — fix config or infrastructure; say where |
| **Transient remote-data access** | HTTP 5xx / 429, `Connection reset`, `SocketTimeout`, S3 `SlowDown`, FTP/SRA fetch timeouts, succeeded on retry or resume with the same inputs | no — retry/resume; at most suggest `errorStrategy`/`maxRetries` on that fetch step if it recurs |
| **Input-specific** | fails for one sample only, domain error from the tool about that file | usually no — fix the input; a code change only if the pipeline should validate or handle it |

Proposing a code change for a transient or environment failure is a false
positive. Say plainly "no code change needed" when that is the answer.

## 2. Recurrence

Normalize each failure to a fingerprint (process name without sample tag +
first real error line, IDs and paths stripped — see
`fingerprints-and-cascades.md`). Group runs by fingerprint and rank:

1. fingerprints seen in several runs with different inputs (likely defect)
2. fingerprints seen repeatedly with the same inputs (defect or input)
3. one-off fingerprints (likely transient; confirm with a later success)

Only use runs the user pointed at or that the workspace listing returns for
this pipeline; state how many runs the ranking is based on.

## 3. Respect the customization

When the pipeline is a fork that diverged from nf-core (or any upstream):

- Before editing, read the git history for the failing module:
  the affected file's Git history, tags, and CHANGELOG entries. Find out *why*
  the local version differs (performance, a vendor tool wrapper, splitting
  a long task into distributed chunks).
- Do not "fix" by reverting to the upstream module. A fix that removes the
  customization (e.g. collapsing distributed processing back into one long
  task) is a regression even if it makes the error go away.
- Fix the root cause, not the symptom: no `errorStrategy 'ignore'`, no
  catching and discarding the exception, no special-casing one sample ID,
  unless the user asks for a workaround and it is labelled as one.
- If the right fix belongs upstream (a shared nf-core module), say so and
  propose the upstream change separately.

## 4. Hand off for human review

The deliverable is a reviewable change, never a merge:

1. Create a branch using the user's selected naming convention.
2. Make the smallest change that fixes the root cause. Add or update a test
   (nf-test with the failing input pattern) that fails before and passes
   after, when test data allows.
3. Run lint and the relevant tests locally or in the sandbox, and report the
   result.
4. Open a pull request (or issue, if the user asked for a ticket or no fix
   is warranted) containing: failing run IDs/URLs, fingerprint, class,
   evidence excerpt, root cause, why the customization is preserved, the
   change, and how it was tested.
5. Do not merge, approve, or push to the default branch.

If a tool for creating the issue or pull request is not available in this
session, write the ticket body and the patch in the response so the user can
file it.
