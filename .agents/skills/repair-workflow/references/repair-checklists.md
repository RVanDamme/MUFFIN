<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Repair Workflow Checklists

Use these compact checklists when repairing an existing Nextflow workflow.

## Validation family chooser

Pick the narrowest validation loop that actually matches the task.

### A. Existing verifier / nf-test / task-local check

Use this when:
- the repo already includes an intended verifier
- the failing surface already has nf-test coverage
- the user is fixing an existing module/subworkflow/workflow surface

Preferred loop:
1. run the intended verifier early
2. inspect the failure on the intended surface
3. make the smallest in-place edit
4. rerun the same verifier

Do not replace this with generic repeated linting.

### B. Explicit nf-core pipeline lint requirement

Use this when:
- instructions explicitly mention `nf-core pipelines lint`
- the repo clearly expects nf-core pipeline metadata/config hygiene
- files like `.nf-core.yml` and `manifest.name` are part of the task surface

Preferred loop:
1. bootstrap only the minimal missing metadata/config
2. run the workflow once if needed
3. run the exact required lint command
4. bucket failures and apply only the mapped repair

### C. Channel-shape / join / fanout / optional-side-channel bug

Use this when:
- outputs are duplicated, missing, or leaking right-only rows
- the bug involves `join`, `combine`, regrouping, fanout/fanin, or optional annotations

Preferred loop:
1. write down example input tuples
2. write down expected output tuples
3. identify the join key
4. identify the authoritative stream
5. then edit operators

### D. Strict syntax / parser migration

Use this when:
- errors mention strict syntax, parser v2, missing type annotation, or unexpected token

Preferred loop:
1. load `migrate-nextflow-code`
2. repair parser issues directly
3. rerun the relevant Nextflow lint or test command

## Deterministic nf-core lint-repair map

When lint failures are metadata/config drift rather than workflow logic, repair these first:

- stale `.nf-core.yml`
- missing `custom_config_version = 'master'`
- missing `custom_config_base = 'https://raw.githubusercontent.com/nf-core/configs/master'`
- missing custom profile include
- wrong `dag.file` extension

Bucket failures into:
- config drift
- metadata drift
- missing validation profile or config
- structural module/subworkflow issue
- stylistic warning only

Rule:
- after one mapped repair, rerun lint once
- if the same bucket remains twice, stop retrying and re-inspect the surface

## Stuck-task checklist

If two repair cycles did not materially change the failure:

1. reread the verifier/test output
2. inspect the exact file or surface the verifier is targeting
3. inspect produced outputs or intermediate artifacts
4. compare expected vs actual tuple/file shapes
5. write a one-paragraph failure hypothesis
6. only then make another edit

Do not keep increasing retry count without changing strategy.

## Anti-pattern rejection checklist

Reject or immediately rewrite edits that introduce:
- helper invocations that assume unverified project-relative locations
- manual copy-back into an assumed output folder
- shadow rewrites of the intended nf-core module or subworkflow
- trivial script-wrapper passthroughs
- broad version bumps or broad snapshot rewrites without specific evidence

## Minimal done criteria

Before declaring the repair done:
- the intended validation loop ran
- the fix landed on the intended existing surface
- no fixture rewrite hid the bug
- no wrapper/copy-back/shadow-surface hack was introduced
- the final failure mode changed from red to green, not merely from one retry to another
