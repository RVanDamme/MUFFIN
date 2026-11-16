<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../../../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Repair Failing nf-tests

Fix nf-test failures using evidence, not guesswork. Every decision must trace back
to concrete test output.

## When to Use

Load this skill when the user:
- reports one or more failing nf-tests
- shares nf-test output showing failures, snapshot mismatches, or assertion errors
- asks to fix, debug, or repair their nf-test suite
- needs help understanding why a test broke after a pipeline change

Do not use this skill when:
- the user wants to write new tests from scratch → use `nf-test`
- the user wants to fix pipeline logic unrelated to tests → use `repair-workflow`

## Cardinal Rule: Evidence Over Speculation

**Do not guess.** Every hypothesis must cite specific evidence from the test output.

Anti-patterns that this skill exists to prevent:
- Changing Nextflow version, profiles, or executor config without evidence they caused the failure
- Rewriting pipeline logic when only test expectations changed
- Modifying `maxForks`, GPU profiles, resource directives, or container versions speculatively
- Updating snapshots blindly without understanding what changed and why
- Making broad changes across multiple files when the failure points to a single assertion

## Repair Protocol

### Step 1: Examine the failure evidence

Read the full nf-test output. Extract these facts before doing anything else:

1. **Which test(s) failed?** — file path, test name, assertion line
2. **What type of failure?**
   - Snapshot mismatch → what changed in the snapshot diff?
   - Assertion error → which assertion, what was expected vs actual?
   - Process failure → exit code, stderr content
   - Infrastructure error → missing files, container issues, timeout
3. **What specifically differs?** — for snapshot mismatches, identify the exact fields/values that changed (e.g., process count went from 5 to 6, MD5 hash differs, output file name changed)

If the user has not provided sufficient test output:
- Ask for the full `nf-test test --verbose` output
- Ask for the `.nf.test.snap` diff if it's a snapshot mismatch
- Ask for `.command.log` and `.command.err` if a process failed
- **Stop and ask.** Do not proceed without evidence.

```
Evidence checklist:
☐ Test file path and name identified
☐ Failure type classified
☐ Specific diff/mismatch values extracted
☐ Root cause hypothesis formed from evidence (not assumption)
```

### Step 2: Classify the failure and determine the correct response

Based on the evidence from Step 1, classify into exactly one category:

#### A. Expected change — update the test

The pipeline was intentionally changed and the test needs to reflect reality.

Evidence patterns:
- A new process was added to the workflow → process count increases
- A tool version was bumped → output file checksums change
- An output channel was renamed or restructured → snapshot keys change
- A new output was added → snapshot gains new entries

Action: update the snapshot or assertion only after checking the intentional contract change against independent expected behavior. Select the affected test using the project's supported selector before passing `--update-snapshot`; the command below shows the flag, not a suite-wide update instruction.

```bash
nf-test test --update-snapshot
```

Then verify the updated snapshot looks correct — do not blindly accept.

#### B. Regression — fix the code

The pipeline change introduced a bug that the test correctly caught.

Evidence patterns:
- A process that should run is now missing from the trace
- Output files are empty or malformed
- Exit code is non-zero with a meaningful error in `.command.err`
- Channel wiring is broken (wrong input, missing input)

Action: fix the pipeline code, not the test. The test is correct.

#### C. Non-deterministic output — stabilize the test

The test is flaky because it asserts on unstable values.

Evidence patterns:
- Timestamps, dates, or UUIDs in output differ between runs
- File ordering is non-deterministic
- Floating-point precision varies

Action: filter out unstable values or use targeted assertions instead of full snapshots.

```groovy
// Instead of snapshotting everything:
assert snapshot(process.out).match()

// Snapshot only stable outputs:
assert snapshot(process.out.versions).match("versions")
assert snapshot(
    process.out.report.collect { path(it[1]).readLines().findAll { !it.contains("timestamp") } }
).match("report")
```

#### D. Environment/infrastructure — fix the setup

The test itself is fine but the environment is wrong.

Evidence patterns:
- "No such file" errors for test data
- Container pull failures
- Profile not found
- Out of memory / timeout

Action: fix configuration, not pipeline code or test assertions.

### Step 3: Plan precise changes

Before editing any file, write out:

1. The exact file(s) to change
2. The specific edit in each file
3. Why this edit addresses the evidence from Step 1
4. What the expected test result will be after the edit

If you cannot articulate all four points, return to Step 1 and gather more evidence.

**Scope rule:** each planned change must address exactly one failure cause. Do not
bundle unrelated fixes.

### Step 4: Apply one change at a time

Make the smallest edit that addresses the identified cause. Then immediately
proceed to Step 5.

Do not:
- Make multiple unrelated changes before testing
- Refactor code that is not related to the failure
- Change Nextflow version, profiles, or executor config unless Step 1 evidence specifically implicates them

### Step 5: Run the tests

```bash
# Run only the failing test(s) first
nf-test test --verbose

# If that passes, run the full suite to check for regressions
nf-test test --verbose
```

Capture the full output. Compare against the Step 1 baseline.

### Step 6: Evaluate the result

Three possible outcomes:

| Outcome | Evidence | Action |
|---------|----------|--------|
| **Fixed** | Previously failing test now passes, no new failures | Proceed to Step 7 |
| **Improved** | Failure changed in kind (different error, fewer failures) | Return to Step 1 with new output |
| **Regressed** | New failures appeared, or original failure unchanged | Revert the change, return to Step 1 |

"Improved" means the failure mode materially changed — not just a different line
number or a slightly different error message for the same root cause.

### Step 7: Decide next action

- If all target tests pass → run full suite, verify no regressions, declare done
- If more failures remain → return to Step 1 with the new test output
- If the next step is ambiguous → stop and ask the user for clarification

## Snapshot-Specific Guidance

Snapshot mismatches are the most common nf-test failure. Handle them carefully:

### Reading a snapshot diff

The `.nf.test.snap` file is JSON. When a mismatch occurs, compare old vs new:

- **Process count changed** → a process was added/removed from the workflow. Check if intentional.
- **MD5 hash changed** → file content differs. Check if a tool version changed or if output is non-deterministic.
- **Channel structure changed** → output channel names or tuple shapes differ. Check workflow `emit:` blocks.
- **New keys appeared** → new outputs were added. Check if intentional.
- **Keys disappeared** → outputs were removed. This is likely a regression.

### When to update vs when to fix

Update the snapshot when:
- You have confirmed the pipeline change was intentional
- The new snapshot values are correct (not just different)
- You can explain why each changed value is expected

Fix the code when:
- The snapshot change reveals missing or malformed output
- Process count decreased unexpectedly
- MD5 checksums differ and you cannot trace the cause to an intentional change

### Never update blindly

Running `--update-snapshot` without understanding the diff is the single most
common mistake. Always read the diff first.

## Anti-Patterns to Reject

Stop immediately if you find yourself:

| Anti-pattern | Why it's wrong |
|-------------|---------------|
| Blaming Nextflow caching for a failed nf-test | nf-test always runs fresh in its own nf-test's managed work directory work directory without `-resume`. Nextflow's resume cache has zero effect on nf-test results. A failed nf-test is never a "cached result" — treat every failure as real. |
| Changing Nextflow version | Version changes cause cascading snapshot diffs unrelated to the bug |
| Modifying executor/profile config | Profiles rarely cause assertion failures |
| Adjusting `maxForks`, `cpus`, `memory` | Resource settings don't change test assertions |
| Rewriting pipeline logic for a snapshot mismatch | The snapshot caught a real change — understand it first |
| Adding/removing test profiles | Test profiles are infrastructure, not the cause of assertion failures |
| Updating tool versions speculatively | Tool bumps change checksums everywhere |
| Making changes to files not referenced in the error | Unrelated edits add noise and risk regressions |

## Escalation

If after two repair cycles (Step 1 → Step 6) the failure persists:

1. Re-read the full test output from scratch
2. Inspect the actual task outputs reported by the affected test
3. Compare the `.command.sh` that ran against what you expected
4. Read the `.nf.test.snap` file directly (not just the diff)
5. Write a one-paragraph hypothesis explaining the failure
6. Only then make another edit — or ask the user for help

Do not keep retrying the same approach. If the failure doesn't change after two
attempts, your hypothesis is wrong.

## Quick Reference

Use nf-test test with verbose or debug output and select the affected tests using
the user's existing test-selection mechanism. Review snapshot differences before
updating. Inspect .command.sh, .command.log, and .command.err for the actual task
reported by the test; do not assume their location.

~~~bash
nf-test test --verbose
nf-test test --debug
~~~

## Related Skills

- `nf-test` for writing new tests, assertion reference, and plugin APIs
- `repair-workflow` for fixing pipeline logic unrelated to test assertions
- `nf-pipeline-design` for structural design rules
