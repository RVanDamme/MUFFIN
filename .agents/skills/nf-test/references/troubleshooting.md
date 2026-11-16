<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nf-test Troubleshooting Guide

## Important: nf-test is NOT affected by Nextflow caching

nf-test runs each test in its own isolated nf-test's managed work directory work directory and never uses `-resume`. Nextflow's resume cache (the managed Nextflow resume cache) has zero effect on nf-test results. A failing nf-test is never a "cached result" — treat every failure as real and investigate the actual test output.

## Common Errors and Solutions

### "Test data file not found"

**Error:**
```
No such file: <user-provided input reference>
```

**Solutions:**

1. Ensure test data config is loaded:
```groovy
// User-provided Nextflow test configuration
includeConfig 'https://raw.githubusercontent.com/nf-core/modules/master/tests/config/test_data.config'
```

2. Use `checkIfExists: true` to catch early:
```groovy
file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
```

3. Check network access (test data is fetched remotely)

---

### "Snapshot mismatch"

**Error:**
```
Snapshot `main` does not match
```

**Diagnosis:**
```bash
# View detailed diff
nf-test test --verbose
```

**Solutions:**

1. If change is expected, update snapshot:
```bash
nf-test test --update-snapshot
```

2. If unexpected, check for non-deterministic output:
   - Timestamps in output files
   - Random seeds
   - Process ordering
   - File paths

3. Exclude unstable values:
```groovy
assert snapshot(
    path(process.out.log[0][1]).readLines()
        .findAll { !it.contains("timestamp") }
        .findAll { !it.contains(workDir) }
).match()
```

---

### "Process/workflow not found"

**Error:**
```
Process `FASTQC` not found in script
```

**Solutions:**

1. Check script path is correct:
```groovy
script "main.nf"  // Use the actual source selected for this test
```

2. Check process/workflow name matches exactly:
```groovy
process "FASTQC"  // Must match: process FASTQC { ... }
workflow "BAM_STATS"  // Must match: workflow BAM_STATS { ... }
```

3. Ensure the script is valid Nextflow:
```bash
nextflow lint
```

---

### "Input channel structure mismatch"

**Error:**
```
No such variable: meta
```

**Solutions:**

Check your input matches the process definition:

```groovy
// Process expects: tuple val(meta), path(reads)
process {
    """
    input[0] = [
        [ id:'test', single_end:true ],  // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)  // reads
    ]
    """
}

// Process expects: tuple val(meta), path(reads), path(index)
process {
    """
    input[0] = [
        [ id:'test' ],
        file('test.bam'),
        file('test.bam.bai')
    ]
    """
}
```

---

### "Out of memory"

**Error:**
```
Process exceeded memory limit
```

**Solutions:**

1. Use smaller test data (sarscov2 instead of homo_sapiens)

2. Override resources in test:
```groovy
test("memory intensive") {
    options "-process.memory=8.GB"
    // ...
}
```

3. Configure in test nextflow.config:
```groovy
// User-provided test configuration
process {
    memory = 4.GB
    cpus = 2
}
```

---

### "Test timeout"

**Error:**
```
Test execution timed out
```

**Solutions:**

1. Increase timeout in test:
```groovy
test("slow process") {
    options "-process.time=2h"
    // ...
}
```

2. Use smaller test data

3. Set global timeout in nf-test.config:
```groovy
config {
    timeout 3600  // seconds
}
```

---

### "Profile not found"

**Error:**
```
Unknown profile: docker
```

**Solutions:**

1. Ensure profile exists in nextflow.config:
```groovy
profiles {
    docker {
        docker.enabled = true
    }
    singularity {
        singularity.enabled = true
    }
}
```

2. Specify profile correctly:
```bash
nf-test test --profile docker
```

3. Set default in nf-test.config:
```groovy
config {
    profile "docker"
}
```

---

### "Container not found"

**Error:**
```
docker: Error response from daemon: manifest not found
```

**Solutions:**

1. Pull container manually:
```bash
docker pull quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
```

2. Check container URL in module meta.yml

3. Use different container registry:
```groovy
process {
    container = "biocontainers/fastqc:0.12.1"
}
```

---

### "Stub test not running stub"

**Error:** Process runs full script instead of stub

**Solutions:**

Enable stub mode:
```groovy
test("stub run") {
    options "-stub"
    // ...
}
```

Or set globally:
```bash
nf-test test --profile test --stub
```

---

## Debugging Workflow

Run nf-test with debug or verbose output, then inspect the task artifacts and
trace records identified by that output. Read .command.sh, .command.log, and
.command.err from the actual affected task, using locations reported by the
execution environment. Do not assume a work-directory layout.

Compare the real invocation, inputs, profile, container, and assertions before
making changes. Any direct task replay must use the same verified environment
and user-provided inputs.

## CI-Specific Issues

### CI Storage and Caching

Inspect the runner's actual disk use and the test's required storage. Use the
project's supported cleanup and caching strategy without deleting assumed system
folders. Cache only the verified resources the test needs, and keep cache keys
sensitive to the relevant inputs and dependency versions.

### Parallel Test Execution

Split the user's test suite into explicit CI jobs or a supported matrix. Pre-download plugins when their download mechanism cannot safely run in parallel.
