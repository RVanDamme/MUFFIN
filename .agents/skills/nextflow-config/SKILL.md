---
name: nextflow-config
description: >
  Generate and explain Nextflow configuration files. Use when user asks to create,
  modify, debug, or understand nextflow.config files. Covers all config scopes,
  process selectors, profiles, executor settings, and container runtimes.
---
<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Nextflow Config Skill

Generate, explain, and troubleshoot `nextflow.config` files.

## When to Use

Load this skill when the user wants to:
- Generate a `nextflow.config` for their pipeline or environment
- Configure executors (local, slurm, sge, aws batch, google batch, etc.)
- Set up container runtimes (docker, singularity, apptainer, wave)
- Define process resource requirements (cpus, memory, time)
- Create config profiles (test, production, cloud, hpc)
- Debug config issues or understand config precedence
- Configure Nextflow Platform (tower) integration
- Set up retry strategies and error handling

## References

Use the host's file-reading tool to load detailed docs:
- `config-guide.md` — Config syntax, includes, profiles, process selectors
- `config-options.md` — Complete reference for all config scopes and options

## Config File Basics

### Precedence (lowest → highest)
1. User-level Nextflow configuration
2. `nextflow.config` in project dir
3. `nextflow.config` in launch dir
4. `-c <file>` CLI option
5. `-profile` selections
6. `--param=value` CLI overrides

### Syntax

```groovy
// Dot syntax
process.executor = 'slurm'
process.queue = 'long'

// Block syntax (equivalent)
process {
    executor = 'slurm'
    queue = 'long'
}

// Include other configs
includeConfig params.additional_config
```

## Common Patterns

### Minimal Local Config

```groovy
params {
    input  = null
    outdir = null  // user supplies an output destination
}

process {
    cpus   = 2
    memory = 4.GB
    time   = 1.h
}

docker.enabled = true
```

### HPC (SLURM) Config

```groovy
process {
    executor = 'slurm'
    queue    = 'normal'
    cpus     = 4
    memory   = 8.GB
    time     = 4.h

    withLabel: big_mem {
        cpus   = 16
        memory = 64.GB
        queue  = 'highmem'
    }

    withLabel: gpu {
        queue       = 'gpu'
        clusterOptions = '--gres=gpu:1'
    }
}
```

### AWS Batch Config

```groovy
process {
    executor = 'awsbatch'
    queue    = 'my-batch-queue'
}

aws {
    region = 'us-east-1'
    batch {
        // Set cliPath only if the user's AWS CLI cannot be found automatically.
    }
}

// Set workDir to the work storage selected for the user's execution environment.
```

### Nextflow Platform (Tower) Config

```groovy
tower {
    enabled   = true
    endpoint  = 'https://api.cloud.seqera.io'
    accessToken = secrets.TOWER_ACCESS_TOKEN
}
```

### Multi-Profile Config

```groovy
profiles {
    standard {
        process.executor = 'local'
        docker.enabled = true
    }

    slurm {
        process.executor = 'slurm'
        process.queue = 'normal'
        singularity.enabled = true
    }

    test {
        params.input  = null  // user selects representative test data
        params.outdir = null  // user selects a test output destination
        process.cpus  = 1
        process.memory = 2.GB
    }
}
```

### Process Selectors

Selectors apply settings to specific processes. Priority (lowest → highest):
1. Global process settings (no selector)
2. Process directives in process definition
3. `withLabel:` matching any label
4. `withName:` matching process name
5. `withName:` matching included alias
6. `withName:` matching fully qualified name (`workflow:process`)

```groovy
process {
    cpus = 2

    withLabel: 'cpu_high' {
        cpus = 16
    }

    withName: 'FASTQC' {
        cpus   = 4
        memory = 8.GB
    }

    // Regex patterns
    withName: 'ALIGN_.*' {
        cpus   = 8
        memory = 32.GB
    }

    // Negation
    withLabel: '!small' {
        memory = 16.GB
    }
}
```

### Retry Strategy

```groovy
process {
    errorStrategy = 'retry'
    maxRetries    = 3
    maxErrors     = '-1'

    // Dynamic resources on retry
    memory = { 8.GB * task.attempt }
    time   = { 4.h * task.attempt }
}
```

### `errorStrategy` placement

Keep `errorStrategy` in config, not in module files. A module-level directive can be overridden via `process { withName: 'FOO' { errorStrategy = ... } }`, but spreading retry policy across many modules makes per-site tuning harder to reason about. Use one global policy in the pipeline's shared configuration:

```groovy
// Shared process configuration
process {
    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'terminate' }
    maxRetries    = 1
}
```

These are Linux signal-derived exit codes: 130 (SIGINT), 137 (SIGKILL / OOM), 143 (SIGTERM), 104 (ECONNRESET / preempted spot instance). Many bioinformatics tools return `1` for all failure modes — this pattern intentionally retries only on signal-derived codes to avoid masking deterministic tool failures. Declare `errorStrategy` in a specific module only when it has a well-understood failure mode genuinely different from the global policy (e.g. a download step that should `ignore` 404s), and document why.

### `resourceLimits`: cap requests to available infrastructure

The `resourceLimits` directive (Nextflow 24.04+) caps resource requests at specified maximums. If a retry expression requests more than the limit, Nextflow silently caps it — the task still runs at the limit. Without this, `memory { 8.GB * task.attempt }` on attempt 9 requests 72 GB and may exceed the largest available node, causing a scheduling failure.

```groovy
// nextflow.config — global limits matching your infrastructure
process {
    resourceLimits = [
        cpus: 72,
        memory: 512.GB,
        time: 168.h
    ]
}
```

Use profiles to vary limits per environment:

```groovy
profiles {
    local {
        process.resourceLimits = [cpus: 8, memory: 16.GB, time: 24.h]
    }
    cloud {
        process.resourceLimits = [cpus: 72, memory: 512.GB, time: 168.h]
    }
}
```

Set to match the largest instance type (or node class) in your compute environment. Always declare in config, not in module files.

### Container Runtimes

```groovy
// Docker
docker {
    enabled    = true
    runOptions = '-u $(id -u):$(id -g)'
}

// Singularity
singularity {
    enabled   = true
    autoMounts = true
    // Set cacheDir only when a user-selected writable cache is required.
}

// Apptainer
apptainer {
    enabled   = true
    autoMounts = true
    // Set cacheDir only when a user-selected writable cache is required.
}

// Wave (Seqera)
wave {
    enabled = true
    strategy = ['conda', 'container', 'dockerfile', 'spack']
}
```

## Generation Guidelines

When generating a config:

1. **Ask about the execution environment** — local, HPC (slurm/sge/pbs), cloud (aws/gcp/azure)?
2. **Ask about container runtime** — docker, singularity, apptainer, conda?
3. **Ask about resource needs** — typical cpu/memory/time for processes?
4. **Ask about Nextflow Platform** — tower integration needed?
5. **Use profiles** — separate environments into profiles (test, dev, production)
6. **Use labels** — group processes by resource needs, not individual withName selectors
7. **Dynamic resources** — use `task.attempt` multiplier for retry strategies
8. **Prefer block syntax** — more readable than dot syntax for complex configs
9. **Include comments** — explain non-obvious settings
10. **Follow nf-core conventions** when the pipeline uses nf-core patterns

## Key Config Scopes

| Scope | Purpose |
|-------|---------|
| `params` | Pipeline parameters |
| `process` | Process directives (executor, resources, containers) |
| `executor` | Executor settings (queue size, poll interval) |
| `docker` | Docker runtime settings |
| `singularity` | Singularity runtime settings |
| `apptainer` | Apptainer runtime settings |
| `podman` | Podman runtime settings |
| `conda` | Conda environment settings |
| `aws` | AWS/S3/Batch settings |
| `azure` | Azure Batch/Storage settings |
| `google` | Google Cloud Batch/Storage settings |
| `wave` | Seqera Wave container settings |
| `tower` | Nextflow Platform settings |
| `mail` | Email notification settings |
| `manifest` | Pipeline metadata |
| `report` / `timeline` / `trace` / `dag` | Execution reports |
| `profiles` | Named config profiles |

For the complete list of options in each scope, load the `config-options.md` reference.
