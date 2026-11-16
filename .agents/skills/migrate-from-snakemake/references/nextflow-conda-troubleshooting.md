<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
<!-- Purpose: Conda/container resolution failures during Snakemake→Nextflow migration. -->

# nextflow-conda-troubleshooting

Quick reference for diagnosing conda and container failures in migrated pipelines.

## PackagesNotFoundError

Most common failure when copying Snakemake conda pins to Nextflow.

**Diagnosis:**
```bash
conda search -c bioconda <package>          # list available versions
conda search -c bioconda <package>=<ver>    # check specific version
```

**Common causes:**
- Old version not built for target platform (especially `osx-arm64`)
- Typo in package name or channel
- Package renamed or moved to a different channel

**Fix:** bump to the latest stable version available for your platform.

## Platform-specific availability

Bioconda builds lag on ARM64 (Apple Silicon). If a version is missing:

1. Try the latest available version: `conda search -c bioconda <pkg> | tail -5`
2. Use a container instead (containers run x86 images via emulation on ARM64)
3. Use Wave containers via `-with-wave` for automatic conda-to-container conversion

## Nextflow conda directive format

```nextflow
process FOO {
    conda 'bioconda::bwa=0.7.19 bioconda::samtools=1.21'
    // or
    conda params.tool_environment  // user-selected environment specification
}
```

Key differences from Snakemake:
- No `name:` field — Nextflow manages env names via content hash
- Multiple packages in one string, space-separated
- Channel prefix required (`bioconda::`, `conda-forge::`)

## Container alternatives

When conda fails, switch to containers:

```nextflow
process FOO {
    container 'biocontainers/bwa:0.7.19--h577a1c6_0'
    // or mulled image for multi-tool processes:
    container 'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:...'
}
```

Use profiles to offer both:
```nextflow
profiles {
    conda  { conda.enabled = true  }
    docker { docker.enabled = true }
}
```

## Common error patterns

| Error | Cause | Fix |
|-------|-------|-----|
| `PackagesNotFoundError` | Version not available for platform | Bump version or use container |
| `ResolvePackageNotFound` | Package doesn't exist in channel | Check channel and package name |
| `UnsatisfiableError` | Version conflict between packages | Split into separate conda envs per process |
| `CondaEnvException` | Corrupt or stale env cache | Identify and clear the affected environment cache with the user's authorization, then rerun |
