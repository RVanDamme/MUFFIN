<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# nft-utils (v0.0.9) — Utility Functions & Path Extensions

A collection of utility functions for nf-test, primarily developed by the nf-core community. Provides snapshot helpers, file collection utilities, dependency management, and output sanitization.

Variables beginning with user in the examples represent files, selections, or
locations supplied by the user or reported by the test environment. Establish
those values before invoking a helper; do not infer fixed folders from examples.

**Source:** https://github.com/nf-core/nft-utils
**Docs:** https://nf-co.re/nft-utils

### Path Extensions

```groovy
// MD5 checksum
assert path(userTextFile).md5 == "d41d8cd98f00b204e9800998ecf8427e"

// JSON parsing
def json = path(userJsonFile).json
assert json.version == "1.0"

// YAML parsing
def yaml = path(userYamlFile).yaml
assert yaml.params != null

// Line-by-line reading
assert path(userTextFile).lines.any { it.contains("PASS") }

// Gzipped file lines
assert path(userCompressedFile).linesGzip.size() > 100
```

### Snapshot Functions

#### `removeNextflowVersion(filePath)`

Removes the Nextflow version entry from nf-core software versions YAML files. Useful because the Nextflow version changes between environments and should not be part of snapshots.

Use the user's verified versions file or matching-file selection with the helper. Do not assume its location.

**Arguments:** Path to a YAML versions file (supports `*` and `?` wildcards). When using wildcards, all matching files are merged.

**Note:** Returned YAML keys are sorted alphabetically for consistent output.

#### `removeFromYamlMap(filePath, sectionKey [, subKey])`

Remove any key or entire section from a YAML file. Supports two patterns:

~~~groovy
// Remove a subkey from a selected section
assert snapshot(removeFromYamlMap(userYamlFile, sectionName, subkeyName)).match()

// Remove a selected section
assert snapshot(removeFromYamlMap(userYamlFile, sectionName)).match()
~~~

Use the user's verified YAML file and actual key names. Do not assume its location.

**Arguments:**
- First: Path to YAML file (supports wildcards)
- Second: Top-level key (section name)
- Third (optional): Subkey to remove. If omitted, entire section is removed.

#### `getAllFilesFromDir()`

Generates a list of all files within a directory (and subdirectories), with filtering support via glob patterns.

**Important:** Requires absolute paths. Use the actual location supplied by the user or reported by nf-test; do not invent an output folder.

~~~groovy
def files = getAllFilesFromDir(userOutputLocation,
    includeDir: false,
    ignore: userIgnorePatterns,
    ignoreFile: userIgnoreRules,
    relative: true)
assert snapshot(files).match()
~~~

**Named parameters** (preferred for clarity):

The positional arguments are the selected location, whether to include directory
names, exclusion patterns, an optional file of additional exclusions, inclusion
patterns, and whether to return relative names. Named equivalents are
includeDir, ignore, ignoreFile, include, and relative. Supply patterns and
selection files from the user's actual test contract.

#### `getRelativePath(fileList, baseDir)`

Converts a list of absolute file paths to paths relative to a base directory. Useful for readable snapshots.

~~~groovy
def files = getAllFilesFromDir(userOutputLocation)
assert snapshot(getRelativePath(files, userBaseLocation)).match()
~~~

**Tip:** You can also use `getAllFilesFromDir()` with `relative: true` instead of calling `getRelativePath()` separately.

#### `getAllFilesFromChannel(channelOutput)`

Extracts absolute file paths from Nextflow channel outputs. Handles nested structures, filters out metadata maps, and returns only file path strings.

```groovy
assert snapshot(
    getAllFilesFromChannel(process.out.html),
    getAllFilesFromChannel(process.out.zip)
).match()

// Combine with .collect() for file names only
assert snapshot(
    getAllFilesFromChannel(process.out.html).collect { f -> file(f).name }
).match()
```

#### `listToMD5(list)`

Converts a list of values to a single MD5 hash. Useful for creating a stable hash after filtering unstable lines from a file.

```groovy
def stableLines = path(userTextFile).lines.findAll { !it.contains("timestamp") }
assert listToMD5(stableLines) == "expected_md5_hash"
```

#### `filterNextflowOutput(output)`

Filters Nextflow stdout/stderr to remove variable content (timestamps, process hashes, paths, versions) for stable snapshots. Lines are sorted where appropriate for reproducibility.

```groovy
// Basic usage with stdout, stderr, or both
def filtered = filterNextflowOutput(workflow.stdout + workflow.stderr)

assert snapshot(
    filterNextflowOutput(workflow.stdout + workflow.stderr)
).match()
```

**Named parameters:**

```groovy
// Keep ANSI escape codes (stripped by default)
filterNextflowOutput(workflow.stdout, keepAnsi: true)

// Ignore lines containing specific strings
filterNextflowOutput(workflow.stdout, ignore: ["Submitted process"])

// Include only lines containing specific strings
filterNextflowOutput(workflow.stdout, include: ["Submitted process"])

// Disable sorting (not recommended — may cause flaky snapshots)
filterNextflowOutput(workflow.stdout, sorted: false)
```

**Automatically filtered patterns:**
- Empty/whitespace-only lines — removed
- Timestamps — replaced with `[TIMESTAMP]`
- Process hashes (e.g., `[57/0d391c]`) — replaced with `[PROCESS_HASH]`
- Absolute paths — replaced with `[PATH]`
- Nextflow version messages — replaced with `[VERSION]` or removed
- `Staging foreign file`, `Submitted process`, `WARN:`, `ERROR:` lines — sorted alphabetically

#### `sanitizeOutput(processOutput [, unstableKeys: [...]])`

Cleans process/workflow outputs by removing numbered keys for more readable snapshots.

```groovy
// Basic usage
assert snapshot(sanitizeOutput(process.out)).match()

// Mark keys with unstable file content (snapshots name only, not md5)
assert snapshot(sanitizeOutput(process.out, unstableKeys: ["zip"])).match()
```

#### `curlAndExtract(url, destDir [, format])`

Downloads an archive from the Internet and extracts it. Supports zip and tar archives (with gzip, bzip2, xz, lz4, lzma, lzop, zstd compression). Format is auto-detected from the URL unless explicitly provided.

```groovy
setup {
    curlAndExtract("https://example.com/database.zip", userExtractionLocation)
    curlAndExtract("https://example.com/db.tar.gz", userExtractionLocation)
    // Explicit format when URL doesn't indicate it
    curlAndExtract("https://example.com/secret/data", userExtractionLocation, "tar.bz2")
}

when {
    params {
        db_path = userDatabaseLocation
    }
}

cleanup {
    new File(userExtractionLocation).deleteDir()
}
```

### Dependency Management (nf-core modules)

These functions manage dependencies on nf-core components for tests in non-nf-core repositories.

#### `nfcoreInitialise(libraryPath)`

Set up a temporary nf-core module library directory.

```groovy
setup {
    nfcoreInitialise(userLibraryLocation)
}
```

#### `nfcoreInstall(libraryPath, modules)`

Install nf-core modules into the temporary library. Accepts a list of module name strings or maps with `name`, `sha` (optional), and `remote` (optional).

```groovy
setup {
    nfcoreInitialise(userLibraryLocation)

    // Simple list of module names
    nfcoreInstall(userLibraryLocation, ["minimap2/index"])

    // With specific SHA and/or custom remote
    nfcoreInstall(userLibraryLocation, [
        [name: "minimap2/align", sha: "5850432aab24a1924389b660adfee3809d3e60a9"],
        [name: "fastqc", remote: "https://github.com/nf-core-test/modules.git"],
    ])
}
```

Tracks installed state to skip redundant re-installs.

#### `nfcoreLink(libraryPath, modulesDir)`

Symlink a temporary library into your modules directory so modules can be referenced normally.

```groovy
setup {
    nfcoreInitialise(userLibraryLocation)
    nfcoreInstall(userLibraryLocation, ["minimap2/index", "minimap2/align"])
    nfcoreLink(userLibraryLocation, userModuleLocation)
}
```

#### `nfcoreUnlink(libraryPath, modulesDir)`

Remove symlinks created by `nfcoreLink()`. Use in `cleanup` block.

```groovy
cleanup {
    nfcoreUnlink(userLibraryLocation, userModuleLocation)
}
```

#### `nfcoreDeleteLibrary(libraryPath)`

Completely delete the temporary library directory.

```groovy
cleanup {
    nfcoreDeleteLibrary(userLibraryLocation)
}
```
