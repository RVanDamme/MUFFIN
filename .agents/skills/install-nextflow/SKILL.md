---
name: install-nextflow
description: Install or upgrade Nextflow on the user's machine. Use when the user wants to install Nextflow, check their current Nextflow version, upgrade Nextflow, or set up the prerequisites (Java 17+) needed to run Nextflow.
---
<!-- Modified for the Nextflow plugin: generic host tools and companion guidance. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](../launch-workflow/references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Install or Upgrade Nextflow

Install Nextflow (or upgrade an existing installation) and ensure the Java prerequisite is in place.

## Step 1: Check Whether Nextflow Is Installed

```bash
command -v nextflow && nextflow -version
```

- **If installed**, parse the version from the output and proceed to Step 2.
- **If not installed**, skip to Step 3.

## Step 2: Check for Upgrades (Existing Install)

If Nextflow is already installed:

1. Read the installed version from `nextflow -version`.
2. If the installed version is **older than 26.04**, the user MUST upgrade — other Nextflow skills depend on it.
3. If the user needs **workflow modules** (registry modules with `Kind: Workflow`), the installed version must be **26.10 or later**.
4. Otherwise, propose an upgrade only if a newer version is available.

Run a self-update (it checks the latest release and updates if newer):

```bash
nextflow self-update
```

Confirm with the user before running it. After the update, re-run `nextflow -version` to verify.

If the install is already up to date, stop here.

## Step 3: Verify Java 17+ (Required for Installing Nextflow)

Nextflow requires Java 17 or later. Check the installed Java version:

```bash
java -version
```

Note: `java -version` writes to **stderr**, not stdout. Capture both streams when parsing:

```bash
java -version 2>&1 | head -1
```

- If Java 17+ is present, proceed to Step 4.
- If Java is missing or older than 17, install Java 21 using SDKMAN (Step 3a).

### Step 3a: Install Java 21 via SDKMAN

First, check whether SDKMAN is installed:

```bash
command -v sdk
```

If SDKMAN is missing, install it:

```bash
curl -s "https://get.sdkman.io" | bash
```

Activate SDKMAN using its installer guidance for the user's shell, then install Java 21:

```bash
sdk install java 21-tem
```

After installation, verify:

```bash
java -version
```

## Step 4: Install Nextflow

Prefer `curl`; fall back to `wget` if `curl` is unavailable.

Detect which downloader is available:

```bash
command -v curl
command -v wget
```

**Using curl:**

```bash
curl -fsSL https://get.nextflow.io | bash
```

**Using wget (if curl is missing):**

```bash
wget -qO- https://get.nextflow.io | bash
```

If neither `curl` nor `wget` is available, stop and ask the user to install one of them.

Make the installed launcher executable and available through the user's selected
installation method and executable search configuration. Respect the operating
system, permissions, and destination chosen by the user; do not prescribe an
installation directory.

## Step 5: Verify the Installation

```bash
nextflow -version
```

Confirm the reported version is **26.04 or later**. Present the version to the user.

## Nextflow Versioning

### Calendar versioning

Stable Nextflow releases use **calendar versioning** in the form `YY.MM.PATCH`:

- `25.10.1` — patch 1 of the October 2025 release
- `26.04.0` — first release of April 2026

Compare versions as calendar dates, not semver — `26.04.0` is newer than `25.10.1`.

### Edge releases

Edge (preview) releases append the `-edge` suffix, e.g. `26.05.0-edge`. They include unreleased features and are not recommended for production.

Switch the active channel via the `NXF_EDGE` environment variable:

```bash
export NXF_EDGE=1   # use edge releases
export NXF_EDGE=0   # back to stable releases (default)
```

This affects which release `nextflow self-update` and the launcher resolve to.

### Pinning a specific version

Pin any version (stable or edge) with `NXF_VER`. The launcher will download and use that exact version on the next run:

```bash
export NXF_VER=25.10.1        # pin a stable release
export NXF_VER=26.05.0-edge   # pin an edge release
```

Unset `NXF_VER` (`unset NXF_VER`) to return to the latest version of the active channel.

When the user asks to install or run a specific version, set `NXF_VER` rather than reinstalling.

## Critical Rules

1. **CHECK FIRST** — always run `nextflow -version` before deciding to install or upgrade.
2. **REQUIRE 26.04+** — other Nextflow skills depend on `nextflow module`, `nextflow auth`, and `nextflow launch`, which require this version. Workflow modules require 26.10 or later.
3. **VERIFY JAVA 17+** before installing Nextflow — install Java 21 via SDKMAN if missing or too old.
4. **CONFIRM destructive actions** — ask the user before running `nextflow self-update`, installing SDKMAN, or moving the `nextflow` launcher to a system directory.
5. **PREFER curl, FALL BACK to wget** — if neither is available, stop and ask the user to install one.
6. **REMEMBER stderr** — `java -version` writes to stderr; redirect with `2>&1` when parsing.
7. **VERIFY at the end** — always run `nextflow -version` after install/upgrade to confirm the result.
8. **USE NXF_VER for pinning** — when the user asks for a specific Nextflow version, set `NXF_VER` instead of reinstalling.
9. **TREAT versions as calendar dates** — Nextflow uses `YY.MM.PATCH`, not semver, so `26.04.0` is newer than `25.10.1`.
