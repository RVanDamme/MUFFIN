---
name: launch-workflow
description: Launch Nextflow pipeline executions on cloud and HPC clusters via Nextflow Platform. Use when the user wants to run/launch/submit a pipeline on a cloud or cluster compute environment, configure a compute environment, push pipeline changes to GitHub before launching, or sign in to Nextflow Platform.
---
<!-- Modified for the Nextflow plugin: generic host tools and companion guidance. -->

## Using this skill

Use the tools and runtimes available in the user's session. Establish the
pipeline, input data and execution environment from their request or current
context. Use user-selected locations for files and results. Consult companion
skills, bundled references and optional helpers by name when available.

For Seqera operations, read the [MCP connection](references/seqera-mcp/README.md) connection reference. Any `platform_*` examples
describe API operations. Discover their exact names and parameter schemas with
`search_seqera_api`, then invoke `call_seqera_api` using those schemas. The host
manages OAuth for the connected MCP server.


# Launch Pipeline on Nextflow Platform

Launch Nextflow pipeline executions on cloud (AWS, Google Cloud, Azure) and HPC clusters (Slurm, LSF, etc.) through Nextflow Platform. Nextflow Platform manages the target compute environment, executes the pipeline from a remote Git repository, and provides monitoring.

**Requires Nextflow 26.04 or later** (for the `nextflow auth` and `nextflow launch` commands).

## Main Flow

1. **Configure the compute environment** — select a Nextflow Platform compute environment for the target cloud or cluster.
2. **Ensure the pipeline is in a remote Git repository** — Nextflow Platform launches pipelines from GitHub (or a compatible Git host such as GitLab or Bitbucket). If the local pipeline is not yet hosted, assist the user in setting up the repository.
3. **Upload local changes** — push any uncommitted local changes to the remote so the launched run reflects the user's current code.
4. **Launch with `nextflow launch`** — submit the pipeline by passing the Git repository URL and the expected parameters.

## Step 1: Authenticate with Nextflow Platform

Check whether the user is already authenticated:

```bash
nextflow auth status
```

If not authenticated, sign in:

```bash
nextflow auth login
```

This opens a browser for authentication. Wait for the user to complete the login flow before proceeding.

## Step 2: Configure the Compute Environment

List available compute environments using the Seqera API:

```
call_seqera_api(
  service: "platform",
  api_name: "platform_list_compute_envs",
  parameters: {}
)
```

Selecting the compute environment:

1. Present the available compute environments to the user (name, platform type, region/cluster).
2. If multiple environments are available, ask the user which one to use.
3. If only one is available, propose it and ask for confirmation.
4. Confirm the selected compute environment with the user before launching.

## Step 3: Ensure the Pipeline Is in a Remote Git Repository

Nextflow Platform launches pipelines from a remote Git URL — it cannot launch directly from a local path. Verify the pipeline directory is a Git repository connected to a remote on GitHub (or a compatible host like GitLab or Bitbucket).

Check the current state:

```bash
git remote -v
git status
```

**If the pipeline is not yet a Git repository or has no remote**, assist the user in setting it up:

Initialize Git if needed, review the user's selected pipeline files, and stage
only the intended changes. Help the user create or select a remote repository,
configure its remote URL, and push the agreed revision using the available tools.

## Step 4: Upload Local Changes

Before launching, push any local changes so the remote reflects the code that should run:

Review the changes, commit the user's intended pipeline files with a descriptive
message, and push the confirmed branch using the available Git tools.

If the working tree is clean and the local branch is in sync with the remote, skip this step.

## Step 5: Launch the Pipeline

Use `nextflow launch` with the user-confirmed public or private Git repository
URL, a pinned revision, and the pipeline parameters from the user's launch
configuration. Obtain input and output destinations from the user rather than
providing example storage locations. Check required parameters before submitting.

A public repository URL such as https://github.com/nf-core/rnaseq identifies the
pipeline; it does not determine the user's inputs or output destination.

### Key options

| Option | Description |
|--------|-------------|
| `-r <revision>` | Pipeline version, branch, or tag |
| `--input` | Input samplesheet or data |
| `--outdir` | Output directory (cloud storage when running on cloud CEs) |
| `-params-file` | User-provided parameter configuration |
| `-resume` | Resume a previous execution |

## Step 6: Monitor the Execution

`nextflow launch` returns a run URL on Nextflow Platform. Present this URL to the user so they can monitor progress in the Platform UI.

## Critical Rules

1. **AUTHENTICATE first** — check `nextflow auth status` before attempting to launch.
2. **CONFIRM the compute environment** — always show the user the selected CE before launching.
3. **REQUIRE a remote Git repository** — the pipeline must be hosted on GitHub or a compatible Git host. If not, help the user set it up before launching.
4. **PUSH local changes first** — the launched run uses the remote code, so local edits must be committed and pushed.
5. **PASS the Git repository URL** to `nextflow launch`, not a local path.
6. **PIN the revision** — always use `-r` to target a specific branch, tag, or commit for reproducibility.
7. **USE cloud storage for `--outdir`** when the target CE runs on a cloud platform (s3://, gs://, az://).
8. **PRESENT the run URL** returned by `nextflow launch` so the user can monitor the run.

## MCP connection reference

When a Seqera connection is missing or operation discovery is unclear, use the shared connection reference. Otherwise rely on the live tool descriptions and returned schemas. Read [MCP connection](references/seqera-mcp/README.md) before proceeding.
