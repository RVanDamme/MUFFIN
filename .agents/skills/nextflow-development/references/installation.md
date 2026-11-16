<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Installation

## Contents
- [Quick install](#quick-install)
- [Docker setup](#docker-setup)
- [Singularity setup (HPC)](#singularity-setup-hpc)
- [nf-core tools (optional)](#nf-core-tools-optional)
- [Verify installation](#verify-installation)
- [Common issues](#common-issues)

## Quick install

Install Nextflow using its [official installation guidance](https://docs.seqera.io/nextflow/install).
Choose a supported Java runtime and make the commands available through the user's
normal command search settings.

```bash
nextflow -version
java -version
```

## Docker setup

### Linux
```bash
sudo apt-get update && sudo apt-get install docker.io
sudo systemctl enable --now docker
sudo usermod -aG docker $USER
# Log out and back in
```

### macOS
Download Docker Desktop: https://docker.com/products/docker-desktop

### Verify
```bash
docker run hello-world
```

## Singularity setup (HPC)

```bash
# Ubuntu/Debian
sudo apt-get install singularity-container

# Or via conda
conda install -c conda-forge singularity
```

### Configure cache

If the runtime requires `NXF_SINGULARITY_CACHEDIR`, set it to a user-selected
writable cache location appropriate for the execution environment.

## nf-core tools (optional)

```bash
pip install nf-core
```

Useful commands:
```bash
nf-core list                    # Available pipelines
nf-core launch rnaseq           # Interactive parameter selection
nf-core download rnaseq -r 3.14.0  # Download for offline use
```

## Verify installation

Run the selected pipeline's small test profile with a user-selected output
destination, then confirm its expected reports and completion status.

## Common issues

**Java version wrong:**
Select the installed supported JDK through `JAVA_HOME` or the user's normal
runtime-management mechanism.

**Docker permission denied:**
```bash
sudo usermod -aG docker $USER
# Log out and back in
```

**Nextflow not found:**
Make the installed Nextflow command available through the user's normal shell
or runtime-management settings.
