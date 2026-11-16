# Running MUFFIN on HPC2N Kebnekaise

The `hpc2n` profile submits every pipeline step as its own Slurm job and runs the tools in Apptainer containers. Nextflow itself is a small, long-lived "head" process that only submits and watches those jobs.

## One-time setup

```bash
# Java 17+ (pick a version listed by `module spider Java`) and Nextflow >= 26.04
module load Java/<version>
mkdir -p ~/.local/bin && curl -fsSL https://get.nextflow.io -o ~/.local/bin/nextflow && chmod +x ~/.local/bin/nextflow
```

Add to `~/.bashrc` (replace the two project names):

```bash
export PATH=$HOME/.local/bin:$PATH
export NXF_VER=26.10.6
export SBATCH_ACCOUNT=hpc2nXXXX-YYY                                        # compute project, charged for the jobs
export NXF_APPTAINER_CACHEDIR=/proj/nobackup/<storage-project>/apptainer_cache   # container images, downloaded once
```

`projinfo` lists your compute projects. Apptainer needs no module on Kebnekaise.

## Launching a run

Work in project storage (`/proj/nobackup/...`), not in your home directory: the `work/` directory gets large.

Save as `run_muffin.sh` next to your samplesheet:

```bash
#!/bin/bash
#SBATCH -A hpc2nXXXX-YYY
#SBATCH -J muffin_head
#SBATCH -n 1
#SBATCH --time=168:00:00
#SBATCH -o muffin_head_%j.log

module load Java/<version>
export PATH=$HOME/.local/bin:$PATH NXF_VER=26.10.6
export SBATCH_ACCOUNT=hpc2nXXXX-YYY
export NXF_APPTAINER_CACHEDIR=/proj/nobackup/<storage-project>/apptainer_cache

nextflow run /path/to/MUFFIN -profile hpc2n \
    --input samplesheet.csv \
    --outdir results \
    -resume
```

```bash
sbatch run_muffin.sh          # submit the head job
squeue -u $USER               # head job + the step jobs it submits
tail -f muffin_head_*.log     # pipeline progress
scancel <jobid>               # cancelling the head job stops the run; -resume continues it later
```

The head job needs one core and little memory but must live as long as the whole run (168 h is the Kebnekaise maximum). Submitting the same script again resumes from the cached steps.

Quick check that the profile works (small test data, a few minutes once analysis steps exist):

```bash
nextflow run /path/to/MUFFIN -profile test,hpc2n --outdir results_test
```

## What the profile sets

| Setting | Value | Why |
|---|---|---|
| Executor | Slurm, at most 50 jobs queued, 10 submissions per minute | be gentle with the scheduler |
| Account | from `SBATCH_ACCOUNT` | Slurm reads it itself; nothing project-specific is stored in the pipeline |
| Containers | Apptainer | available on all nodes |
| Per-job ceiling | 28 cores, 180 GB, 168 h | one regular Skylake node and the walltime cap |

## Not yet verified on Kebnekaise

- Memory requests: Nextflow asks Slurm for memory with `--mem`. Kebnekaise ties memory to cores (about 6.8 GB per Skylake core), so a job asking for much memory on few cores may be rejected or billed for more cores. The first real run will show; the fix, if needed, is one line in `conf/hpc2n.config`.
- Steps needing more than 180 GB (possibly GTDB-Tk or large assemblies) need the large-memory nodes (`-p largemem`), which require a separate allocation.
- No node type is forced; add `process.clusterOptions = '-C skylake'` in a custom config to pin one.
