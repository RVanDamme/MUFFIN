<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# snakemake-nextflow-playbook

Practical, incremental flow for migrating Snakemake workflows to Nextflow.

## 1) Inventory

- **Generate the DAG** — run visualization commands to understand workflow structure before reading code:
  - `snakemake --rulegraph | dot -Tsvg > rulegraph.svg` — rule-level dependency graph (best starting point)
  - `snakemake --dag | dot -Tsvg > dag.svg` — full job-level DAG (shows per-sample fan-out)
  - `snakemake --filegraph | dot -Tsvg > filegraph.svg` — file-level dependencies (reveals implicit inputs)
  - `--rulegraph` and `--dag` also accept `mermaid-js` format: `snakemake --rulegraph mermaid-js`
- Find entry points: `Snakefile`, `*.smk`, included files
- Find wrappers/scripts/envs used by rules
- Find existing tests and fixtures
- Identify **implicit prerequisites**: files used as inputs but never produced by a rule (indexes, databases, genome dicts) — these need explicit Nextflow processes. The `--filegraph` output makes these visible.

## 2) Verify conda/container availability

Before running anything, confirm packages resolve on the target platform:

- `conda search -c bioconda <pkg>` for each tool — verify versions exist
- On `osx-arm64`: bioconda lags; prefer latest stable over old tutorial pins
- Multi-tool envs (bwa + samtools): verify all packages coexist in one solve
- Have container fallbacks ready (BioContainers/mulled images)

See `nextflow-conda-troubleshooting.md` reference for common errors and fixes.

## 3) Baseline run (Snakemake)

Run on a small deterministic dataset and capture:
- Produced files (manifest)
- Hashes/checksums for key artifacts
- Summary metrics (rows, sample counts, report stats)

Use this baseline as migration oracle.

## 4) Create Nextflow frame

- Build `main.nf` graph skeleton matching the rulegraph from step 1 — each rule node → process, each edge → channel
- Add `nextflow.config` params/resources/profile defaults
- Stub key processes so `nextflow run <pipeline> -preview` compiles early
- Add explicit processes for any implicit prerequisites found in step 1

## 5) Rule-by-rule conversion

For each Snakemake rule:
1. Search nf-core modules first
2. If match exists, install/reuse and patch behavior deltas
3. Else implement local module with explicit I/O and deterministic script
4. Convert `script:` blocks — replace `snakemake.input`/`output`/`params` with CLI args
5. Add/adjust `nf-test` coverage

## 6) Included files and wrappers

- `include: '*.smk'` blocks often map to subworkflows
- Wrapper commands map to modules (nf-core or local)
- Keep modules small; compose in workflows

## 7) Verification loop

After each chunk:
1. Run targeted `nf-test`
2. Run `nextflow run <pipeline> -preview`
3. Run `nextflow lint`
4. Run fixture pipeline and compare against Snakemake baseline

## 8) Equivalence checks

Compare migrated outputs to baseline using:
- Filename-level parity for required outputs
- Checksums for deterministic artifacts
- Metric tolerances for non-deterministic steps

If mismatch: add regression test first, then patch module/workflow.
