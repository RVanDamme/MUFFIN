<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# Auditing tool and helper-code availability

Two questions per dependency: **was the code provided?** and **can the tool be
containerized and run?** A conversion is not startable until both are answered
for every step.

## 1. Resolve helper code against what you were given

Build the list of code the source calls, then check each entry exists in the
material you actually have.

| Source pattern | Resolves to |
|---|---|
| A Python helper invocation | The referenced helper must be provided |
| An R helper invocation | Check the user-provided helper and its dependencies |
| `source("helpers.R")`, `import mylab.utils` | Sibling module or an internal package |
| A modified import search location | An external library that must be available |
| Snakemake `script:` / `wrapper:` | Helper script or a pinned wrapper repo revision |
| A Java JAR invocation | A binary artifact rather than a Python or Conda package |
| A compiled or vendored helper | Required binary or source artifact |
| A versioned cluster module load | Required tool and validated version |

Anything that does not resolve is a blocker. Say precisely what is missing and
what you need: the file, the internal package name and version, or permission to
reimplement the logic from the description.

Two specifics worth calling out when you see them:

- **Internal/private packages** — an `import` of a package that is not on PyPI,
  CRAN, Bioconductor, or conda needs either the source or a wheel. Note that it
  will also need to be installed into whatever container the process uses.
- **Vendored binaries with no source** — record the exact filename and
  architecture. A Linux container cannot run a macOS-only binary, and this
  surfaces late and confusingly if you skip it.

## 2. Check external tool availability

For each bare command the source invokes (`bwa`, `samtools`, `salmon`, `gatk`,
custom vendor CLIs), work in this order and stop as soon as you have an answer:

1. **Existing nf-core module** — `search_nfcore_module` / `describe_nfcore_module`.
   A module means the container, versions topic, and stub are already solved. This
   is the cheapest possible outcome, so check it first.
2. **Bioconda / biocontainers** — a bioconda recipe implies a
   Biocontainers image. Note the recipe name, since it is often not
   the command name.
3. **Official upstream image** — vendor-published images on Docker Hub, quay.io,
   or `nvcr.io`. Prefer a version-pinned tag over `latest`.
4. **Build required** — nothing usable exists. Hand off to `create-container`
   rather than assuming a Dockerfile is trivial.

In the report, name the image (or the absence of one) per tool. "Container: needs
build" is a schedule fact the user should see during the audit, not later.

## 3. Categories that block a conversion outright

Flag these explicitly, even when a container exists — they change the plan, not
just the Dockerfile:

- **License-gated** — commercial or academic-only tools that cannot ship in a
  public image (some aligners, structural-variant callers, and vendor secondary
  analysis suites). Ask how the license is provisioned before designing around
  the tool.
- **Registration-gated downloads** — reference bundles and databases behind a
  login (certain annotation databases, some model checkpoints). The pipeline
  cannot fetch these unattended; they become a staged input or a data link.
- **GPU-only** — a GPU requirement constrains the compute environment. Record the
  CUDA/driver expectation. If the target environment has no GPU, route to
  `find-alternative-tools` for a CPU path.
- **Interactive or GUI-only** — cannot be a pipeline step at all. It needs either
  a headless mode or removal from scope.
- **Abandoned with no image and no source release** — treat like a missing tool
  and offer alternatives.

## 4. Versions

Unpinned versions are a silent reproducibility gap, not a detail to settle later.
Collect whatever evidence the source gives you — `module load` lines, lockfiles
(`environment.yml`, `requirements.txt`, `renv.lock`, `Pipfile.lock`), container
tags, `--version` output captured in logs, methods sections — and for anything
still unresolved, ask once.

Two version findings that are worth surfacing on their own:

- A tool whose output format changed across major versions (VCF-producing callers,
  annotation tools). Getting this wrong makes downstream steps fail confusingly.
- A language runtime with a hard floor (`python>=3.10` syntax, R package that
  needs R 4.3+). This constrains the base image.

## 5. What good output looks like

One row per dependency, sorted with blockers first:

| Dependency | Kind | Provided? | Version | Container / module | Status |
|---|---|---|---|---|---|
| User-referenced helper | helper script | no | — | n/a | **blocker** — not in handover |
| `mylab.qc` | internal package | no | unknown | needs install into image | **blocker** — need source or wheel |
| `samtools` | external tool | n/a | unpinned | nf-core module `samtools/sort` | ok — confirm version |
| `novoalign` | external tool | n/a | 4.03 | no public image | **blocker** — license |
| `picard.jar` | vendored binary | yes | 2.27.5 | biocontainers image available | ok |
