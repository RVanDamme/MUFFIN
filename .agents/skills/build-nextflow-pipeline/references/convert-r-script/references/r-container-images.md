<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
# R Container Images and Conda Environments

## Rocker images (Docker)

| Image | Contents | Use when |
|-------|----------|----------|
| `rocker/r-ver:4.3` | Base R only | Minimal scripts, no tidyverse |
| `rocker/tidyverse:4.3` | R + tidyverse + devtools | General data analysis |
| `rocker/verse:4.3` | tidyverse + TeX | Scripts producing PDF/LaTeX |
| `rocker/rstudio:4.3` | RStudio Server | Interactive development |
| `rocker/ml:4.3` | tidyverse + ML packages | Machine learning |

## Bioconductor images (Docker)

| Image | Contents |
|-------|----------|
| `bioconductor/bioconductor_docker:RELEASE_3_18` | R 4.3 + BiocManager |
| `bioconductor/bioconductor_docker:devel` | Latest development |

For specific Bioconductor tools, prefer biocontainers:

```
quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_0
quay.io/biocontainers/bioconductor-edger:4.0.2--r43hf17093f_0
quay.io/biocontainers/bioconductor-limma:3.58.1--r43hf17093f_0
```

## Conda environments

### Base R + tidyverse

```nextflow
conda 'conda-forge::r-base=4.3 conda-forge::r-tidyverse=2.0'
```

### Bioconductor packages

```nextflow
conda 'conda-forge::r-base=4.3 bioconda::bioconductor-deseq2=1.42'
```

### Common R packages by domain

**Data wrangling:**
```
r-tidyverse r-data.table r-readxl r-jsonlite r-xml2
```

**Statistics:**
```
r-survival r-lme4 r-mass r-car r-multcomp
```

**Visualization:**
```
r-ggplot2 r-pheatmap r-complexheatmap r-rcolorbrewer r-scales r-patchwork
```

**Bioinformatics (bioconda channel):**
```
bioconductor-deseq2 bioconductor-edger bioconductor-limma
bioconductor-genomicranges bioconductor-biostrings
bioconductor-rtracklayer bioconductor-clusterprofiler
bioconductor-org.hs.eg.db bioconductor-enrichplot
```

**Single-cell (bioconda channel):**
```
bioconductor-singlecellexperiment bioconductor-scran
bioconductor-scater r-seurat
```

**Spatial (conda-forge + bioconda):**
```
r-seurat r-spatialexperiment bioconductor-spatiallibd
```

## Matching R packages to containers

When converting, map `library()` calls:

```r
library(DESeq2)        → bioconda::bioconductor-deseq2
library(tidyverse)     → conda-forge::r-tidyverse
library(ggplot2)       → conda-forge::r-ggplot2
library(Seurat)        → conda-forge::r-seurat
library(pheatmap)      → conda-forge::r-pheatmap
library(data.table)    → conda-forge::r-data.table
library(optparse)      → conda-forge::r-optparse
library(argparse)      → conda-forge::r-argparse
```

### Conda channel priority

```
channels:
  - conda-forge    # Most R packages
  - bioconda       # Bioconductor + bioinformatics
  - defaults       # Fallback
```

## Tips

- Pin R version explicitly: `r-base=4.3` not just `r-base`
- For reproducibility, pin all package versions
- Use mulled containers for multiple biocontainers dependencies
- When in doubt, use `rocker/tidyverse` as a base and install extras
