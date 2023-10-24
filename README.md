# Snakemake workflow: smallRNA Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.9.0-brightgreen.svg)](https://snakemake.github.io/)
[![DOI](https://zenodo.org/badge/631045084.svg)](https://zenodo.org/badge/latestdoi/631045084)

Workflow for the analysis of small RNA-seq data for Yan et. al, "An endogenous
small RNA-binding protein safeguards prime editing" (in press).

The workflow is written using [Snakemake](https://snakemake.github.io/) and
[Quarto](https://quarto.org/).

Dependencies are installed using [Bioconda](https://bioconda.github.io/) where
possible.

The workflow consists of two pieces, one written in Snakemake, the other is
composed of Quarto notebooks.

## Snakemake workflow

### Setup environment and run workflow

1. Clone workflow into working directory

    ```bash
    git clone <repository> <dir>
    cd <dir>
    ```

2. Download input data

    Copy data to `data` directory

3. Edit config as needed

    ```bash
    nano config.yaml
    ```

4. Install dependencies into isolated environment

    ```bash
    conda env create -n <project> --file environment.yaml
    ```

5. Activate environment

    ```bash
    source activate <project>
    ```

6. Execute main workflow

    ```bash
    snakemake --cores 1
    ```

## Quarto notebooks

The Quarto notebooks utilize [R](https://www.r-project.org/) and are run
separately.

1. Run the workflow as above

2. Load the Rproject `pe-small-rna-seq-analysis.Rproj` in RStudio.

3. This project uses
   [`renv`](https://rstudio.github.io/renv/articles/renv.html) to keep track of
   installed packages. Install `renv` if not installed and load dependencies
   with `renv::restore()`.

4. Load one of the quarto notebooks below and notebook and run all of the cells
   or use the "Render" button in RStuido.

   * `biotype-comparison.qmd`
   * `fragment-size-distributions.qmd`
   * `alignment_statistics.qmd`
   * `coverage-plots.qmd`

5. Some of the notebooks use parameters to generate a few different versions of
   the plots. If Quarto and all of the required R packages are installed, you
   can use the `render_quarto_reports.sh` script to render all of the quarto
   notebooks.
