# Snakemake workflow: smallRNA Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.6.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/troycomi/adamson_smRNA.svg?branch=master)](https://travis-ci.org/snakemake-workflows/adamson_smRNA)

Workflow for analysis of smRNA

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description, covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules` and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

The workflow is written using [Snakemake](https://snakemake.readthedocs.io/).

Dependencies are installed using [Bioconda](https://bioconda.github.io/) where possible.

## Setup environment and run workflow

1.  Clone workflow into working directory

    ```
    git clone <repo> <dir>
    cd <dir>
    ```

2.  Download input data

    Copy data from [URL]() to `data` directory

3.  Edit config as needed

    ```
    nano config.yaml
    ```

4.  Install dependencies into isolated environment

    ```
    conda env create -n <project> --file environment.yaml
    ```

5.  Activate environment

    ```
    source activate <project>
    ```

6.  Execute workflow

    ```
    snakemake -n
    ```

7. Investigate results

    After successful execution, you can create a self-contained interactive HTML report with all results via:

        snakemake --report report.html

    This report can, e.g., be forwarded to your collaborators.


## Running workflow on `gen-comp1`

```
snakemake --cluster-config cluster_config.cetus.yaml \
          --drmaa " --cpus-per-task={cluster.n} --mem={cluster.memory} --qos={cluster.qos}" \
          --use-conda -w 60 -rp -j 1000
```

## Generating the R report

Currently, the R report has to be run separately.

1. Run the workflow as above
2. Load the Rproject `adamson-smallRNA.Rproj` in RStudio
3. Load `exogenous-rna-profiles.Rmd`, run all cells, and save to create the `exogenous_targets_paired_read_summary.tsv` as well as the report `exogenous-rna-profiles.nb.html`

### R packages

1. If needed, install the `renv` package
    ```r
    install.packages("renv")
    ```
2. Load required packages
    ```r
    renv::activate()
    ```

## Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/troycomi/adamson_smRNA) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.

## Testing

Tests cases are in the subfolder `.test`. They should be executed via continuous integration with Travis CI.
