#!/usr/bin/env bash
set -euxo pipefail

quarto render biotype-comparison.qmd

quarto render differential-expression.qmd

quarto render coverage-plots.qmd \
    --to html \
    -P normalization:human_small_rna \
    --output coverage-plots-norm-human-small-rna.html

quarto render coverage-plots.qmd \
    --to html \
    -P normalization:exogenous_rna \
    --output coverage-plots-norm-exogenous-rna.html

quarto render coverage-plots.qmd \
    --to html \
    -P normalization:exogenous_rna_category \
    --output coverage-plots-norm-exogenous-rna-category.html
