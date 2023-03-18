#!/usr/bin/env bash
set -euxo pipefail

quarto render biotype-comparison.qmd

quarto render differential-expression.qmd

quarto render coverage-plots.qmd \
    --to html \
    -P plot_other:FALSE \
    --output coverage-plots-without-other.html

quarto render coverage-plots.qmd \
    --to html \
    -P plot_other:TRUE \
    --output coverage-plots-with-other.html
