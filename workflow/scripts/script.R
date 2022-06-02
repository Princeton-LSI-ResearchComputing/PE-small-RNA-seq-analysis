#!/usr/bin/env Rscript --vanilla

library(docopt, quietly=TRUE)

'Example script.

Usage:
  script.R <input> [options] 
  script.R -h | --help
  script.R --version

Options:
  -h --help                    Show this screen.
  --version                    Show version.
  --output=<output>            Output filename [default: output.csv].
  --threads=<threads>          Number off threads [default: 1].
  --threshold=<threshold>      Threshold parameter [default: 25].
  --conifgparam=<configparam>  Configuration parameter [default: 1].

' -> doc

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  args_to_parse <- c(
    snakemake@input[["input"]],
    "--output", snakemake@output[["output_filename"]],
    "--threads", snakemake@threads,
    "--threshold", snakemake@params[["threshold"]],
    "--configparam", snakemake@config[["param"]],
  )
} else if (interactive()) {
    # Arguments supplied inline (for debug/testing when running interactively)
    print("Running interactively...")
  input <- "results/input.csv"
  output <- "results/my_output.csv"
  args_to_parse <- c(input, 
                     "--output", output)
} else {
  # Arguments from command line
  args_to_parse <- commandArgs(trailingOnly=TRUE) 
}

arguments <- docopt(doc, args = args_to_parse, version = 'Script 0.1')
print(arguments)


################
# Do Something #
################

do_something <- function(data_path, out_path, threads, threshold, myparam) {
  # R code
}


########
# MAIN #
########
do_something(arguments$input, arguments$output, arguments$threads, arguments$threshold, arguments$configparam)

