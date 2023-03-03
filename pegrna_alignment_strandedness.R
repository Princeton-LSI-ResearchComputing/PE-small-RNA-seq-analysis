source("pegrna_plots.R")

pegrna_alignment_strandedness <- function(sequence_name,
                                          mix = NA,
                                          start_offset = NA) {
  plot_data <- get_pegrna_plot_data(
    sequence_name = sequence_name,
    normalization_factor = "grch38_mapped_reads",
    mix = mix
  )

  counts_by_strand <-
    tibble(
      sequence_name = character(),
      sample_unit = character(),
      type = character(),
      strand = character(),
      count = numeric()
    )
  for (s in names(plot_data)) {
    if (is.na(start_offset)) {
      primary_granges <- plot_data[[s]][["concordant_granges"]]
      primary_type <- "concordant"
      secondary_granges <- plot_data[[s]][["discordant_granges"]]
      secondary_type <- "discordant"
    } else {
      concordant_granges <- plot_data[[s]][["concordant_granges"]]
      start_site <- GRanges(
        seqnames = rna_species,
        ranges = IRanges(
          start = plot_range$start,
          end = plot_range$start + start_offset
        )
      )
      primary_granges <- subsetByOverlaps(concordant_granges, start_site)
      primary_type <- "concordant_full"
      secondary_granges <-
        concordant_granges[concordant_granges %outside% start_site]
      secondary_type <- "concordant_partial"
    }
    counts_by_strand <- rows_append(
      counts_by_strand,
      tibble(
        sequence_name = sequence_name,
        sample_unit = s,
        type = primary_type,
        strand = "+",
        count = sum(strand(primary_granges) == "+")
      )
    )
    counts_by_strand <- rows_append(
      counts_by_strand,
      tibble(
        sequence_name = sequence_name,
        sample_unit = s,
        type = primary_type,
        strand = "-",
        count = sum(strand(primary_granges) == "-")
      )
    )
    counts_by_strand <- rows_append(
      counts_by_strand,
      tibble(
        sequence_name = sequence_name,
        sample_unit = s,
        type = secondary_type,
        strand = "+",
        count = sum(strand(primary_granges) == "+")
      )
    )
    counts_by_strand <- rows_append(
      counts_by_strand,
      tibble(
        sequence_name = sequence_name,
        sample_unit = s,
        type = secondary_type,
        strand = "-",
        count = sum(strand(primary_granges) == "-")
      )
    )
  }
  return(counts_by_strand)
}
