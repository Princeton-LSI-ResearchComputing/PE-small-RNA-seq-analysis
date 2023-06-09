rna_species_plot_range <- function(sequence_name) {
  plot_range <- rna_mixes %>%
    filter(.data$rna_species == sequence_name) %>%
    select(start, end) %>%
    unique()
  return(plot_range)
}

### GRanges from BAM
granges_from_bam <- function(sample_unit, sequence_name, is_proper_pair) {
  plot_range <- rna_species_plot_range(sequence_name)
  which <- GRanges(sprintf(
    "%s:%i-%i", sequence_name, plot_range$start, plot_range$end
  ))

  param <- ScanBamParam(
    flag = scanBamFlag(isProperPair = is_proper_pair),
    mapqFilter = 1,
    which = which
  )

  if (is_proper_pair) {
    aligned_fragments <- granges(
      readGAlignmentPairs(
        sprintf("%s/%s.bam", bam_dir, sample_unit),
        param = param
      ),
      on.discordant.seqnames = "drop"
    )
  } else {
    aligned_fragments <- granges(
      readGAlignments(
        sprintf("%s/%s.bam", bam_dir, sample_unit),
        param = param
      )
    )
  }
  return(aligned_fragments)
}


###
get_pegrna_plot_data <- function(sequence_name,
                                 normalization_factor,
                                 mix = NA) {
  # Get list of samples containing this pegrna species
  samples <- sample_units %>%
    full_join(rna_mixes, by = "exogenous_rna") %>%
    arrange(.data$exogenous_rna, .data$day, .data$cell_line) %>%
    filter(.data$rna_species == sequence_name)
  # Limit by mix if specified (needed for the control sample in both mixes)
  if (!is.na(mix)) {
    samples <- samples %>% filter(.data$exogenous_rna == mix)
  }
  samples <- samples %>% pull(.data$sample_unit)

  # Setup list to hold data
  plot_data <- list()

  # Iterate through samples, aggregate data
  for (s in samples) {
    # Look up normalization read count
    if (normalization_factor == "exogenous_rna_mapped_reads") {
      norm_seqname <- sequence_name
      read_norm_factor <- 1
    } else if (normalization_factor == "grch38_mapped_reads") {
      norm_seqname <- "grch38_mapped_reads"
      read_norm_factor <- 1000000
    }
    normalization_read_count <- idxstats %>%
      filter(sample == s) %>%
      filter(sequence_name == norm_seqname) %>%
      pull(.data$mapped_reads)
    normalization_read_count <- normalization_read_count / read_norm_factor

    # Look up cell line
    cell_line <- sample_units %>%
      filter(.data$sample_unit == s) %>%
      pull(cell_line)

    # Look up day
    day <- sample_units %>%
      filter(.data$sample_unit == s) %>%
      pull(day)

    # Discordant
    discordant <- granges_from_bam(s, sequence_name, is_proper_pair = FALSE)
    cov_discordant <- coverage(discordant)
    cov_discordant_norm <- cov_discordant / normalization_read_count

    # Concordant
    concordant <- granges_from_bam(s, sequence_name, is_proper_pair = TRUE)
    cov_concordant <- coverage(concordant)
    cov_concordant_norm <- cov_concordant / normalization_read_count

    # Gather plot data for this sample into list
    plot_data[[s]] <- list(
      cell_line = cell_line,
      day = day,
      discordant_granges = discordant,
      discordant = cov_discordant_norm,
      concordant_granges = concordant,
      concordant = cov_concordant_norm,
      normalization_read_count = normalization_read_count
    )
  }

  return(plot_data)
}

###
pegrna_plots <- function(sequence_name,
                         normalization_factor,
                         mix = NA,
                         ylim = NA,
                         ylab,
                         vlines = NA,
                         start_offset = NA) {
  plot_data <- get_pegrna_plot_data(
    sequence_name = sequence_name,
    normalization_factor = normalization_factor,
    mix = mix
  )
  plot_range <- rna_species_plot_range(sequence_name)

  # Find largest value for rna_species
  if (is.na(ylim)) {
    m <- 0
    for (s in names(plot_data)) {
      m <- max(
        m,
        max(plot_data[[s]][["concordant"]][[sequence_name]]),
        max(plot_data[[s]][["discordant"]][[sequence_name]])
      )
    }
    ylim <- c(0, m)
  }

  last_day <- 0

  par(mfrow = c(2, 1))

  for (s in names(plot_data)) {
    if (is.na(start_offset)) {
      series_data_concordant <- plot_data[[s]][["concordant"]][[sequence_name]]
      series_data_discordant <- plot_data[[s]][["discordant"]][[sequence_name]]
      concordant_legend_text <- "Concordant"
      discordant_legend_text <- "Discordant"
    } else {
      concordant_granges <- plot_data[[s]][["concordant_granges"]]
      normalization_read_count <- plot_data[[s]][["normalization_read_count"]]
      start_site <- GRanges(
        seqnames = rna_species,
        ranges = IRanges(
          start = plot_range$start,
          end = plot_range$start + start_offset
        )
      )
      full_length_granges <- subsetByOverlaps(concordant_granges, start_site)
      partial_granges <-
        concordant_granges[concordant_granges %outside% start_site]
      series_data_concordant <-
        coverage(full_length_granges)[[sequence_name]] /
          normalization_read_count
      series_data_discordant <-
        coverage(partial_granges)[[sequence_name]] / normalization_read_count
      concordant_legend_text <- sprintf(
        "Start within %s bp of 5'",
        start_offset
      )
      discordant_legend_text <- sprintf(
        "Start after %s bp of 5'",
        start_offset
      )
    }


    day <- plot_data[[s]][["day"]]
    cell_line <- plot_data[[s]][["cell_line"]]
    concordant_color <- concordant_cell_line_colors[[cell_line]]
    discordant_color <- discordant_cell_line_colors[[cell_line]]

    if (day != last_day) {
      plot(series_data_concordant,
        type = "l",
        col = concordant_color,
        xlim = c(plot_range$start, plot_range$end),
        ylim = ylim,
        main = sprintf("Day %s", day),
        xlab = sprintf("%s position", sequence_name),
        ylab = ylab
      )
      lines(series_data_discordant,
        type = "l",
        col = discordant_cell_line_colors[[plot_data[[s]][["cell_line"]]]]
      )
    } else {
      lines(series_data_concordant,
        type = "l",
        col = concordant_color
      )
      lines(series_data_discordant,
        type = "l",
        col = discordant_color
      )
    }
    legend("top",
      legend = c(
        sprintf("Parental:%s", concordant_legend_text),
        sprintf("P1E10:%s", concordant_legend_text),
        sprintf("Parental:%s", discordant_legend_text),
        sprintf("P1E10:%s", discordant_legend_text)
      ),
      col = unlist(c(
        concordant_cell_line_colors,
        discordant_cell_line_colors
      )),
      lty = 1,
      cex = 0.75
    )
    for (vline in vlines) {
      abline(v = vline, lty = 2, lwd = 0.25)
      text(
        x = vline, y = 0, as.character(vline),
        pos = 2, offset = 0.25, cex = 0.75,
      )
    }
    last_day <- day
  }
}
