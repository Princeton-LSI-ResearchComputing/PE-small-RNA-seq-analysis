library("UpSetR")
library("readr")

ncrna_intersections <- read_tsv(
  paste0(
    "data/references/rna_central/genome_coordinates/",
    "homo_sapiens.GRCh38.multiinter.tsv"
  ),
  col_types = list(
    chrom = col_character(),
    start = col_integer(),
    end = col_integer(),
    num = col_integer(),
    list = col_character(),
    .default = col_integer()
  )
)

ncrna_upset_data <- ncrna_intersections %>%
  dplyr::select(-chrom, -start, -end, -num, -list)

upset(ncrna_upset_data, main.bar.color = "black")

upset(as.data.frame(ncrna_upset_data),
  nsets = 30,
  main.bar.color = "black",
  order.by = "degree",
  decreasing = TRUE,
  group.by = "sets",
  cutoff = 5
)
