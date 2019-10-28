library(tidyverse)

# 1 - SV BED FILE
# args <- "transfer/combined.bed"
args <- commandArgs(TRUE)

sv_bed <- data.table::fread(args[1]) %>%
  dplyr::filter(CHROM != "CHROM") %>%
  dplyr::mutate(START = as.numeric(START),
                END = as.numeric(END)) %>%
  dplyr::arrange(CHROM, START, END) %>%
  dplyr::select(CHROM, START, END, SVTYPE_CLEAN, SAMPLE) %>%
  dplyr::mutate(n_samples = length(unique(SAMPLE))) %>%
  dplyr::group_by(CHROM, START, END) %>%
  dplyr::distinct() %>%
  dplyr::mutate(GT = 1,
                AF = n()/n_samples) %>%
  tidyr::spread(SAMPLE, GT)

sv_bed[is.na(sv_bed)] <- -1

write.table(x = sv_bed,
            file = glue::glue("Population_SV_GenotypeMatrix.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")

