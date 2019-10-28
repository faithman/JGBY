#!/usr/bin/env Rscript

library(tidyverse)
library(rrBLUP)

args <- commandArgs(trailingOnly = TRUE)

manta_repeat_matrix <- readr::read_table2(file = args[1], col_names = T)

refseq_repeat_matrix <- readr::read_table2(file = args[2], col_names = T)

combined_repeat_matrix <- dplyr::bind_rows(manta_repeat_matrix,
                                           refseq_repeat_matrix) %>%
  dplyr::filter(REF != "N", ALT != "N") %>%
  dplyr::arrange(CHROM, POS) %>%
  dplyr::mutate(ALT=strsplit(ALT,",")) %>%
  tidyr::unnest(ALT) %>%
  dplyr::group_by(MARKER) %>%
  dplyr::mutate(ALT_ALLELE = 1:n()) %>%
  dplyr::ungroup() %>%
  tidyr::gather(SAMPLE, GT, -CHROM, -POS, -MARKER, -REF, -ALT, -ALT_ALLELE) %>%
  dplyr::arrange(CHROM, POS) 

# remove sites with only HET calls
has_hom_alt <- combined_repeat_matrix%>%
  dplyr::group_by(CHROM,POS,MARKER, GT)%>%
  dplyr::tally() %>%
  dplyr::mutate(no_hom_alt = ifelse(!"1/1" %in% unique(GT) , "NO_HOM_ALT", "HOM_ALT"))%>%
  dplyr::filter(no_hom_alt != "NO_HOM_ALT")

combined_repeats_hom <- combined_repeat_matrix %>%
  dplyr::filter(MARKER %in% has_hom_alt$MARKER) 

# frequency cutoff
af_cutoff <- combined_repeats_hom %>%
  dplyr::group_by(CHROM,POS,MARKER, GT)%>%
  dplyr::tally() %>%
  dplyr::mutate(AF = n/sum(n)) %>%
  dplyr::filter(GT == "0/0") %>%
  dplyr::mutate(perc5_cut = ifelse(AF < 0.05 | AF > 0.95, "RARE", "COMMON")) %>%
  dplyr::filter(perc5_cut != "RARE")

combined_repeats_hom_af <- combined_repeats_hom %>%
  dplyr::filter(MARKER %in% af_cutoff$MARKER)

# pull REF genotypes
combined_repeats_ref <- combined_repeats_hom %>%
  dplyr::filter(MARKER %in% af_cutoff$MARKER) %>%
  dplyr::filter(GT == "0/0")%>%
  dplyr::mutate(GT1 = GT)

# convert multiallelic alts
combined_repeats_alt <- combined_repeats_hom_af %>%
  dplyr::filter(GT != "0/0") %>%
  dplyr::rowwise() %>%
  dplyr::filter(grepl(glue::glue("/{ALT_ALLELE}$"), GT)) %>%
  dplyr::mutate(GT1 = ifelse(grepl("0/", GT), glue::glue("0/1"), "1/1"))
  
combined_repeats_fixed <- dplyr::bind_rows(combined_repeats_ref, combined_repeats_alt) %>%
  dplyr::arrange(CHROM,POS,SAMPLE) %>%
  dplyr::mutate(GT_n = recode(GT1, "1/1" = 1, "0/1" = 0, "0/0" = -1)) %>%
  dplyr::select(-GT, -GT1) %>%
  tidyr::spread(SAMPLE, GT_n, fill = -1) 

combined_repeat_matrix_meta <- combined_repeats_fixed %>%
  tidyr::separate(MARKER, into = c("chrom", "pos", "alt", "repeat_name", "repeat_class", "temp1"), sep = "_") %>%
  tidyr::unite(repeat_class1, repeat_class, temp1, sep = "_") %>%
  dplyr::rename(repeat_class = repeat_class1) %>%
  dplyr::select(-chrom, -pos, -alt) %>%
  dplyr::select(-ALT_ALLELE)

combined_repeat_matrix_meta$repeat_class <- gsub("_NA", "", combined_repeat_matrix_meta$repeat_class)

repeat_kinship <- rrBLUP::A.mat(t(combined_repeat_matrix_meta[,7:ncol(combined_repeat_matrix_meta)]))

repeat_pca <- prcomp(repeat_kinship)

repeat_pca_df <- data.frame(SAMPLE = row.names(repeat_pca$x), repeat_pca$x)

write.table(repeat_kinship,
            file = glue::glue("{args[3]}_Kinship.tsv"),
            sep = "\t", 
            quote = F, 
            row.names = F, 
            col.names = T)

write.table(repeat_pca_df,
            file = glue::glue("{args[3]}_Kinship_PCA.tsv"),
            sep = "\t", 
            quote = F, 
            row.names = F, 
            col.names = T)

