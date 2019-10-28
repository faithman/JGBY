#!/usr/bin/env Rscript

library(tidyverse)
library(rrBLUP)

args <- commandArgs(trailingOnly = TRUE)

manta_repeat_matrix <- readr::read_table2(file = args[1], col_names = T)

refseq_repeat_matrix <- readr::read_table2(file = args[2], col_names = T)

combined_repeat_matrix <- dplyr::bind_rows(manta_repeat_matrix,
                                           refseq_repeat_matrix) %>%
  dplyr::arrange(CHROM, POS)

combined_repeat_matrix <- dplyr::bind_rows(manta_repeat_matrix,
                                           refseq_repeat_matrix) %>%
  dplyr::arrange(CHROM, POS) %>%
  dplyr::mutate(ALT=strsplit(ALT,",")) %>%
  tidyr::unnest(ALT) %>%
  dplyr::group_by(MARKER) %>%
  dplyr::mutate(ALT_ALLELE = 1:n()) %>%
  dplyr::ungroup() %>%
  tidyr::gather(SAMPLE, GT, -CHROM, -POS, -MARKER, -REF, -ALT, -ALT_ALLELE) %>%
  dplyr::arrange(CHROM, POS) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(make_biallelic = ifelse(grepl(ALT_ALLELE, GT) | 
                                          GT == "0/0" | GT == "0/1",
                                        "CorrectALT", "WrongALT")) %>%
  dplyr::mutate(GT = replace(GT, 
                             make_biallelic == "CorrectALT" & 
                               GT != "0/0" & 
                               GT != "0/1", "1/1")) %>%
  dplyr::filter(make_biallelic != "WrongALT") %>%
  dplyr::select(-ALT_ALLELE, -make_biallelic) %>%
  dplyr::mutate(GT_n = recode(GT, "1/1" = 1, "0/1" = 0, "0/0" = -1)) %>%
  dplyr::select(-GT) %>%
  tidyr::spread(SAMPLE, GT_n, fill = -1) 

combined_repeat_matrix_meta <- combined_repeat_matrix %>%
  tidyr::separate(MARKER, into = c("chrom", "pos", "alt", "repeat_name", "repeat_class", "temp1"), sep = "_") %>%
  tidyr::unite(repeat_class1, repeat_class, temp1, sep = "_") %>%
  dplyr::rename(repeat_class = repeat_class1) %>%
  dplyr::select(-chrom, -pos, -alt)

combined_repeat_matrix_meta$repeat_class <- gsub("_NA", "", combined_repeat_matrix_meta$repeat_class)

TE <- combined_repeat_matrix_meta %>%
  dplyr::filter(repeat_name != ")n",
                !repeat_class %in% c("Satellite", "Simple_repeat", 
                                     "Low_complexity", "Unknown")) %>%
  dplyr::distinct(repeat_class, .keep_all = T) %>%
  dplyr::pull(repeat_class)

small_repeats <- combined_repeat_matrix_meta %>%
  dplyr::filter(!(repeat_class %in% TE),
                !(repeat_class %in% c("Low_complexity", "Unknown"))) %>%
  dplyr::distinct(repeat_class, .keep_all = T) %>%
  dplyr::pull(repeat_class)

unk_lowcomp <- combined_repeat_matrix_meta %>%
  dplyr::filter(!(repeat_class %in% TE),
                !(repeat_class %in% small_repeats)) %>%
  dplyr::distinct(repeat_class, .keep_all = T) %>%
  dplyr::pull(repeat_class)

dir.create("Plots")

analyze_g_matrix <- function(df, repeats, name){
  
  geno_matrix <- df %>%
    dplyr::filter(repeat_class %in% repeats)
  
  repeat_kinship <- rrBLUP::A.mat(t(geno_matrix[,7:ncol(geno_matrix)]))
  
  repeat_pca <- prcomp(repeat_kinship)
  
  repeat_pca_df <- data.frame(SAMPLE = row.names(repeat_pca$x), repeat_pca$x)
  
  pc_plt <- ggplot(repeat_pca_df) +
    aes(x = PC1, y = PC2)+
    geom_point()+
    theme_bw(15) + 
    labs(title = glue::glue("{name} Kinship PCA"))
  
  ggsave(pc_plt, filename = glue::glue("Plots/{name}_Kinship_PCA.pdf"), 
         height = 6, width = 8)
  
  write.table(geno_matrix,
              file = glue::glue("{name}_Combined_Gmatrix.tsv"),
              sep = "\t", 
              quote = F, 
              row.names = F, 
              col.names = T)
  
  write.table(repeat_kinship,
              file = glue::glue("{name}_Kinship.tsv"),
              sep = "\t", 
              quote = F, 
              row.names = F, 
              col.names = T)
  
  write.table(repeat_pca_df,
              file = glue::glue("{name}_Kinship_PCA.tsv"),
              sep = "\t", 
              quote = F, 
              row.names = F, 
              col.names = T)
}

analyze_g_matrix(combined_repeat_matrix_meta, TE, "Transposon")
analyze_g_matrix(combined_repeat_matrix_meta, small_repeats, "SimpleRepeats")
analyze_g_matrix(combined_repeat_matrix_meta, unk_lowcomp, "Unknown")
analyze_g_matrix(combined_repeat_matrix_meta, c(TE,small_repeats,unk_lowcomp), "All_Repeats")


