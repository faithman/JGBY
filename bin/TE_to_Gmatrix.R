#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

processed_repeats <- readr::read_table2(file = args[1], col_names = F)

colnames(processed_repeats) <- c("CHROM","POS","REF","ALT","SVTYPE",
                                 "CONTIG","group_GT","SAMPLES","Contig_ID","Contig_range","query",
                                 "bit_score","perc_div","perc_del","perc_ins",
                                 "query_match_start","query_match_end","remaining_query",
                                 "repeat_strand","repeat_name","repeat_class","repeat_seq_before_match",
                                 "repeat_match_start","repeat_match_end","rm_id","multi_match")


processed_repeats_long <- processed_repeats %>%
  dplyr::mutate(SAMPLES=strsplit(SAMPLES,",")) %>%
  tidyr::unnest(SAMPLES)%>%
  tidyr::separate(SAMPLES,c("SAMPLE","FT"),":")%>%
  dplyr::mutate(FT = replace(FT, FT == ".", "PASS"))

TE <- processed_repeats_long %>%
  dplyr::filter(repeat_name != ")n",
                !repeat_class %in% c("Satellite", "Simple_repeat", "Low_complexity", "Unknown")) %>%
  dplyr::distinct(repeat_class, .keep_all = T) %>%
  dplyr::pull(repeat_class)

repeats <- processed_repeats_long %>%
  dplyr::filter(!(repeat_class %in% TE),
                !(repeat_class %in% c("Low_complexity", "Unknown"))) %>%
  dplyr::distinct(repeat_class, .keep_all = T) %>%
  dplyr::pull(repeat_class)

unk_lowcomp <- processed_repeats_long %>%
  dplyr::filter(!(repeat_class %in% TE),
                !(repeat_class %in% repeats)) %>%
  dplyr::distinct(repeat_class, .keep_all = T) %>%
  dplyr::pull(repeat_class)

df_to_Gmatrix <- function(df, repeats, name) {
  
  geno_matrix <- df %>%
    dplyr::group_by(CHROM, POS, ALT) %>%
    dplyr::filter(bit_score == max(bit_score)) %>%
    dplyr::ungroup()%>%
    dplyr::filter(repeat_class %in% repeats) %>%
    tidyr::unite(MARKER, CHROM, POS, ALT, repeat_name, repeat_class, sep = "_", remove = F) %>%
    dplyr::filter(FT == "PASS") %>%
    dplyr::select(CHROM, POS, MARKER, REF, ALT, SAMPLE, GT = group_GT) %>%
    dplyr::distinct() %>%
    tidyr::spread(SAMPLE, GT) %>%
	dplyr::arrange(CHROM, POS)
  
  write.table(geno_matrix,
              file = glue::glue("{name}_Gmatrix.tsv"),
              sep = "\t", 
              quote = F, 
              row.names = F, 
              col.names = T)
}

df_to_Gmatrix(processed_repeats_long, TE, "Transposon")
df_to_Gmatrix(processed_repeats_long, repeats, "SimpleRepeats")
df_to_Gmatrix(processed_repeats_long, unk_lowcomp, "Unknown")
df_to_Gmatrix(processed_repeats_long, c(TE,repeats,unk_lowcomp), "All_Repeats")


