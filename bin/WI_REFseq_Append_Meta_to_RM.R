#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args = commandArgs(trailingOnly=TRUE)

chunk_name <- strsplit(args[2], split = "\\.fa")[[1]][1]

# load meta data file
contig_meta <- readr::read_table2(args[1], col_names = F) %>%
  na.omit()

colnames(contig_meta) <- c("CHROM","POS","END","ID","SVTYPE",
                           "SVLEN","CALLER", "REF", "ALT", 
                           "ANN", "GT_FT","SAMPLES")

contig_meta <- contig_meta %>%
  tidyr::unite(query1, CHROM, POS, END, sep = "-", remove = F) %>%
  tidyr::unite(query, query1, SVTYPE, sep = "_", remove = F) %>%
  dplyr::select(-query1)

# repeatmasker output
if(args[3] == "none"){
  empty_df <- data.frame()
  
  write.table(empty_df, file = glue::glue("{chunk_name}_RM_processed_{args[3]}.txt"), 
              quote = F, 
              col.names = F, 
              row.names = F)
  
} else {
  rm_output <- readr::read_table2(args[2], col_names = F, skip =3) 
  
  colnames(rm_output) <- c("bit_score","perc_div","perc_del","perc_ins",
                           "query","query_match_start","query_match_end",
                           "remaining_query","repeat_strand","repeat_name",
                           "repeat_class","repeat_seq_before_match",
                           "repeat_match_start","repeat_match_end",
                           "rm_id","multi_match")
  
  rm_meta <- dplyr::inner_join(contig_meta,rm_output, by = "query")
  
  ubit <- length(unique(rm_meta$bit_score))
  uclass <- length(unique(rm_meta$repeat_class))
  uname <- length(unique(rm_meta$repeat_name))
  query_length <- nrow(rm_output)
  
  id <- strsplit(as.character(runif(100)[1]), split = "\\.")[[1]][[2]]

  write.table(rm_meta, file = glue::glue("{chunk_name}_RM_processed_{ubit}_{uclass}_{uname}_{query_length}_{id}.txt"), 
              quote = F, 
              col.names = F, 
              row.names = F)
}


