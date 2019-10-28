#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args = commandArgs(trailingOnly=TRUE)

chunk_name <- strsplit(args[2], split = "\\.")[[1]][1]

# load meta data file
contig_meta <- readr::read_table2(args[1], col_names = F) %>%
  na.omit()

contig_meta$X11 <- gsub(">", "", contig_meta$X11)

colnames(contig_meta) <- c("CHROM","POS","REF","ALT","SVTYPE",
                           "CONTIG","group_GT","SAMPLES","Contig_ID","Contig_range","query")

# repeatmasker output
if(args[3] == "none"){
  empty_df <- data.frame()
  
  write.table(empty_df, file = glue::glue("{chunk_name}_RM_processed.txt"), 
              quote = F, 
              col.names = F, 
              row.names = F)
  
} else {
  rm_output <- readr::read_table2(args[2], col_names = F, skip =3) 
  
  colnames(rm_output) <- c("bit_score","perc_div","perc_del","perc_ins",
                           "query","query_match_start","query_match_end","remaining_query",
                           "repeat_strand","repeat_name","repeat_class","repeat_seq_before_match",
                           "repeat_match_start","repeat_match_end","rm_id","multi_match")
  
  
  rm_meta <- dplyr::inner_join(contig_meta,rm_output, by = "query")
  
  write.table(rm_meta, file = glue::glue("{chunk_name}_RM_processed.txt"), 
              quote = F, 
              col.names = F, 
              row.names = F)
}