#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args = commandArgs(trailingOnly=TRUE)

# load file
contigs <- readr::read_table2(args[1], col_names = F) 

# make unique identifier for contigs
pr_contigs <- contigs %>%
  dplyr::mutate(contig_number = 1:n()) %>%
  dplyr::group_by(X1, X2, X3, X4) %>%
  dplyr::mutate(contig_range = paste(min(contig_number), 
                                     max(contig_number), 
                                     sep = "-"))%>%
  dplyr::ungroup() %>%
  dplyr::mutate(contig_fasta = glue::glue(">CONTIGS={contig_range}={X1}_{X2}\n{X6}"))

# save file

write.table(pr_contigs, file = "Processed_Manta_Contigs.tsv", 
            quote = F, 
            col.names = F, 
            row.names = F)

distinct_contigs <- dplyr::distinct(pr_contigs, X1, X2, X3, X4, .keep_all=T)

for(chrom in unique(distinct_contigs$X1)){
  
  contig_output <- dplyr::filter(distinct_contigs, X1 == chrom) %>%
    dplyr::mutate(bin = ntile(X2, 20))
  
  for(chunk in unique(contig_output$bin)){
    
    chunked_contig_output <- dplyr::filter(contig_output, bin == chunk)
    
    writeLines(
      c(
        chunked_contig_output$contig_fasta
      ),
      con = glue::glue("Chrom{chrom}_Chunk-{chunk}-Contig.fa"), 
      sep = "\n"
    )
  }
}
