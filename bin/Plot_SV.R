library(tidyverse)

# 1 - SV BED FILE
# 2 - SAMPLE NAME

args <- commandArgs(TRUE)

sv_bed <- data.table::fread(args[1],
                            col.names = c("CHROM", "START", "END", "SUPPORT", "SVTYPE", "STRAND", "SVTYPE_CALLER", "SVpos_CALLER", "CALLER", "GT", "SNPEFF_TYPE", "SNPEFF_PRED", "SNPEFF_EFF", "TRANSCRIPT", "WBGeneID")) 

sv_bed_pr <- sv_bed %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  tidyr::separate(CALLER, into = c("SAMPLE", "CALLER"), sep = "_") %>%
  dplyr::group_by(CHROM, START, END) %>%
  dplyr::mutate(SVTYPE = ifelse(CALLER == "manta" & SVTYPE != "BND", SVTYPE_CALLER, SVTYPE)) %>%
  dplyr::mutate(SVTYPE = ifelse(grepl(",", SVTYPE_CALLER), SVTYPE_CALLER, SVTYPE)) %>%
  dplyr::mutate(SVTYPE_CLEAN = ifelse(grepl(",", SVTYPE), "COMPLEX", SVTYPE)) %>%
  dplyr::mutate(TRANSCRIPT = ifelse(TRANSCRIPT == "", "Intergenic", TRANSCRIPT),
                WBGeneID = ifelse(WBGeneID == "", "Intergenic", WBGeneID)) %>%
  transform(TRANSCRIPT = strsplit(TRANSCRIPT, "\\&")) %>%
  tidyr::unnest(TRANSCRIPT) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SIZE = END-START,
                HIGH_EFF = ifelse(SNPEFF_EFF == "HIGH", "HIGH", "LOW-MOD")) %>%
  dplyr::filter(SVTYPE != "BND", !grepl("TRA", SVTYPE)) %>%
  transform(WBGeneID = strsplit(WBGeneID, "\\&")) %>%
  tidyr::unnest(WBGeneID) %>%
  dplyr::distinct(CHROM, START, END, WBGeneID, .keep_all = T) %>%
  dplyr::filter(!(CHROM == "MtDNA" & (SVTYPE_CLEAN == "DUP" | SVTYPE_CLEAN == "INV" | SVTYPE_CLEAN == "COMPLEX"))) 

write.table(x = sv_bed_pr,
            file = glue::glue("{args[2]}_processed_SV.bed"), quote = F, col.names = T, row.names = F, sep = "\t")

sv_bed_pr %>%
  dplyr::filter(SIZE < 1e5) %>%
  ggplot() +
  aes(x = START/1e6, xend = END/1e6, y=SIZE, yend=SIZE, color = SVTYPE_CLEAN) +
  geom_segment(size =1, arrow = arrow(length = unit(0.1,"cm")))+
  facet_grid(CHROM~HIGH_EFF, scales = "free")+
  theme_bw(18) +
  theme(strip.background = element_blank(),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray70"),
        panel.grid = element_blank())+
  scale_color_manual(values =c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(x = "Genomic Position", color = "SVTYPE") 

ggsave(glue::glue("{args[2]}_SVplot.pdf"), height = 10, width = 12)
