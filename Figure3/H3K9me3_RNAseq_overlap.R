# load libraries

library(tibble)
library(dplyr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)
library(DESeq2)

# load in peaks
sig_peaks_p05 <- import("path_to_H3K9me3_peaks/k9_tiling90_k9_seacrP05_reducedSet_DOWNonly.bed")
sig_peaks_p01 <- import("path_to_H3K9me3_peaks/k9_tiling90_k9_seacrP01_reducedSet_DOWNonly.bed")
sig_peaks_p001 <- import("path_to_H3K9me3_peaks/k9_tiling90_k9_seacrP001_reducedSet_DOWNonly.bed")

k9_peak_list <- list(`p < 0.05` = sig_peaks_p05, 
                     `p < 0.01` = sig_peaks_p01, 
                     `p < 0.001` = sig_peaks_p001)


# load in genes tht go up in rna seq
rna_up_p05FC2 <- read.table("path_to_DE_genes/G7_B10_overlap_p05_FC2_UPgenes_ID.txt") %>%
  pull(V1)

# annotate K9me3 peaks with genomic features and nearest gene
cs_all <- lapply(k9_peak_list, 
                 annotatePeak,
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                 annoDb = "org.Mm.eg.db",
                 tssRegion=c(-3000, 3000))

# gene list for all sig peaks 
cs_all_symbol <- lapply(cs_all, function(x){
  df_temp <- as.data.frame(x)
  df_temp$SYMBOL %>%
    unique
})

cs_all_id <- lapply(cs_all, function(x){
  df_temp <- as.data.frame(x)
  df_temp$geneId %>%
    unique
})

# make venn diagram
# the venn diagram for K9me3 peaks at p05 threshold is used in Figure 3D
dir.create("cs_venn_diagrams")
for(i in seq_along(cs_all_id)){
  p_val <- gsub(" < 0.", "", names(cs_all_id)[i])
  venn_list <- list(k9_down = cs_all_id[[i]], rna_up_p05FC2)
  venn.diagram(venn_list, filename = paste0("cs_venn_diagrams/cs_k9me3_down_", p_val, "_RNAseqFC2_bothUP.png"), 
               disable.logging = T, 
               main = paste0("k9me3 down (",p_val, ") over RNAseq UP (P05 + FC2)")
  )
}
