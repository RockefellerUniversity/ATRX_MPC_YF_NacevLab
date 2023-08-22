# load libraries

library(profileplyr)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rGREAT)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

# get up and down genes from RNAseq analysis
up_geneID_p05FC2 <- read.table("/Users/douglasbarrows/Desktop/BRC/collaborations/Allis/Yan/20220124_10T_AtrxKO_B10_G7_RNAseq/doug_analysis/DE_gene_lists/G7_B10_overlap_p05_FC2_UPgenes_ID.txt") %>%
  pull(V1) 

down_geneID_p05FC2 <- read.table("/Users/douglasbarrows/Desktop/BRC/collaborations/Allis/Yan/20220124_10T_AtrxKO_B10_G7_RNAseq/doug_analysis/DE_gene_lists/G7_B10_overlap_p05_FC2_DOWNgenes_ID.txt") %>%
  pull(V1) 

# make bed files for TSS of up and down RNAseq genes
mm10_genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
mm10_tss <- promoters(mm10_genes, upstream = 1, downstream = 1)
mm10_tss$name <- mm10_tss$gene_id

mm10_tss_up <- mm10_tss[mm10_tss$gene_id %in% up_geneID_p05FC2]
rtracklayer::export(mm10_tss_up, "mm10_TSS_up.bed")

mm10_tss_down <- mm10_tss[mm10_tss$gene_id %in% down_geneID_p05FC2]
rtracklayer::export(mm10_tss_down, "mm10_TSS_down.bed")

################################################################################
################################################################################
# quantify H3K4me3 signal over these TSS and make heatmap
################################################################################
################################################################################



bigwigs <- c("path_to_bigwigs/Sorted_WT_IgG_R1Normalised.bw",
             "path_to_bigwigs/Sorted_WT_IgG_R2Normalised.bw",
             "path_to_bigwigs/Sorted_WT_IgG_R4Normalised.bw",
             "path_to_bigwigs/Sorted_B10_IgG_R1Normalised.bw", 
             "path_to_bigwigs/Sorted_B10_IgG_R2Normalised.bw",
             "path_to_bigwigs/Sorted_B10_IgG_R4Normalised.bw",
             "path_to_bigwigs/Sorted_G7_IgG_R1Normalised.bw",
             "path_to_bigwigs/Sorted_G7_IgG_R2Normalised.bw",
             "path_to_bigwigs/Sorted_G7_IgG_R4Normalised.bw",
             "path_to_bigwigs/Sorted_WT_k4me3_R1Normalised.bw",
             "path_to_bigwigs/Sorted_WT_k4me3_R2Normalised.bw",
             "path_to_bigwigs/Sorted_WT_k4me3_R4Normalised.bw",
             "path_to_bigwigs/Sorted_B10_k4me3_R2Normalised.bw", 
             "path_to_bigwigs/Sorted_B10_k4me3_R3Normalised.bw", 
             "path_to_bigwigs/Sorted_B10_k4me3_R4Normalised.bw",
             "path_to_bigwigs/Sorted_G7_k4me3_R1Normalised.bw",
             "path_to_bigwigs/Sorted_G7_k4me3_R2Normalised.bw",
             "path_to_bigwigs/Sorted_G7_k4me3_R4Normalised.bw")


ranges_h33 <- c("mm10_TSS_up.bed",
                "mm10_TSS_down.bed")

cp <- BamBigwig_to_chipProfile(signalFiles = bigwigs,
                               testRanges = ranges_h33,
                               format = "bigwig", 
                               style = "point", 
                               distanceUp = 3000,
                               distanceDown = 3000)

pp <- as_profileplyr(cp)

# clean up profileplyr object
rownames(sampleData(pp)) <- gsub("Normalised.bw", "",sampleData(pp)$sample_labels)
rownames(sampleData(pp)) <- gsub("Sorted_", "",rownames(sampleData(pp)))
sampleData(pp)$chip <- gsub("WT_|G7_|B10_", "", rownames(sampleData(pp)))
sampleData(pp)$chip <- gsub("_R1|_R2|_R3", "", sampleData(pp)$chip)

mcols(pp)$peak_set <- ifelse(mcols(pp)$sgGroup == "mm10_TSS_down.bed", 
                             "DE_down_peaks", 
                             "DE_up_peaks")

pp <- groupBy(pp, group = "peak_set")

### get mean signal for replicates

antibodies <- c("IgG", "k4me3")
mean_pp_all <- list()
for (j in seq_along(antibodies)){
  pp_temp <- pp[,,grepl(antibodies[j], sampleData(pp)$chip)]
  
  groups <- paste(c("WT", "B10", "G7"), antibodies[j], sep = "_")
  mean_matrices <- lapply(groups, function(x){
    mats <- assays(pp_temp)[grepl(x, names(assays(pp_temp)))]
    mats_cbind <- do.call(cbind, mats)
    mats_array <- array(mats_cbind, dim=c(dim(mats[[1]]), length(mats)))
    
    apply(mats_array, c(1, 2), mean, na.rm = TRUE)
    
  })
  
  names(mean_matrices) <- groups
  
  tempSE <- SummarizedExperiment(assays= mean_matrices,
                                 rowRanges=rowRanges(pp_temp))
  sampleData <- sampleData(pp_temp)[c(1,4,7), ]
  rownames(sampleData) <- gsub("_R1|_R2" , "", rownames(sampleData))
  mean_pp_all[[j]] <- new("profileplyr", tempSE,
                          params=params(pp_temp),
                          sampleData=sampleData)
  
}
names(mean_pp_all) <- antibodies

mean_pp_k4 <- mean_pp_all[["k4me3"]]
mean_pp_k4_sum <- profileplyr::summarize(mean_pp_k4, fun = rowMeans, output = "matrix")
mean_pp_k4 <- mean_pp_k4[!rowSums(mean_pp_k4_sum) == 0,,]

mean_pp_k4_up <- mean_pp_k4[mcols(mean_pp_k4)$sgGroup == "mm10_TSS_up.bed"]

up_peaks_k4 <- rtracklayer::import("path_to_DE_peaks/k4me3_stringent0.01_wRep4_B10andG7_overlap_sigPeaks_up.bed")
mcols(mean_pp_k4_up)$over_sig <- ifelse(mean_pp_k4_up %over% up_peaks_k4, "sig_up", "no_change")
mcols(mean_pp_k4_up)$over_sig <- ordered(mcols(mean_pp_k4_up)$over_sig, levels = c("sig_up", "no_change"))
mean_pp_k4_up <- groupBy(mean_pp_k4_up, group = "over_sig")

# this heatmap is Figure 4A
pdf("K4me3_over_DEproms_up_meanSignal_groupByupPeaks.pdf", height = 8, width =5)
generateEnrichedHeatmap(mean_pp_k4_up,
                        ylim = "common_max",
                        matrices_pos_line = FALSE,
                        include_group_annotation = T, use_raster = T 
)
dev.off()


################################################################################
################################################################################
# Gene ontology for K4me3 peaks over promoters
################################################################################
################################################################################


mean_pp_k4_up_geneID <- mcols(mean_pp_k4_up)$name[mcols(mean_pp_k4_up)$over_sig == "sig_up"]

# get universe from differential analysis
universe_G7 <- read.table("/Users/douglasbarrows/Desktop/BRC/collaborations/Allis/Yan/20220124_10T_AtrxKO_B10_G7_RNAseq/salmon/Group_Atrx_sg5_G7_minus_Atrx_sg6_B10GOIonlyDEG.xls") %>%
  pull(V1) 
universe_B10 <- read.table("/Users/douglasbarrows/Desktop/BRC/collaborations/Allis/Yan/20220124_10T_AtrxKO_B10_G7_RNAseq/salmon/Group_Atrx_sg6_B10_minus_WTGOIonlyDEG.xls") %>%
  pull(V1) 

universe_both <- c(universe_G7, universe_B10) %>% unique

# GOBP
up_peaks_GOBP <- enrichGO(mean_pp_k4_up_geneID,
                          ont = "BP", 
                          universe = universe_both,
                          OrgDb = "org.Mm.eg.db", 
                          pvalueCutoff = 0.5, 
                          qvalueCutoff = 0.5)

up_peaks_GOBP_df <- data.frame(up_peaks_GOBP)

# this table is Supplementary Table 18 and is used for Figure 4B
rio::export(up_peaks_GOBP_df, "K4me3_over_DEpromsUpGenes_upPeaksOnly_GOBP.xlsx")


