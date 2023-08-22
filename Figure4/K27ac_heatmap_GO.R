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

# function to collages great GRanges output
collapse_great_id <- function(peaks, species){
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  great <- submitGreatJob(peaks, species = species, request_interval = 0)
  great_df <- plotRegionGeneAssociationGraphs(great) %>%
    data.frame 
  if(grepl("mm", species)){
    symbol2id <- AnnotationDbi::select(org.Mm.eg.db, keys = great_df$gene, keytype = "SYMBOL", columns = "ENTREZID") %>%
      dplyr::filter(!duplicated(SYMBOL))
  }else if (grepl("hg", species)){
    symbol2id <- AnnotationDbi::select(org.Hs.eg.db, keys = great_df$gene, keytype = "SYMBOL", columns = "ENTREZID") %>%
      dplyr::filter(!duplicated(SYMBOL))
  }
  
  great_df_id <- great_df %>%
    left_join(symbol2id, by = c("gene" = "SYMBOL")) %>%
    dplyr::rename(gene_id = ENTREZID)
  great_nested <- nest(great_df_id, data = c(gene, gene_id, distTSS)) %>%
    mutate(gene_id = lapply(.[["data"]], function(x) x$gene_id)) %>%
    mutate(gene = lapply(.[["data"]], function(x) x$gene)) %>%
    mutate(distTSS = lapply(.[["data"]], function(x) x$distTSS)) %>%
    dplyr::select(-data)
  
  great_nested_gr <- GRanges(seqnames = great_nested$seqnames,
                             IRanges(start = great_nested$start,
                                     end = great_nested$end),
                             strand = great_nested$strand,
                             mcols = great_nested %>% dplyr::select(gene, distTSS, gene_id))
  colnames(mcols(great_nested_gr)) <- gsub("mcols.", "", colnames(mcols(great_nested_gr)))
  great_nested_gr$score <- seq(length(great_nested_gr))
  great_nested_gr$name <- lapply(great_nested_gr$gene, paste, collapse = "/")
  great_nested_gr$name <- as.character(great_nested_gr$name)
  great_nested_gr$name_id <- lapply(great_nested_gr$gene_id, paste, collapse = "/")
  great_nested_gr$name_id <- as.character(great_nested_gr$name_id)
  return(great_nested_gr)
}

# get consensus set of H3K27ac peaks

# import high confidence K27ac peaks (peak exists in at least 2 replicates in at least one cell line)
k27ac_hc <- rtracklayer::import("k27ac_HC_Peaks_stringent0.01.bed")

k27ac_hc_great_gr <- collapse_great_id(k27ac_hc, species = "mm10")
mm10_bl <- rtracklayer::import("~/Desktop/BRC/blacklists/mm10-blacklist.v2.bed")
k27ac_hc_great_gr_noBL <- k27ac_hc_great_gr[!k27ac_hc_great_gr %over% mm10_bl]

### categorize peaks as overlapping up genes
# get k27ac peaks over up genes
up_geneID_p05FC2 <- read.table("G7_B10_overlap_p05_FC2_UPgenes_ID.txt") %>%
  pull(V1) 
up_geneSym_p05FC2 <-  AnnotationDbi::select(org.Mm.eg.db, keys = as.character(up_geneID_p05FC2), keytypes = "ENTREZID", columns = "SYMBOL") %>%
  pull(SYMBOL) 

great_up_gene_cat_vec <- vector(length = length(k27ac_hc_great_gr_noBL))
for(i in seq_along(k27ac_hc_great_gr_noBL)){
  temp_genes <- k27ac_hc_great_gr_noBL$name_id[i] %>%
    strsplit(split = "/") %>%
    unlist
  if(length(intersect(up_geneID_p05FC2, temp_genes)) > 0){
    great_up_gene_cat_vec[i] <- "gene_sig_up"
  }else{
    great_up_gene_cat_vec[i] <- "gene_no_change"
  }
}
k27ac_hc_great_gr_noBL$great_up_gene_cat <- great_up_gene_cat_vec

k27ac_hc_great_UPgenes <- k27ac_hc_great_gr_noBL[k27ac_hc_great_gr_noBL$great_up_gene_cat  == "gene_sig_up"]

mm10_proms <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, upstream = 2000, downstream = 50)

k27ac_hc_great_UPgenes_noProm <- k27ac_hc_great_UPgenes[!k27ac_hc_great_UPgenes %over% mm10_proms]
k27ac_hc_great_UPgenes_noProm$name <- paste(seqnames(k27ac_hc_great_UPgenes_noProm),
                                            start(k27ac_hc_great_UPgenes_noProm),
                                            end(k27ac_hc_great_UPgenes_noProm),
                                            sep = "_")


rtracklayer::export(k27ac_hc_great_UPgenes_noProm, "k27ac_HC_Peaks_overUP_DEgenes_noProm.bed")

### categorize peaks as overlapping down genes
# get k27ac peaks over down genes
down_geneID_p05FC2 <- read.table("G7_B10_overlap_p05_FC2_DOWNgenes_ID.txt") %>%
  pull(V1) 
down_geneSym_p05FC2 <-  AnnotationDbi::select(org.Mm.eg.db, keys = as.character(down_geneID_p05FC2), keytypes = "ENTREZID", columns = "SYMBOL") %>%
  pull(SYMBOL) 

great_down_gene_cat_vec <- vector(length = length(k27ac_hc_great_gr_noBL))
for(i in seq_along(k27ac_hc_great_gr_noBL)){
  temp_genes <- k27ac_hc_great_gr_noBL$name_id[i] %>%
    strsplit(split = "/") %>%
    unlist
  if(length(intersect(down_geneID_p05FC2, temp_genes)) > 0){
    great_down_gene_cat_vec[i] <- "gene_sig_down"
  }else{
    great_down_gene_cat_vec[i] <- "gene_no_change"
  }
}
k27ac_hc_great_gr_noBL$great_down_gene_cat <- great_down_gene_cat_vec

k27ac_hc_great_DOWNgenes <- k27ac_hc_great_gr_noBL[k27ac_hc_great_gr_noBL$great_down_gene_cat  == "gene_sig_down"]

mm10_proms <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, upstream = 2000, downstream = 50)

k27ac_hc_great_DOWNgenes_noProm <- k27ac_hc_great_DOWNgenes[!k27ac_hc_great_DOWNgenes %over% mm10_proms]
k27ac_hc_great_DOWNgenes_noProm$name <- paste(seqnames(k27ac_hc_great_DOWNgenes_noProm),
                                              start(k27ac_hc_great_DOWNgenes_noProm),
                                              end(k27ac_hc_great_DOWNgenes_noProm),
                                              sep = "_")


rtracklayer::export(k27ac_hc_great_DOWNgenes_noProm, "k27ac_HC_Peaks_overDOWN_DEgenes_noProm.bed")

# quantify H3K4me3 signal over these peaks that overlap up and down genes (non-promoter)

bigwigs <- c("path_to_bigwigs/Sorted_WT_IgG_R1Normalised.bw",
             "path_to_bigwigs/Sorted_WT_IgG_R2Normalised.bw",
             "path_to_bigwigs/Sorted_WT_IgG_R3Normalised.bw",
             "path_to_bigwigs/Sorted_B10_IgG_R1Normalised.bw", 
             "path_to_bigwigs/Sorted_B10_IgG_R2Normalised.bw",
             "path_to_bigwigs/Sorted_B10_IgG_R3Normalised.bw", 
             "path_to_bigwigs/Sorted_G7_IgG_R1Normalised.bw",
             "path_to_bigwigs/Sorted_G7_IgG_R2Normalised.bw",
             "path_to_bigwigs/Sorted_G7_IgG_R3Normalised.bw",
             "path_to_bigwigs/Sorted_WT_k27ac_R1Normalised.bw",
             "path_to_bigwigs/Sorted_WT_k27ac_R2Normalised.bw",
             "path_to_bigwigs/Sorted_WT_k27ac_R3Normalised.bw",
             "path_to_bigwigs/Sorted_B10_k27ac_R1Normalised.bw",
             "path_to_bigwigs/Sorted_B10_k27ac_R2Normalised.bw",
             "path_to_bigwigs/Sorted_B10_k27ac_R3Normalised.bw",
             "path_to_bigwigs/Sorted_G7_k27ac_R1Normalised.bw",
             "path_to_bigwigs/Sorted_G7_k27ac_R2Normalised.bw",
             "path_to_bigwigs/Sorted_G7_k27ac_R3Normalised.bw")


ranges_h33 <- c("k27ac_HC_Peaks_overUP_DEgenes_noProm.bed",
                "k27ac_HC_Peaks_overDOWN_DEgenes_noProm.bed")

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

mcols(pp)$peak_set <- ifelse(mcols(pp)$sgGroup == "k27ac_HC_Peaks_overDOWN_DEgenes_noProm.bed", 
                             "DE_down_peaks", 
                             "DE_up_peaks")

pp <- groupBy(pp, group = "peak_set")

antibodies <- c("IgG", "k27ac")


### mean signal for replicates

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

# subset object to just the K27 over up genes
mean_pp_k27 <- mean_pp_all[["k27ac"]]
mean_pp_k27_sum <- profileplyr::summarize(mean_pp_k27, fun = rowMeans, output = "matrix")
mean_pp_k27 <- mean_pp_k27[!rowSums(mean_pp_k27_sum) == 0,,]

mean_pp_k27_up <- mean_pp_k27[mcols(mean_pp_k27)$sgGroup == "k27ac_HC_Peaks_overUP_DEgenes_noProm.bed"]
mcols(mean_pp_k27_up)$over_sig <- ifelse(rowRanges(mean_pp_k27_up)$names %in% k27ac_hc_sigUP$range_id, "sig_up", "no_change")
mcols(mean_pp_k27_up)$over_sig <- ordered(mcols(mean_pp_k27_up)$over_sig, levels = c("sig_up", "no_change"))
mean_pp_k27_up <- groupBy(mean_pp_k27_up, group = "over_sig")

# this heatmap is Figure 4C
pdf("k27_over_DE_sig_peaks_notProm_UP_meanSignal_groupByUpPeaks.pdf", height = 8, width =5)
generateEnrichedHeatmap(mean_pp_k27_up,
                        ylim = "common_max",
                        matrices_pos_line = FALSE,
                        include_group_annotation = T, use_raster = T 
)
dev.off()


################################################################################
################################################################################
# Gene ontology for K27ac peaks not over promoters
################################################################################
################################################################################

# get the peaks that are associated wit hthe up peaks for gene ontology
mean_pp_k27_up_peaks <- rowRanges(mean_pp_k27_up)[mcols(mean_pp_k27_up)$over_sig == "sig_up"]
# extract actual peaks, not  profileplyr windows, to get genes for GO
original_up_peaks <- data.frame(peaks = mcols(mean_pp_k27_up_peaks)$name) %>%
  tidyr::separate(col = peaks, into = c("seqnames", "start", "end"))
original_up_peaks_gr <- GRanges(seqnames = original_up_peaks$seqnames, 
                                IRanges(start = as.numeric(original_up_peaks$start), end = as.numeric(original_up_peaks$end)))
original_up_peaks_great <- submitGreatJob(original_up_peaks_gr, species = "mm10", request_interval = 0) %>%
  plotRegionGeneAssociationGraphs()

# get universe of all genes in differential expression analysis
universe_G7 <- read.table("Group_Atrx_sg5_G7_minus_Atrx_sg6_B10GOIonlyDEG.xls") %>%
  pull(V1) 
universe_B10 <- read.table("Group_Atrx_sg6_B10_minus_WTGOIonlyDEG.xls") %>%
  pull(V1) 

universe_both <- c(universe_G7, universe_B10) %>% unique

# up genes
up_peaks_great_genes_id <- AnnotationDbi::select(org.Mm.eg.db, keys = up_peaks_great_genes, keytype = "SYMBOL", columns = "ENTREZID") %>%
  pull(ENTREZID) %>%
  .[!is.na(.)]

# GOBP
up_peaks_great_GOBP <- enrichGO(up_peaks_great_genes_id,
                                ont = "BP", 
                                universe = universe_both,
                                OrgDb = "org.Mm.eg.db", 
                                pvalueCutoff = 0.5, 
                                qvalueCutoff = 0.5)

up_peaks_great_GOBP_df <- data.frame(up_peaks_great_GOBP)

# this table is Supplementary Table 19 and is used for Figure 4D
rio::export(up_peaks_great_GOBP_df, "k27ac_nonProm_overDEup_upPeaksOnly_GOBP.xlsx")
