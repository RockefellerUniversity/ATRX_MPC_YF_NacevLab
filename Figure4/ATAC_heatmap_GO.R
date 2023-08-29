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
up_geneID_p05FC2 <- read.table("DE_gene_lists/G7_B10_overlap_p05_FC2_UPgenes_ID.txt") %>%
  pull(V1) 

down_geneID_p05FC2 <- read.table("DE_gene_lists/G7_B10_overlap_p05_FC2_DOWNgenes_ID.txt") %>%
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


# ATACseq bigwig files
bigwigs <- c("path_to_bigwigs/Sorted_WT_R1Normalised.bw",
             "path_to_bigwigs/Sorted_WT_R2Normalised.bw",
             "path_to_bigwigs/Sorted_WT_R3Normalised.bw",
             "path_to_bigwigs/Sorted_B10_R1Normalised.bw",
             "path_to_bigwigs/Sorted_B10_R2Normalised.bw",
             "path_to_bigwigs/Sorted_B10_R3Normalised.bw",
             "path_to_bigwigs/Sorted_G7_R1Normalised.bw",
             "path_to_bigwigs/Sorted_G7_R2Normalised.bw",
             "path_to_bigwigs/Sorted_G7_R3Normalised.bw")


ranges_tss <- c("mm10_TSS_up.bed",
                "mm10_TSS_down.bed")

cp <- BamBigwig_to_chipProfile(signalFiles = bigwigs,
                               testRanges = ranges_tss,
                               format = "bigwig", 
                               style = "point", 
                               distanceUp = 3000,
                               distanceDown = 3000)

pp <- as_profileplyr(cp)

# clean up profileplyr object
rownames(sampleData(pp)) <- gsub("Normalised.bw", "",sampleData(pp)$sample_labels)
rownames(sampleData(pp)) <- gsub("Sorted_", "",rownames(sampleData(pp)))

mcols(pp)$peak_set <- ifelse(mcols(pp)$sgGroup == "mm10_TSS_down.bed", 
                             "DE_down_peaks", 
                             "DE_up_peaks")

pp <- groupBy(pp, group = "peak_set")

### get mean signal for replicates
pp_temp <- pp

groups <- c("WT", "B10", "G7")
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
mean_pp_atac <- new("profileplyr", tempSE,
                    params=params(pp_temp),
                    sampleData=sampleData)

### get promoters that overlap significant peaks

# read in ATAC peaks that go up
atac_up <- rtracklayer::import("bothClones_UP_P05_ATAC.bed")

# filter out promoter regions with no ATAC signal
mean_pp_atac_sum <- profileplyr::summarize(mean_pp_atac, fun = rowMeans, output = "matrix")
mean_pp_atac <- mean_pp_atac[!rowSums(mean_pp_atac_sum) == 0,,]

# subset to only the promoters of genes that go up in RNAseq
mean_pp_atac_up <- mean_pp_atac[mcols(mean_pp_atac)$sgGroup == "mm10_TSS_up.bed"]

# further subset to those promoters that also overlap an ATACseq peak that has increased signal
mcols(mean_pp_atac_up)$over_sig <- ifelse(mean_pp_atac_up %over% atac_up, "sig_up", "no_change")
mcols(mean_pp_atac_up)$over_sig <- ordered(mcols(mean_pp_atac_up)$over_sig, levels = c("sig_up", "no_change"))
mean_pp_atac_up <- groupBy(mean_pp_atac_up, group = "over_sig")

# this is the heatmap in Figure 4F
pdf("ATAC_over_DEproms_up_meanSignal_groupByupPeaks.pdf", height = 8, width =5)
generateEnrichedHeatmap(mean_pp_atac_up,
                        ylim = "common_max",
                        matrices_pos_line = FALSE,
                        include_group_annotation = T, use_raster = T 
)
dev.off()

################################################################################
################################################################################
# Gene ontology for ATAC peaks over promoters of genes that go up
################################################################################
################################################################################


mean_pp_atac_up_geneID <- mcols(mean_pp_atac_up)$name[mcols(mean_pp_atac_up)$over_sig == "sig_up"]

# get universe from differential analysis
universe_G7 <- read.table("path_to_DEG_output/Group_Atrx_sg5_G7_minus_WTGOIonlyDEG.xls") %>%
  pull(V1) 
universe_B10 <- read.table("path_to_DEG_output/Group_Atrx_sg6_B10_minus_WTGOIonlyDEG.xls") %>%
  pull(V1) 

universe_both <- c(universe_G7, universe_B10) %>% unique


# GOBP
up_peaks_GOBP <- enrichGO(mean_pp_atac_up_geneID,
                          ont = "BP", 
                          universe = universe_both,
                          OrgDb = "org.Mm.eg.db", 
                          pvalueCutoff = 0.5, 
                          qvalueCutoff = 0.5)

up_peaks_GOBP_df <- data.frame(up_peaks_GOBP)

# this table is used for dot plot in figure 4F and is Supplementary Table 21
rio::export(up_peaks_GOBP_df, "ATAC_over_DEpromsUpGenes_upPeaksOnly_GOBP.xlsx")


################################################################################
################################################################################
# Gene ontology for all ATAC peaks that go up in both clones
################################################################################
################################################################################

### read in G7
# get GRanges of peaks at various cutoffs 

G7overWT <- read.table("../../../../DE_Genes/salmon/Antibody_ATAC_fromMACSisBlacklisted/ATAC/Antibody_ATAC_fromMACSisBlacklisted___Group_G7_ATAC_minus_WT_ATACDEG.xls",
                       header = TRUE) %>%
  na.omit

G7overWT_UP_P05 <- G7overWT[G7overWT$padj < 0.05 & G7overWT$log2FoldChange > 0, ]
head(G7overWT_UP_P05)
G7overWT_UP_P05_peaks <- strsplit(G7overWT_UP_P05$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
G7overWT_UP_P05_gr <- GRanges(seqnames = G7overWT_UP_P05_peaks[,2], 
                              IRanges(start = as.numeric(G7overWT_UP_P05_peaks[,3]),
                                      end = as.numeric(G7overWT_UP_P05_peaks[,4])))

G7overWT_UP_P01 <- G7overWT[G7overWT$padj < 0.01 & G7overWT$log2FoldChange > 0, ]
head(G7overWT_UP_P01)
G7overWT_UP_P01_peaks <- strsplit(G7overWT_UP_P01$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
G7overWT_UP_P01_gr <- GRanges(seqnames = G7overWT_UP_P01_peaks[,2], 
                              IRanges(start = as.numeric(G7overWT_UP_P01_peaks[,3]),
                                      end = as.numeric(G7overWT_UP_P01_peaks[,4])))

G7overWT_UP_P001 <- G7overWT[G7overWT$padj < 0.001 & G7overWT$log2FoldChange > 0, ]
head(G7overWT_UP_P001)
G7overWT_UP_P001_peaks <- strsplit(G7overWT_UP_P001$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
G7overWT_UP_P001_gr <- GRanges(seqnames = G7overWT_UP_P001_peaks[,2], 
                               IRanges(start = as.numeric(G7overWT_UP_P001_peaks[,3]),
                                       end = as.numeric(G7overWT_UP_P001_peaks[,4])))

G7overWT_UP_P001FC2 <- G7overWT[G7overWT$padj < 0.001 & G7overWT$log2FoldChange > 1, ]
head(G7overWT_UP_P001FC2)
G7overWT_UP_P001FC2_peaks <- strsplit(G7overWT_UP_P001FC2$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
G7overWT_UP_P001FC2_gr <- GRanges(seqnames = G7overWT_UP_P001FC2_peaks[,2], 
                                  IRanges(start = as.numeric(G7overWT_UP_P001FC2_peaks[,3]),
                                          end = as.numeric(G7overWT_UP_P001FC2_peaks[,4])))

### read in B10
# get GRanges of peaks at various cutoffs

B10overWT <- read.table("../../../../DE_Genes/salmon/Antibody_ATAC_fromMACSisBlacklisted/ATAC/Antibody_ATAC_fromMACSisBlacklisted___Group_B10_ATAC_minus_WT_ATACDEG.xls",
                        header = TRUE) %>%
  na.omit

B10overWT_UP_P05 <- B10overWT[B10overWT$padj < 0.05 & B10overWT$log2FoldChange > 0, ]
head(B10overWT_UP_P05)
B10overWT_UP_P05_peaks <- strsplit(B10overWT_UP_P05$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
B10overWT_UP_P05_gr <- GRanges(seqnames = B10overWT_UP_P05_peaks[,2], 
                               IRanges(start = as.numeric(B10overWT_UP_P05_peaks[,3]),
                                       end = as.numeric(B10overWT_UP_P05_peaks[,4])))

B10overWT_UP_P01 <- B10overWT[B10overWT$padj < 0.01 & B10overWT$log2FoldChange > 0, ]
head(B10overWT_UP_P01)
B10overWT_UP_P01_peaks <- strsplit(B10overWT_UP_P01$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
B10overWT_UP_P01_gr <- GRanges(seqnames = B10overWT_UP_P01_peaks[,2], 
                               IRanges(start = as.numeric(B10overWT_UP_P01_peaks[,3]),
                                       end = as.numeric(B10overWT_UP_P01_peaks[,4])))

B10overWT_UP_P001 <- B10overWT[B10overWT$padj < 0.001 & B10overWT$log2FoldChange > 0, ]
head(B10overWT_UP_P001)
B10overWT_UP_P001_peaks <- strsplit(B10overWT_UP_P001$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
B10overWT_UP_P001_gr <- GRanges(seqnames = B10overWT_UP_P001_peaks[,2], 
                                IRanges(start = as.numeric(B10overWT_UP_P001_peaks[,3]),
                                        end = as.numeric(B10overWT_UP_P001_peaks[,4])))

B10overWT_UP_P001FC2 <- B10overWT[B10overWT$padj < 0.001 & B10overWT$log2FoldChange > 1, ]
head(B10overWT_UP_P001FC2)
B10overWT_UP_P001FC2_peaks <- strsplit(B10overWT_UP_P001FC2$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
B10overWT_UP_P001FC2_gr <- GRanges(seqnames = B10overWT_UP_P001FC2_peaks[,2], 
                                   IRanges(start = as.numeric(B10overWT_UP_P001FC2_peaks[,3]),
                                           end = as.numeric(B10overWT_UP_P001FC2_peaks[,4])))

# merge peaks and gpeaks that exist in both clones
bothClones_consensus <- c(G7overWT_gr, B10overWT_gr) %>%
  reduce
bothClones_UP_P05_gr <- bothClones_consensus[bothClones_consensus %over% G7overWT_UP_P05_gr & bothClones_consensus %over% B10overWT_UP_P05_gr]
bothClones_UP_P01_gr <- bothClones_consensus[bothClones_consensus %over% G7overWT_UP_P01_gr & bothClones_consensus %over% B10overWT_UP_P01_gr]
bothClones_UP_P001_gr <- bothClones_consensus[bothClones_consensus %over% G7overWT_UP_P001_gr & bothClones_consensus %over% B10overWT_UP_P001_gr]
bothClones_UP_P001FC2_gr <- bothClones_consensus[bothClones_consensus %over% G7overWT_UP_P001FC2_gr & bothClones_consensus %over% B10overWT_UP_P001FC2_gr]

## annotate peaks using GREAT

# get universe of genes for great
# we will use the genes from all the high confidence peaks 

great_universe_bothClones <- submitGreatJob(bothClones_consensus, 
                                            species = "mm10")
great_universe_bothClones_gr <- plotRegionGeneAssociationGraphs(great_universe_bothClones)
great_universe_bothClones_geneID <- great_universe_bothClones_gr$gene %>%
  unique() %>%
  AnnotationDbi::select(org.Mm.eg.db, keys = ., keytype = "SYMBOL", columns = "ENTREZID") %>%
  pull(ENTREZID) %>%
  .[!is.na(.)]

# annotate differntial peak sets with GREAT
great_bothClones <- lapply(list(`p < 0.05` = bothClones_UP_P05_gr, 
                                `p < 0.01` = bothClones_UP_P01_gr, 
                                `p < 0.001` = bothClones_UP_P001_gr,
                                `p < 0.001 and FC > 2` = bothClones_UP_P001_gr), 
                           submitGreatJob,
                           species = "mm10",
                           request_interval = 0)

great_bothClones_regions <- lapply(great_bothClones, rGREAT::plotRegionGeneAssociationGraphs, request_interval = 0)

# get gene ids for these genes
great_bothClones_geneid <- lapply(great_bothClones_regions, function(x){
  sym <- x$gene
  ids <- AnnotationDbi::select(org.Mm.eg.db, keys = sym, keytype = "SYMBOL", columns = "ENTREZID") %>%
    pull(ENTREZID)
  ids <- ids[!is.na(ids)] %>%
    unique
  ids
})

# run gene ontology for these gene sets
allgenes_great_bothClones_GOBP <- lapply(great_bothClones_geneid, 
                                         enrichGO,
                                         ont = "BP",
                                         OrgDb = "org.Mm.eg.db", 
                                         universe = great_universe_bothClones_geneID,
                                         pvalueCutoff = 0.5, 
                                         qvalueCutoff = 0.5)

names(allgenes_great_bothClones_GOBP) <- gsub(" < | > | and ", "", names(allgenes_great_bothClones_GOBP))

# this table is Supplementary Table 20 and is used to make Figure 4E
lapply(seq_along(allgenes_great_bothClones_GOBP), function(x) write.csv(allgenes_great_bothClones_GOBP[[x]], paste0("bothClonesoverWT_GOBP_result_", names(allgenes_great_bothClones_GOBP)[x], "_UPpeaks_GREAT.csv")))
