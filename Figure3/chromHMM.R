
# load libraries

library(readr)
library(matrixStats)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(BiocParallel)
library(profileplyr)

library(rtracklayer)
library(clusterProfiler)
library(msigdbr)
library(rGREAT)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rio)
library(ChIPseeker)



################################################################################
################################################################################
# run chromHMM
################################################################################
################################################################################

# downloded Chromhmm v1.23 zip and put on HPC here: path_to_software/ChromHMM_v1.23

# make tab separated file for the binarizebam function
setwd("chromHMM_results_directory")

bams <- list.files("path_to_BAMs")
bams_igg <- bams[grepl("IgG", bams)]
bams_igg_noR4 <- bams_igg[!grepl("R4.bam", bams_igg)]

bams_k4me3 <- bams[grepl("k4me3", bams)]
bams_k4me3_select <- bams_k4me3[!bams_k4me3 %in% c("Sorted_B10_k4me3_R1.bam", "Sorted_G7_k4me3_R3.bam", "Sorted_WT_k4me3_R3.bam")]
bams_igg_fork4me3 <- bams_igg[!bams_igg %in% c("Sorted_B10_IgG_R1.bam", "Sorted_G7_IgG_R3.bam", "Sorted_WT_IgG_R3.bam")]

bams_k27ac <- bams[grepl("k27ac", bams)]
bams_k27me3 <- bams[grepl("k27me3", bams)]
bams_k9me3 <- bams[grepl("k9me3", bams)]

chromHMM_design <- data.frame(cell = rep("10T", 36),
                              mark = c(rep("H3k4me3", 9),
                                       rep("H3k27ac", 9),
                                       rep("H3k27me3", 9),
                                       rep("H3k9me3", 9)),
                              file1 = c(bams_k4me3_select,
                                        bams_k27ac,
                                        bams_k27me3,
                                        bams_k9me3),
                              file2 = c(bams_igg_fork4me3,
                                        rep(bams_igg_noR4, 3))
)

# only use WT samples
chromHMM_design_WTonly <- chromHMM_design[grepl("WT", chromHMM_design$file1), ]
write_tsv(chromHMM_design_WTonly, "chromHMM_design_WTonly.txt", col_names = F, quote = "none")

system(paste("java -mx4000M -jar path_to_software/ChromHMM_v1.23/ChromHMM.jar BinarizeBam",
             "path_to_software/ChromHMM_v1.23/CHROMSIZES/mm10.txt",
             "path_to_BAMs", 
             "chromHMM_design_WTonly.txt", 
             "output"))

learn_model <- function(num_states){
  output <- file.path("output", paste0("learn_model_output_",num_states, "states"))
  dir.create(output)
  system(paste("java -mx4000M -jar path_to_software/ChromHMM_v1.23/ChromHMM.jar LearnModel",
               "output",
               output, 
               num_states, 
               "mm10"))
}

# Figure 3a was an output for the 10 state model from this function
bplapply(2:20, learn_model)

# move all the emissions files to one folder for compare models function
dir.create("output/emissions_files_forCompare")
system("cp ./output/**/emissions* output/emissions_files_forCompare/")


compare_model_16 <- function(path_to_emissions){
  output <- file.path(path_to_emissions, "compare_output")
  dir.create(output)
  system(paste("java -mx4000M -jar path_to_software/ChromHMM_v1.23/ChromHMM.jar CompareModels",
               "output/learn_model_output_16states/emissions_16.txt",
               path_to_emissions, 
               file.path(output, paste0("compare_outputvs16"))))
}

compare_model_16("output/emissions_files_forCompare")

################################################################################
################################################################################
# optimize number of states
################################################################################
################################################################################

# Repeat the methods in paper (https://www.nature.com/articles/s41586-020-2093-3) to optimize number of states

emissions_all <- data.frame()
for (i in 2:16){
  emissions_temp <- read_tsv(file.path("output",
                                       paste0("learn_model_output_", i, "states"), 
                                       paste0("emissions_", i, ".txt")))
  colnames(emissions_temp)[1] <- "State"
  emissions_temp$State <- paste(i, emissions_temp$State, sep = "_")
  
  emissions_all <- rbind(emissions_all, emissions_temp)
}

pheatmap(emissions_all %>% column_to_rownames(var = "State"))

# do 100 samples of kmeans clustering and take average between sum squares to total SS ratio
set.seed(1234)
seeds <- sample(x = 1:1000000, size = 100, replace = F)

ss_ratio_df <- data.frame(row.names = 2:16)
for (x in seq_along(seeds)){
  set.seed(seeds[x])
  
  k16 <- kmeans(emissions_all %>% column_to_rownames(var = "State"), 
                centers = 16, )
  
  ss_ratio_temp <- vector()
  for (i in 2:16){
    set.seed(seeds[x])
    k_temp <- kmeans(emissions_all %>% column_to_rownames(var = "State"), 
                     centers = i)
    ratio <- k_temp$betweenss/k_temp$totss
    ss_ratio_temp <- c(ss_ratio_temp, ratio)
  }
  ss_ratio_df[, x] <- ss_ratio_temp
  
}


ss_meanRatio_df <- data.frame(num_clusters = 2:16, 
                              ss_ratio_mean = rowMeans(ss_ratio_df))
ss_meanRatio_df <- ss_meanRatio_df %>%
  mutate(ss_ratio_mean_relative = round(ss_ratio_mean/max(ss_ratio_mean), 3)*100)

# plot for Supplementary Figure 4
ggplot(ss_meanRatio_df, aes(x = num_clusters, y = ss_ratio_mean))  +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_text(aes(label = ss_ratio_mean_relative), hjust = 1,  vjust = -1) +
  ylim(c(0,1.1)) +
  ylab("Between SS / Total SS")
ggsave("compareModels_sumSquareRatio_linePlot_WTonly.pdf", width = 8, height = 5)


################################################################################
################################################################################
# make ranged heatmaps for histone marks over each state of the 10 state model
################################################################################
################################################################################

# get bed files for each state of the 10 state model
state_10_bed <- rtracklayer::import("../output/learn_model_output_10states/10T_10_dense.bed")
mcols(state_10_bed) <- state_10_bed$name # for some reason the IRanges metadata column throws off profileplyr
colnames(mcols(state_10_bed))[1] <- "name"

for (i in seq(10)){
  state_temp <- state_10_bed[state_10_bed$name == i]
  rtracklayer::export(state_temp, paste0("state_bed_files/chromHMM_10states_onlyState", i, ".bed"))
}

# using the bed files from the 10 state model, make a ranged heatmap over these regions for each histone mark
setwd("profileplyr_10states")

bigwigs_k4 <- c(
  "path_to_bigWigs/Sorted_WT_k4me3_R1Normalised.bw",
  "path_to_bigWigs/Sorted_WT_k4me3_R2Normalised.bw",
  "path_to_bigWigs/Sorted_WT_k4me3_R4Normalised.bw",
  "path_to_bigWigs/Sorted_B10_k4me3_R2Normalised.bw",
  "path_to_bigWigs/Sorted_B10_k4me3_R3Normalised.bw",
  "path_to_bigWigs/Sorted_B10_k4me3_R4Normalised.bw",
  "path_to_bigWigs/Sorted_G7_k4me3_R1Normalised.bw",
  "path_to_bigWigs/Sorted_G7_k4me3_R2Normalised.bw",
  "path_to_bigWigs/Sorted_G7_k4me3_R4Normalised.bw")

bigwigs_k27ac <- c("path_to_bigWigs/Sorted_WT_k27ac_R1Normalised.bw",
                   "path_to_bigWigs/Sorted_WT_k27ac_R2Normalised.bw",
                   "path_to_bigWigs/Sorted_WT_k27ac_R3Normalised.bw",
                   "path_to_bigWigs/Sorted_B10_k27ac_R1Normalised.bw",
                   "path_to_bigWigs/Sorted_B10_k27ac_R2Normalised.bw",
                   "path_to_bigWigs/Sorted_B10_k27ac_R3Normalised.bw",
                   "path_to_bigWigs/Sorted_G7_k27ac_R1Normalised.bw",
                   "path_to_bigWigs/Sorted_G7_k27ac_R2Normalised.bw",
                   "path_to_bigWigs/Sorted_G7_k27ac_R3Normalised.bw")

bigwigs_k27me <- c("path_to_bigWigs/Sorted_WT_k27me3_R1Normalised.bw",
                   "path_to_bigWigs/Sorted_WT_k27me3_R2Normalised.bw",
                   "path_to_bigWigs/Sorted_WT_k27me3_R3Normalised.bw",
                   "path_to_bigWigs/Sorted_B10_k27me3_R1Normalised.bw",
                   "path_to_bigWigs/Sorted_B10_k27me3_R2Normalised.bw",
                   "path_to_bigWigs/Sorted_B10_k27me3_R3Normalised.bw",
                   "path_to_bigWigs/Sorted_G7_k27me3_R1Normalised.bw",
                   "path_to_bigWigs/Sorted_G7_k27me3_R2Normalised.bw",
                   "path_to_bigWigs/Sorted_G7_k27me3_R3Normalised.bw")

bigwigs_k9me <- c("path_to_bigWigs/Sorted_WT_k9me3_R1Normalised.bw",
                  "path_to_bigWigs/Sorted_WT_k9me3_R2Normalised.bw",
                  "path_to_bigWigs/Sorted_WT_k9me3_R3Normalised.bw",
                  "path_to_bigWigs/Sorted_B10_k9me3_R1Normalised.bw",
                  "path_to_bigWigs/Sorted_B10_k9me3_R2Normalised.bw",
                  "path_to_bigWigs/Sorted_B10_k9me3_R3Normalised.bw",
                  "path_to_bigWigs/Sorted_G7_k9me3_R1Normalised.bw",
                  "path_to_bigWigs/Sorted_G7_k9me3_R2Normalised.bw",
                  "path_to_bigWigs/Sorted_G7_k9me3_R3Normalised.bw")

all_bw <- list(k4me3 = bigwigs_k4,
               k27ac = bigwigs_k27ac,
               k27me3 = bigwigs_k27me,
               k9me3 = bigwigs_k9me)

# remove state 6 as this seems to be mostly artifact regions
beds <- list.files("state_bed_files/", pattern = ".bed", full.names = T)
beds_no6 <- beds[!beds == "state_bed_files//chromHMM_10states_onlyState6.bed"]

# quantify signal over these regions for each 
# do this in parallel for each type of histone mark  
bplapply(seq(length(all_bw)), function(x, bw, beds){
  library(profileplyr)
  base_bw <- basename(bw[[x]])
  base_bw <- gsub("Sorted_", "", base_bw)
  base_bw <- gsub("Normalised.bw", "", base_bw)
  temp_cp <- BamBigwig_to_chipProfile(bw[[x]],
                                      testRanges = beds,
                                      format = "bigwig",
                                      style = "percentOfRegion",
                                      distanceAround = 100)
  temp_pp <- as_profileplyr(temp_cp)
  saveRDS(temp_pp, paste0("proplyr_", names(bw)[x], "_10states_no6.rds"))
}, bw = all_bw, beds = beds_no6)

# separate all the objects with quantified signal into each state, and we will eventually combine the marks of each state to make a heatmap for each state

setwd("profileplyr_10states")

dir.create("sep_into_states")

# make directories for each state
dir_names <- paste(paste0("chromHMM_10states_onlyState", 1:11))
for(i in seq_along(dir_names)){
  dir.create(file.path("sep_into_states",
                       dir_names[i]))
}

# separate k4me3 object and save each state
k4me3_pp <- readRDS("proplyr_k4me3_10states_no6.rds")
state_beds <- names(table(mcols(k4me3_pp)$sgGroup))
for (i in seq(11)){
  k4me3_pp_temp <- k4me3_pp[mcols(k4me3_pp)$sgGroup == paste0("chromHMM_10states_onlyState", i, ".bed"),]
  saveRDS(k4me3_pp_temp,
          file.path("sep_into_states",
                    paste0("chromHMM_10states_onlyState", i),
                    paste0("k4me3_10states_state", i, ".rds")))
}

# separate k27me3 object and save each state
k27me3_pp <- readRDS("proplyr_k27me3_10states_no6.rds")
state_beds <- names(table(mcols(k27me3_pp)$sgGroup))
for (i in seq(11)){
  k27me3_pp_temp <- k27me3_pp[mcols(k27me3_pp)$sgGroup == paste0("chromHMM_10states_onlyState", i, ".bed"),]
  saveRDS(k27me3_pp_temp,
          file.path("sep_into_states",
                    paste0("chromHMM_10states_onlyState", i),
                    paste0("k27me3_10states_state", i, ".rds")))
}

# separate k27ac object and save each state
k27ac_pp <- readRDS("proplyr_k27ac_10states_no6.rds")
state_beds <- names(table(mcols(k27ac_pp)$sgGroup))
for (i in seq(11)){
  k27ac_pp_temp <- k27ac_pp[mcols(k27ac_pp)$sgGroup == paste0("chromHMM_10states_onlyState", i, ".bed"),]
  saveRDS(k27ac_pp_temp,
          file.path("sep_into_states",
                    paste0("chromHMM_10states_onlyState", i),
                    paste0("k27ac_10states_state", i, ".rds")))
}

# separate k9me3 object and save each state
k9me3_pp <- readRDS("proplyr_k9me3_10states_no6.rds")
state_beds <- names(table(mcols(k9me3_pp)$sgGroup))
for (i in seq(11)){
  k9me3_pp_temp <- k9me3_pp[mcols(k9me3_pp)$sgGroup == paste0("chromHMM_10states_onlyState", i, ".bed"),]
  saveRDS(k9me3_pp_temp,
          file.path("sep_into_states",
                    paste0("chromHMM_10states_onlyState", i),
                    paste0("k9me3_10states_state", i, ".rds")))
}

# make range heatmap for each state
for (i in seq(10)){
  if(!i == 6){
    state_dir <- file.path("sep_into_states",
                           paste0("chromHMM_10states_onlyState", i)) 
    
    #if(!file.exists(file.path(state_dir, paste0("heatmap_meanSignal_allMarks_10states_onlyState", i, ".pdf")))){
    obj_paths <- list.files(state_dir, full.names = T, pattern = paste0("state", i, ".rds"))
    objs <- lapply(obj_paths, readRDS)
    objs_combined <- do.call(c, objs)
    rownames(sampleData(objs_combined)) <- gsub("Sorted_", "", rownames(sampleData(objs_combined)))
    rownames(sampleData(objs_combined)) <- gsub("Normalised.bw", "", rownames(sampleData(objs_combined)))
    sampleData(objs_combined)$chip <- gsub("WT_|G7_|B10_", "", rownames(sampleData(objs_combined)))
    sampleData(objs_combined)$chip <- gsub("_R1|_R2|_R3|_R4", "", sampleData(objs_combined)$chip)
    
    mm10_blacklist <- rtracklayer::import("/rugpfs/fs0/brc/scratch/dbarrows/blacklists/mm10-blacklist.v2.bed")
    objs_combined_noBL <- objs_combined[!objs_combined %over% mm10_blacklist]
    saveRDS(objs_combined_noBL, file.path(state_dir, paste0("proplyr_allMarks_10states_onlyState", i, "_noBL.rds")))
    
    antibodies <- c("k4me3", "k27ac", "k9me3", "k27me3")
    
    # get mean signal across replicates
    mean_pp_all <- list()
    for (j in seq_along(antibodies)){
      pp_temp <- objs_combined_noBL[,,grepl(antibodies[j], sampleData(objs_combined_noBL)$chip)]
      
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
    
    mean_pp_all_combine <- do.call(c, mean_pp_all)
    
    saveRDS(mean_pp_all_combine, file.path(state_dir, paste0("proplyr_meanSignal_allMarks_10states_onlyState", i, "_noBL.rds")))
    
    # this plot for state 3 and state 4 are used for Figure 3B and 3C
    pdf(file.path(state_dir, paste0("heatmap_meanSignal_allMarks_10states_onlyState", i, ".pdf")), width = 12)
    generateEnrichedHeatmap(mean_pp_all_combine,
                            ylim = "chip",
                            color_by_sample_group = "chip",
                            matrices_pos_line = FALSE,
                            include_group_annotation = FALSE, 
                            use_raster = TRUE
    )
    dev.off()
    #}
    
  }
  
  
}

################################################################################
################################################################################
# gene ontology for annotated regions from State 3 of the 10 state model
# only for those regions that overlap with K9me3 down peaks
################################################################################
################################################################################

state_bed <- rtracklayer::import("../../state_bed_files/chromHMM_10states_onlyState3.bed")

# subset chromHMM regions to those that overlap k9me3 differential regions
dir.create("k9_down_overlap")
k9down <- rtracklayer::import("/Users/douglasbarrows/Desktop/BRC/collaborations/Allis/Yan/20211223_YF_CnR_WT_G7_B10/doug_analysis/tiling_seacr_union_k9me3/tiling90/k9_tiling90_k9_seacrP05_reducedSet_DOWNonly.bed")
state_k9down_bed <- k9down[k9down %over% state_bed]

# annotate regions with GREAT
state_k9down_great <- submitGreatJob(state_k9down_bed, 
                                     species = "mm10",
                                     request_interval = 0)

state_k9down_great_gr <- rGREAT::plotRegionGeneAssociationGraphs(state_k9down_great, request_interval = 0)              
sym_k9down <- state_k9down_great_gr$gene
ids_k9down <- AnnotationDbi::select(org.Mm.eg.db, keys = sym_k9down, keytype = "SYMBOL", columns = "ENTREZID") %>%
  pull(ENTREZID) %>%
  .[!is.na(.)] %>%
  unique

gobp_res_k9down <- enrichGO(ids_k9down, 
                            OrgDb = "org.Mm.eg.db",
                            ont = "BP",
                            pvalueCutoff = 0.5,
                            qvalueCutoff = 0.5
)

# This table is used for Supplementary Table 8 and is the basis for the dot plot in Figure 3B
rio::export(as.data.frame(gobp_res_k9down), "k9_down_overlap/GOBP_State3regions_overK9me3DOWN_GREATanno.xlsx")

################################################################################
################################################################################
# gene ontology for annotated regions from State 4 of the 10 state model
# only for those regions that overlap with K9me3 down peaks
################################################################################
################################################################################

state_bed <- rtracklayer::import("../../state_bed_files/chromHMM_10states_onlyState4.bed")

# subset chromHMM regions to those that overlap k9me3 differential regions
dir.create("k9_down_overlap")
k9down <- rtracklayer::import("/Users/douglasbarrows/Desktop/BRC/collaborations/Allis/Yan/20211223_YF_CnR_WT_G7_B10/doug_analysis/tiling_seacr_union_k9me3/tiling90/k9_tiling90_k9_seacrP05_reducedSet_DOWNonly.bed")
state_k9down_bed <- k9down[k9down %over% state_bed]

# annotate regions with GREAT
state_k9down_great <- submitGreatJob(state_k9down_bed, 
                                     species = "mm10",
                                     request_interval = 0)

state_k9down_great_gr <- rGREAT::plotRegionGeneAssociationGraphs(state_k9down_great, request_interval = 0)              
sym_k9down <- state_k9down_great_gr$gene
ids_k9down <- AnnotationDbi::select(org.Mm.eg.db, keys = sym_k9down, keytype = "SYMBOL", columns = "ENTREZID") %>%
  pull(ENTREZID) %>%
  .[!is.na(.)] %>%
  unique

gobp_res_k9down <- enrichGO(ids_k9down, 
                            OrgDb = "org.Mm.eg.db",
                            ont = "BP",
                            pvalueCutoff = 0.5,
                            qvalueCutoff = 0.5
)

# This table is used for Supplementary Table 9 and is the basis for the dot plot in Figure 3C
rio::export(as.data.frame(gobp_res_k9up), "k9_up_overlap/GOBP_State4regions_overK9me3UP_GREATanno.xlsx")
