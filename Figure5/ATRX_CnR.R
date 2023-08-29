
# load libraries

library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)
library(tibble)
library(dplyr)
library(tidyr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(Herper)
library(profileplyr)


################################################################################
################################################################################
# Peak differential analysis
################################################################################
################################################################################

# use DESeq2 to get peaks that go down in the ATRX KO samples

# load in ATRX counts over high confidence peaks
load("ATRX_counts.RData")
myATRXCounts

blacklist <- import("mm10-blacklist.v2.bed")
myATRXCounts <- myATRXCounts[!myATRXCounts %over% blacklist]

groups <- gsub("Sorted_", "", colnames(myATRXCounts)) %>%
  gsub("_ATRX_r.*\\.bam", "", .) 

metaDataFrame <- data.frame(CellLine = groups)
rownames(metaDataFrame) <- colnames(myATRXCounts)
metaDataFrame

deseqATRX <- DESeqDataSetFromMatrix(countData = assay(myATRXCounts), 
                                    colData = metaDataFrame,
                                    design = ~CellLine, 
                                    rowRanges = rowRanges(myATRXCounts))
rownames(deseqATRX) <- paste(seqnames(rowRanges(myATRXCounts)), 
                             start(rowRanges(myATRXCounts)), 
                             end(rowRanges(myATRXCounts)), 
                             sep = "_")

colData(deseqATRX2)$CellLine <- relevel(colData(deseqATRX2)$CellLine, ref = "WT")
deseqATRX2 <- DESeq(deseqATRX2)

res_both <- results(deseqATRX2, contrast = c(0, 1/2, 1/2), format = "GRanges")
res_both <- res_both[order(res_both$pvalue), ]
res_both_anno <- annotatePeak(res_both, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db") %>%
  data.frame
rownames(res_both_anno) <- paste(res_both_anno$seqnames, res_both_anno$start, res_both_anno$end, sep = "_")
res_both_toSave <- as.data.frame(res_both_anno) %>%
  dplyr::select(baseMean:padj, SYMBOL, geneId, distanceToTSS, annotation)

# this is Supplementary Table 25 and is used for Figure 5A 
write.csv(res_both_toSave, "ATRX_stringent_IgGnorm_bothClones_DESeqModel_allPeaks.csv")

both_sig <- res_both[res_both$padj < 0.05]

################################################################################
################################################################################
# Overlap of ATRX binding sites with chromHMM
################################################################################
################################################################################

ATRX_down_DESeq_padj05 <- both_sig[both_sig$log2FoldChange < 0]

state_beds <- list.files("state_bed_files_10states", 
                         full.names = T)

padjP05_state_overlap <- data.frame()
for(i in seq_along(state_beds)){
  
  state_bed_temp <- rtracklayer::import(state_beds[i])
  state_temp <- gsub("state_bed_files_10states/chromHMM_10states_onlyState", "state_", state_beds[i])
  state_temp <- gsub(".bed", "", state_temp)
  
  state_numOverlap <- ATRX_down_DESeq_padj05 %over% state_bed_temp %>% sum
  state_num_peaks <- length(state_bed_temp)
  state_normOverlap <- state_numOverlap/state_num_peaks
  temp_repeat_df <- data.frame(row.names = state_temp,
                               num_overlap = state_numOverlap,
                               num_regions_in_state = state_num_peaks,
                               num_overlap_perRegion = state_normOverlap) 
  padjP05_state_overlap <- rbind(padjP05_state_overlap, temp_repeat_df)
  
}
padjP05_state_overlap_bar <- rownames_to_column(padjP05_state_overlap, var = "state")
padjP05_state_overlap_bar$state <- ordered(padjP05_state_overlap_bar$state, levels = paste0("state_", 1:10))

# this is Figure 5B
ggplot(padjP05_state_overlap_bar %>% dplyr::filter(!state == "state_5"), aes(x = state, y = num_overlap_perRegion)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ggtitle("Overlap of ATRX peaks (padj < 0.05) and chromHMM states")
ggsave("ATRX_down_DESeq_padj05/barplot_ATRX_down_DESeq_padj05_overPerRegion_chromHMM_noState5.pdf")

################################################################################
################################################################################
# Overlap of ATRX binding sites with ATAC, H3K4me3, H3K27ac
################################################################################
################################################################################

# make list of high confidence peak sets (in at least 2 of three replicates for one group)

K27ac_WT_seacr_HC <- rtracklayer::import("individual_HC_peaks/k27ac_WTonly_HC_Peaks_stringent0.01.bed")
K27ac_WT_seacr_HC_anno <- annotatePeak(K27ac_WT_seacr_HC, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")

K27ac_WT_seacr_HC_annoGR <- as.GRanges(K27ac_WT_seacr_HC_anno)
K27ac_WT_seacr_HC_prom1kb <- K27ac_WT_seacr_HC_annoGR[K27ac_WT_seacr_HC_annoGR$annotation == "Promoter (<=1kb)"]
K27ac_WT_seacr_HC_NOTprom1kb <- K27ac_WT_seacr_HC_annoGR[!K27ac_WT_seacr_HC_annoGR$annotation == "Promoter (<=1kb)"]

K4me3_WT_seacr_HC <- rtracklayer::import("individual_HC_peaks/k4me3_WTonly_HC_Peaks_stringent0.01_wRep4.bed")

# this object is a GRanges that shows presence of peaks for each sample across consensus peaks
# here we grab the peaks that overlap at least two WT samples

load("peakOverlaps.RData")
central
ATAC_WT_HC <- central[rowSums(as.data.frame(mcols(central)[, c("WT_R1","WT_R2", "WT_R3")])) >= 2] 

peak_list <- list(
  K27ac_WT_seacr_HC = K27ac_WT_seacr_HC,
  K27ac_NOTprom1kb = K27ac_WT_seacr_HC_NOTprom1kb,
  K4me3_WT_seacr_HC = K4me3_WT_seacr_HC,
  ATAC_WT_HC = ATAC_WT_HC
)

# these venn diagrams are Figure 5C-E
# the tables in this loop are Supplementary Tables 27, 28, 29, and 30
for(i in seq_along(peak_list)){
  
  overlap <- findOverlapsOfPeaks(ATRX_down_DESeq_padj05, peak_list[[i]], maxgap=0)
  
  output_path <- file.path("ATRX_down_DESeq_padj05", names(peak_list)[i])
  dir.create(output_path)
  
  pdf(file.path(output_path, paste0("ATRX_down_DESeq_padj05_over_", names(peak_list)[i],"_venn.pdf")))
  makeVennDiagram(overlap, NameOfPeaks = c("ATRX_down_DESeq_padj05", names(peak_list)[i]), cat.pos = 0, disable.logging = T)
  dev.off()
  
  overlap_peaks <- overlap$peaklist[["ATRX_down_DESeq_padj05///peak_list..i.."]]
  
  if(length(overlap_peaks) > 0){

    overlap_peaks_great <- submitGreatJob(overlap_peaks, request_interval = 0, species = "mm10")
    overlap_peaks_great_df <- plotRegionGeneAssociationGraphs(overlap_peaks_great) %>%
      as.data.frame %>%
      mutate(unique_id = paste(seqnames, start, end, sep = "_")) %>%
      group_by(seqnames, start, end) %>%
      summarise(all_genes = paste(gene, collapse = "/"),
                all_distTSS = paste(distTSS, collapse = "/"))
    
    
    rio::export(overlap_peaks_great_df, file.path(output_path, paste0("ATRX_down_DESeq_padj05_over_", names(peak_list)[i],"_withAnnoGREAT.xlsx")))
  }
  
}



################################################################################
################################################################################
# ATRX signal over genes grouped by expression level
################################################################################
################################################################################

tpm <- read.csv("TPM_PerSampleTable.csv", header = TRUE)
tpm_wt <- tpm %>%
  dplyr::select(IDs, Symbols, contains("WT_")) %>%
  mutate(mean_wt_exp = dplyr::select(tpm, contains("WT_")) %>% rowMeans)

# get geenes with no expression
geneIDs_no_expression <- tpm_wt %>%
  dplyr::filter(mean_wt_exp == 0) %>%
  pull(IDs)

# break the non-zero values into quantiles
quant_no_zero <- quantile(tpm_wt$mean_wt_exp[!tpm_wt$mean_wt_exp == 0], probs = seq(0.2, 1, by = .2))

tpm_wt$exp_cat <- ifelse(tpm_wt$mean_wt_exp == 0, "no_expression",
                         ifelse(tpm_wt$mean_wt_exp > 0 & tpm_wt$mean_wt_exp <= quant_no_zero[1], "quant_20", 
                                ifelse(tpm_wt$mean_wt_exp > quant_no_zero[1] & tpm_wt$mean_wt_exp <= quant_no_zero[2], "quant_40", 
                                       ifelse(tpm_wt$mean_wt_exp > quant_no_zero[2] & tpm_wt$mean_wt_exp <= quant_no_zero[3], "quant_60", 
                                              ifelse(tpm_wt$mean_wt_exp > quant_no_zero[3] & tpm_wt$mean_wt_exp <= quant_no_zero[4], "quant_80", 
                                                     ifelse(tpm_wt$mean_wt_exp > quant_no_zero[4] & tpm_wt$mean_wt_exp <= quant_no_zero[5], "quant_100", "other"))))))

# get gene ids
quant_genes_list <- lapply(levels(as.factor(tpm_wt$exp_cat)), function(x) tpm_wt %>% dplyr::filter(exp_cat == x) %>% pull(IDs))
names(quant_genes_list) <- levels(as.factor(tpm_wt$exp_cat))

# get genomic locations of these genes
genes_mm10 <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
quant_genes_gr_list <- lapply(quant_genes_list, function(x) genes_mm10[genes_mm10$gene_id %in% x])

# write bed files for deeptools 
dir.create("bed_files")
lapply(seq_along(quant_genes_gr_list), function(x){
  rtracklayer::export(quant_genes_gr_list[[x]], paste0("bed_files/mm10_genes_", names(quant_genes_gr_list)[x], ".bed" ))
})

# quantify signal over these genes with deeptools

Herper::local_CondaEnv(new = "env_deeptools")
dir.create("deeptools_logs")

bws <- list.files("path_to_input_norm_bigwigs", pattern = ".bw$", full.names = T)

beds <- list.files("bed_files/", full.names = TRUE)

mat_name <- "dtMat_ATRX_IgGnormalized_3kOut_9kIn.mat"
system2(command = "computeMatrix", 
        args = c("scale-regions",
                 "-R", paste(beds, collapse = " "),
                 "-S", paste(bws, collapse = " "),
                 "-out", mat_name,
                 "--regionBodyLength 9000",
                 "--upstream 3000",
                 "--downstream 3000",
                 "--binSize 50",
                 "--numberOfProcessors 4"), 
        stdout = paste0("deeptools_logs/stdout_", Sys.Date(), ".txt"),
        stderr = paste0("deeptools_logs/stderr_", Sys.Date(), ".txt"))

# read deeptools matrix into profileplyr object to eventually plot
pp_in <- import_deepToolsMat("dtMat_ATRX_IgGnormalized_3kOut_9kIn.mat")

bl_mm10 <- rtracklayer::import("mm10-blacklist.v2.bed")
pp <- pp_in[!pp_in %over% bl_mm10]
rownames(sampleData(pp)) <- gsub("Sorted_|Normalised_log2input", "", rownames(sampleData(pp)))

# reorder samples
pp <- pp[,,c(7:9, 1:6)]

mcols(pp)$dpGroup <- ordered(mcols(pp)$dpGroup, levels = c("mm10_genes_no_expression.bed",
                                                           "mm10_genes_quant_20.bed",
                                                           "mm10_genes_quant_40.bed",
                                                           "mm10_genes_quant_60.bed",
                                                           "mm10_genes_quant_80.bed",
                                                           "mm10_genes_quant_100.bed"))

# the original generateProfilePlot function from profileplyr only facets by sample, we want to facet by peak set and color by sample
# this is a modified function

generateProfilePlot_colorByGroup <- function(object, facetGroup, reps){
  binSize <- sampleData(object)$bin.size[1]
  upstream <- sampleData(object)$upstream[1]
  downstream <- sampleData(object)$downstream[1]
  body <- sampleData(object)$body[1]
  generation_method <- sampleData(object)$generation.method[1]
  
  tidy_pp_dtlyr_data_list <- list()
  for (i in seq_along(assays(object))) {
    pp_dtlyr_matrixdata <- as.data.frame(assay(object[, , i]))
    colnames(pp_dtlyr_matrixdata) <- c(1:ncol(assay(object[, 
                                                           , i])))
    pp_dtlyr_rowdata <- as.data.frame(mcols(object))
    pp_dtlyr_data <- merge(pp_dtlyr_matrixdata, pp_dtlyr_rowdata, 
                           by = "row.names", all = T)
    pp_dtlyr_summary <- data.frame(bin = as.numeric(colnames(pp_dtlyr_data)[2:(ncol(assay(object[, 
                                                                                                 , i])) + 1)]))
    facetGroup_column <- pp_dtlyr_data[[facetGroup]] %>% as.factor
    facetGroup_column <- droplevels(facetGroup_column)
    for (x in seq_along(levels(facetGroup_column))) {
      pp_dtlyr_summary[, x + 1] <- colMeans(pp_dtlyr_data[facetGroup_column %in% 
                                                            levels(facetGroup_column)[x], 2:(ncol(assay(object[, , i])) + 1)])
    }
    colnames(pp_dtlyr_summary) <- c("bin", levels(facetGroup_column))
    pp_dtlyr_summary$sample <- rep(rownames(sampleData(object))[i], nrow(pp_dtlyr_summary))
    tidy_pp_dtlyr_data_list[[i]] <- pivot_longer(pp_dtlyr_summary, 
                                                 cols = 2:(length(levels(facetGroup_column)) + 1), 
                                                 names_to = "groups", values_to = "signal")
  }
  tidy_pp_dtlyr_data_list_all <- do.call(rbind, tidy_pp_dtlyr_data_list)
  tidy_pp_dtlyr_data_list_all$groups <- ordered(tidy_pp_dtlyr_data_list_all$groups,
                                                levels = levels(droplevels(mcols(object)[[facetGroup]])))
  
  tidy_pp_dtlyr_data_list_all$groups <- droplevels(tidy_pp_dtlyr_data_list_all$groups)
  tidy_pp_dtlyr_data_list_all$sample <- ordered(tidy_pp_dtlyr_data_list_all$sample, 
                                                levels = rownames(sampleData(object)))
  
  if(reps){
    tidy_pp_dtlyr_data_list_all <- tidy_pp_dtlyr_data_list_all %>%
      tidyr::separate(sample, into = c("sample", "rep"), sep="_(?=[^_]+$)") %>%
      group_by(bin, sample, groups) %>%
      summarise(signal = mean(signal))
    tidy_pp_dtlyr_data_list_all$sample <- ordered(tidy_pp_dtlyr_data_list_all$sample, 
                                                  levels = strsplit(rownames(sampleData(object)), split ="_") %>% lapply(function(x) paste(x[1], x[2], sep = "_")) %>% unlist %>% unique)
  }
  
  
  
  if (body == 0) {
    labels <- c(-upstream, 0, downstream)
    breaks <- c(0, upstream/binSize, (upstream + downstream)/binSize)
    hjust <- c(0, 0.5, 1)
  }else {
    if (upstream == 0 & downstream == 0) {
      labels <- c("start", "end")
      breaks = c(0, body * binSize)
      hjust <- c(0.5, 0.5)
    }
    else if (upstream == 0) {
      labels <- c("start", "end", downstream)
      breaks = c(0, body * binSize, downstream + body * 
                   binSize)
      hjust <- c(0.5, 0.5, 1)
    }
    else if (downstream == 0) {
      labels <- c(-upstream, "start", "end")
      breaks = c(0, upstream, upstream + 2 * (body * 
                                                binSize))
      hjust <- c(0, 0.5, 0.5)
    }
    else {
      labels <- c("-3kb", "TSS", "TES", "+3kb")
      #labels <- c(-(upstream), "TSS", "TES", downstream)
      # breaks = c(0, upstream, upstream + body * binSize, 
      #            upstream + body * binSize + downstream)
      
      # had to edit this, didnt quite work for deeptools input here
      breaks = c(0, upstream/binSize, (upstream + body)/binSize, 
                 (upstream + body + downstream)/binSize)
      hjust <- c(0, 0.5, 0.5, 1)
    }
  }
  
  ggplot(tidy_pp_dtlyr_data_list_all, aes(x = .data$bin, y = .data$signal, color = .data$groups)) + 
    geom_line(size = 0.7) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), 
          axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10, hjust = hjust), 
          legend.title = element_blank(), legend.text = element_text(size = 10), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          strip.text.x = element_text(size = 10)) + 
    #scale_color_manual(values = c("#8F2320", "#1C6D10")) + 
    scale_x_continuous(labels = labels, breaks = breaks) + 
    expand_limits(y = 0.06) + ylab(NULL) + xlab(NULL) +
    geom_vline(xintercept = c(breaks[2], breaks[3]), linetype = "dashed") +
    geom_vline(xintercept = c(breaks[1], breaks[4]), linetype = "dashed", color = "gray") + 
    facet_wrap(~sample, nrow = 1, scales = "fixed")
  
}

# run this function on the profileplyr object 
generateProfilePlot_colorByGroup(object = pp, facetGroup = "dpGroup", reps = T)

# this is Figure 5F
ggsave("profileplots/profileplot_ATRX_IgGnormSignal_over_geneBody_2kbOut_8kbIn_geneExpCats.pdf", width = 12, height = 4)

################################################################################
################################################################################
# TE transcripts
################################################################################
################################################################################

#### run TE transcripts

# set up conda environment
Herper::install_CondaTools(tools = c("python=3.7.7", "pysam"),
                           env = "yan_tetranscripts",
                           pathToMiniConda = "/rugpfs/fs0/brc/scratch/dbarrows/miniconda")

Herper::install_CondaTools(tools = c("star"),
                           env = "yan_tetranscripts",
                           pathToMiniConda = "/rugpfs/fs0/brc/scratch/dbarrows/miniconda", updateEnv = TRUE)
                           
local_CondaEnv("yan_tetranscripts",pathToMiniConda="/rugpfs/fs0/brc/scratch/dbarrows/miniconda")
                           
# make index 

star_index <- function(fasta, gtf, overhang, threads, genomeDir){
  system2("STAR", 
          args = c("--runMode genomeGenerate",
                   "--runThreadN", threads,
                   "--genomeDir", genomeDir,
                   "--genomeFastaFiles", fasta,
                   "--sjdbGTFfile", gtf,
                   "--sjdbOverhang", overhang) # this is supposed to be read length -1
  )
}


# run on brc node 
fasta_bs_mm10 <- "BSgenome.Mmusculus.UCSC.mm10.fa"
gtf_known_mm10 <- "TxDb.Mmusculus.UCSC.mm10.knownGene.gtf"
genomeDir <- "genomes_indexes/mm10/BSgenome_STAR_RL75"

star_index(fasta_bs_mm10, gtf = gtf_known_mm10, overhang = 74, threads = 10, genomeDir = genomeDir)

# align
star_align <- function(fastq_r1, threads, genomeDir, out_dir_parent){
  
  r1_temp <- fastq_r1
  
  out_dir <- gsub("FQs/", "", r1_temp)
  out_dir <- gsub("_R1_001.fastq.gz", "", out_dir)
  dir.create(file.path(out_dir_parent, out_dir))
  out_prefix <- paste0(file.path(out_dir_parent, out_dir), "/", out_dir, "_")
  system2("STAR", 
          args = c("--runMode alignReads",
                   "--runThreadN", threads,
                   "--genomeDir", genomeDir,
                   "--readFilesIn", r1_temp,
                   "--readFilesCommand", "zcat",
                   "--winAnchorMultimapNmax 100",
                   "--outFilterMultimapNmax 100",
                   "--outFileNamePrefix", out_prefix,
                   "--outSAMtype BAM SortedByCoordinate"
          ), 
          stdout = paste0("logs/STAR_stdout_", 
                          basename(out_dir), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
          stderr = paste0("logs/STAR_stderr_", 
                          basename(out_dir), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S")))
  
}

genomeDir <- "genomes_indexes/mm10/BSgenome_STAR_RL75"
fastqs_R1 <- list.files("FQs", 
                        full.names = T)
out_dir_parent <- "STAR_align"

bplapply(fastqs_R1, star_align, threads = 1, genomeDir = genomeDir, out_dir_parent = out_dir_parent)

# run te transcripts function
te_transcripts_run <- function(treat_bams,control_bams, gtf, te_gtf, project_name){
  treat_bams_str <- paste(treat_bams, collapse = " ")
  control_bams_str <- paste(control_bams, collapse = " ")
  system2("TEtranscripts", 
          args = c("-t", treat_bams_str,
                   "-c",control_bams_str,
                   "--GTF", gtf,
                   "--TE", te_gtf,
                   "--format BAM",
                   "--sortByPos",
                   "--mode multi",
                   "--project", project_name
          ), 
          stdout = paste0("logs/TRtr_stdout_", 
                          basename(project_name), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
          stderr = paste0("logs/TEtr_stderr_", 
                          basename(project_name), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S")))
}

gtf_known_mm10 <- "TxDb.Mmusculus.UCSC.mm10.knownGene.gtf"
te_ucsc <- "mm10_rmsk_TE.gtf"

# B10 vs WT
treat_bams <- list.files("STAR_align/", recursive=T, full.names = T) %>% 
  .[grepl("B10_ATRX", .)] %>% 
  .[grepl(".bam", .)]
control_bams <- list.files("STAR_align/", recursive=T, full.names = T) %>% 
  .[grepl("WT_ATRX", .)] %>% 
  .[grepl(".bam", .)]
te_transcripts_run(treat_bams = treat_bams,
                   control_bams= control_bams,
                   gtf = gtf_known_mm10,
                   te_gtf = te_ucsc,
                   project_name = "B10vWT_TE")

# G7 vs WT
treat_bams <- list.files("/rugpfs/fs0/brc/scratch/dbarrows/Allis/Yan/20221028_YF_CnR_ATRX_rep134_fix_nucExt/te_transcripts/STAR_align/", recursive=T, full.names = T) %>% 
  .[grepl("G7_ATRX", .)] %>% 
  .[grepl(".bam", .)]
control_bams <- list.files("/rugpfs/fs0/brc/scratch/dbarrows/Allis/Yan/20221028_YF_CnR_ATRX_rep134_fix_nucExt/te_transcripts/STAR_align/", recursive=T, full.names = T) %>% 
  .[grepl("WT_ATRX", .)] %>% 
  .[grepl(".bam", .)]
te_transcripts_run(treat_bams = treat_bams,
                   control_bams= control_bams,
                   gtf = gtf_known_mm10,
                   te_gtf = te_ucsc,
                   project_name = "/rugpfs/fs0/brc/scratch/dbarrows/Allis/Yan/20221028_YF_CnR_ATRX_rep134_fix_nucExt/te_transcripts/te_transcripts_family/G7vWT_TE")


Herper::local_CondaEnv("yan_tetranscripts",pathToMiniConda="miniconda")

########### Make heatmaps with TEtranscripts results 

# read in TEtranscripts output counts table
atrx_cts_B10_wt <- read.table("B10vWT_TE.cntTable",
                              header = T)

# clean up messy column names
colnames(atrx_cts_B10_wt) <- gsub("X.rugpfs.fs0.brc.scratch.dbarrows.Allis.Yan.20221028_YF_CnR_ATRX_rep134_fix_nucExt.te_transcripts.STAR_align..10T_", "", colnames(atrx_cts_B10_wt))
colnames(atrx_cts_B10_wt) <- gsub("_ab_0_1_1min_NE", "",colnames(atrx_cts_B10_wt))
colnames(atrx_cts_B10_wt) <- gsub("_ab_0.1_1min_NE_", "_1st_",colnames(atrx_cts_B10_wt))
colnames(atrx_cts_B10_wt) <- gsub("_S.*", "",colnames(atrx_cts_B10_wt))

atrx_cts_B10_wt_rep <- dplyr::filter(atrx_cts_B10_wt, grepl(":", gene.TE))

# read in TEtranscripts output counts table
atrx_cts_G7_wt <- read.table("G7vWT_TE.cntTable",
                             header = T)

# clean up messy column names
colnames(atrx_cts_G7_wt) <- gsub("X.rugpfs.fs0.brc.scratch.dbarrows.Allis.Yan.20221028_YF_CnR_ATRX_rep134_fix_nucExt.te_transcripts.STAR_align..10T_", "", colnames(atrx_cts_G7_wt))
colnames(atrx_cts_G7_wt) <- gsub("_ab_0_1_1min_NE", "",colnames(atrx_cts_G7_wt))
colnames(atrx_cts_G7_wt) <- gsub("_ab_0.1_1min_NE_", "_1st_",colnames(atrx_cts_G7_wt))
colnames(atrx_cts_G7_wt) <- gsub("_S.*", "",colnames(atrx_cts_G7_wt))

# select TEs which have colon in name
atrx_cts_G7_wt_rep <- dplyr::filter(atrx_cts_G7_wt, grepl(":", gene.TE))

atrx_rep_merge <- full_join(atrx_cts_B10_wt_rep, 
                            atrx_cts_G7_wt_rep, 
                            by = c("gene.TE",
                                   "WT_ATRX_3rd", 
                                   "WT_ATRX_4th", 
                                   "WT_ATRX_1st")) %>%
  column_to_rownames(var = "gene.TE") %>%
  as.matrix

atrx_rep_merge <- atrx_rep_merge[, order(colnames(atrx_rep_merge))]

# group the TEs by subfamily
ctsSum_groupedSF <- as.data.frame(cts) %>% 
  rownames_to_column(var = "TErepeat") %>%
  pivot_longer(names_to = "sample", values_to = "cts",  WT_ATRX_1st:G7_ATRX_4th) %>% 
  separate(TErepeat, into = c("repeat", "subfamily", "family"), sep = ":", remove = FALSE) %>%
  group_by(subfamily, sample) %>%
  mutate(n = n()) %>%
  dplyr::filter(n > 1) %>% # there are a lot of repeat families with only one repeat and don't include groups we are interested in
  dplyr::select(-n) %>% 
  dplyr::summarise(sum_cts = as.integer(sum(cts))) %>%
  pivot_wider(names_from = "sample", values_from = "sum_cts") %>%
  column_to_rownames(var = "subfamily")

ctsSum_groupedSF <- ctsSum_groupedSF[, match(rownames(colData_atrx_cts), colnames(ctsSum_groupedSF))]

colData_atrx_cts$cell_line <- as.character(colData_atrx_cts$cell_line)
dds_ctsSum <- DESeqDataSetFromMatrix(countData = ctsSum_groupedSF,
                                     colData = colData_atrx_cts,
                                     ~cell_line)

keep <- rowSums(counts(dds_ctsSum)) >= 10
table(keep)
dds_ctsSum <- dds_ctsSum[keep,]

# control for batch
rlog_sum_cts <- rlog(dds_ctsSum)
rlog_sum_cts_mat <- assay(rlog_sum_cts)
rlog_sum_cts_mat_batchRm <- limma::removeBatchEffect(rlog_sum_cts_mat, colData_atrx_cts$replicate)
rlog_sum_cts_mat_batchRm_scale <- t(scale(t(rlog_sum_cts_mat_batchRm)))

# this is Figure 5H
pdf("sum_cts/TE_subfamily_heatmapCts_thenRlog_10T_batch_TEtransc.pdf", height = 10)
Heatmap(rlog_sum_cts_mat_batchRm_scale,
        cluster_columns = FALSE)
dev.off()
