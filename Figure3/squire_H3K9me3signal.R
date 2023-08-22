library(Herper)
library("vsn")
library("pheatmap")
library("DESeq2")
library("BiocParallel")
library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library("hexbin")
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(ComplexHeatmap)

################################################################################
################################################################################
# set up squire environment
################################################################################
################################################################################

Herper::install_CondaTools(c("python=2.7.18"), 
                           env = "squire_yan",  
                           pathToMiniConda = "path_to_miniconda")

Herper::install_CondaTools(tools = c("star=2.5.3a",
                                     "bedtools=2.25.0",
                                     "samtools=1.1",
                                     "stringtie=1.3.3"),
                           env = "squire_yan",
                           updateEnv = TRUE,
                           pathToMiniConda = "path_to_miniconda")

Herper::install_CondaTools(tools = c("ucsc-genepredtobed",
                                     "ucsc-gtftogenepred",
                                     "ucsc-genepredtogtf",
                                     "ucsc-bedgraphtobigwig"),
                           env = "squire_yan",
                           updateEnv = TRUE,
                           pathToMiniConda = "path_to_miniconda")

Herper::install_CondaTools(tools = c("git"),
                           env = "squire_yan",
                           updateEnv = TRUE,
                           pathToMiniConda = "path_to_miniconda")


Herper::install_CondaTools(tools = c("pyfaidx"),
                           env = "squire_yan",
                           updateEnv = TRUE,
                           pathToMiniConda = "path_to_miniconda")

Herper::local_CondaEnv(new = "squire_yan", 
                       pathToMiniConda="path_to_miniconda")

#system("git clone https://github.com/wyang17/SQuIRE; cd SQuIRE; pip install -e .")

################################################################################
################################################################################
# run squire - same workflow used for H3K9me3 and RNAseq
################################################################################
################################################################################

### Squire Fetch funtion
system(paste("path_to_squire/Fetch.py",
             "-b mm10 -r -g -p 15 -k -f -c -x -v "))


### Squire Clean function
system(paste("path_to_squire/Clean.py",
             "-b mm10 -v"))

fastqs <- list.files("path_to_FQs", full.names = T)

### Squire Map function
squire_map <- function(fastq){
  system2(command  = "path_to_squire/Map.py",
          args = paste("--read1", fastq,
                       "--read_length 75",
                       "-f path_to_squire_output/squire_fetch",
                       "-b mm10",
                       "-p 20",
                       "-v"),
          stdout = paste0("path_to_map_logs/out_", 
                          basename(fastq), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
          stderr = paste0("path_to_map_logs/err_", 
                          basename(fastq), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
  )
}

bplapply(fastqs, squire_map)

### Squire Count function

bams <- list.files("squire_map", pattern = ".bam$")
bams2 <- gsub(".fastq.bam", "", bams)

count_squire <- function(bam){
  system2(command  = "path_to_squire/Count.py",
          args = paste("--read_length 75",
                       "-b mm10",
                       "-p 1",
                       "-f path_to_squire_output/squire_fetch",
                       "-c path_to_squire_output/squire_clean",
                       "--name", bam,
                       "-v"),
          stdout = paste0("path_to_count_logs", 
                          bam, "_count_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
          stderr = paste0("path_to_count_logs", 
                          bam, "_count_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
  )
}

bplapply(bams2, count_squire)

### Squire Call function

#####
#k9me3
#####

#B10 vs Wt

system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 B10_k9me3*",
                     "--group2 WT_k9me3*",
                     "--condition1 B10_k9me3",
                     "--condition2 WT_k9me3",
                     "-s True",
                     "-N B10_k9me3_v_WT_k9me3",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_B10_v_WT_k9me3_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_B10_v_WT_k9me3_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_B10_v_WT_k9me3")

#G7 vs Wt

system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 G7_k9me3*",
                     "--group2 WT_k9me3*",
                     "--condition1 G7_k9me3",
                     "--condition2 WT_k9me3",
                     "-s True",
                     "-N G7_k9me3_v_WT_k9me3",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_G7_v_WT_k9me3_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_G7_v_WT_k9me3_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_G7_v_WT_k9me3") 

#G7 vs B10

system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 G7_k9me3*",
                     "--group2 B10_k9me3*",
                     "--condition1 G7_k9me3",
                     "--condition2 B10_k9me3",
                     "-s True",
                     "-N G7_k9me3_v_B10_k9me3",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_G7_v_B10_k9me3_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_G7_v_B10_k9me3_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_G7_v_B10_k9me3") 


################################################################################
################################################################################
# make heatmap for K9me3 signal over TEs
################################################################################
################################################################################

k9_cts_B10_wt <- read.table("squire_call_B10_v_WT_k9me3/B10_k9me3_v_WT_k9me3_gene_subF_counttable.txt",
                            header = T)

k9_cts_B10_wt_rep <- dplyr::filter(k9_cts_B10_wt, !grepl(",", gene_id))

k9_cts_G7_wt <- read.table("squire_call_G7_v_WT_k9me3/G7_k9me3_v_WT_k9me3_gene_subF_counttable.txt",
                           header = T)

k9_cts_G7_wt_rep <- dplyr::filter(k9_cts_G7_wt, !grepl(",", gene_id))

k9_rep_merge <- full_join(k9_cts_B10_wt_rep, 
                          k9_cts_G7_wt_rep, 
                          by = c("gene_id",
                                 "WT_k9me3_R1_S3_R1_001", 
                                 "WT_k9me3_R3_S3_R1_001", 
                                 "WT_k9me3_R2_S24_R1_001")) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix

k9_rep_merge <- k9_rep_merge[, order(colnames(k9_rep_merge))]

colnames(k9_rep_merge) <- gsub("_S.*", "",colnames(k9_rep_merge) )

colData_k9 <- data.frame(row.names = colnames(k9_rep_merge),
                         strsplit(colnames(k9_rep_merge), split = "_") %>%
                           do.call(rbind.data.frame, .))
colnames(colData_k9) <- c("cell_line", "antibody", "replicate")
colData_k9$cell_line <- ordered(colData_k9$cell_line, levels = c("WT", "B10", "G7"))
colData_k9 <- dplyr::arrange(colData_k9, cell_line, replicate)


colData_k9$cell_line <- as.character(colData_k9$cell_line)
k9_rep_merge <- k9_rep_merge[, match(rownames(colData_k9), colnames(k9_rep_merge))]

# make heatmap by summing up counts for each subfamily, then taking rlog 

cts <- k9_rep_merge
cts <- cts[!grepl("\\?", rownames(cts)), ]

colnames(cts) <- gsub("_S.*", "",colnames(cts) )

colData_k9_cts <- data.frame(row.names = colnames(cts),
                             strsplit(colnames(cts), split = "_") %>%
                               do.call(rbind.data.frame, .))
colnames(colData_k9_cts) <- c("cell_line", "antibody", "replicate")
colData_k9_cts$cell_line <- ordered(colData_k9_cts$cell_line, levels = c("WT", "B10", "G7"))
colData_k9_cts <- dplyr::arrange(colData_k9_cts, cell_line, replicate)

cts <- cts[, match(rownames(colData_k9_cts), colnames(cts))]

cts_groupedSF <- as.data.frame(cts) %>% 
  rownames_to_column(var = "TErepeat") %>%
  pivot_longer(names_to = "sample", values_to = "cts", WT_k9me3_R1:G7_k9me3_R3) %>% 
  separate(TErepeat, into = c("repeat", "subfamily", "family"), sep = ":", remove = FALSE) %>%
  group_by(subfamily, sample) %>%
  dplyr::summarise(sum_cts = as.integer(sum(cts))) %>%
  pivot_wider(names_from = "sample", values_from = "sum_cts") %>%
  column_to_rownames(var = "subfamily")

cts_groupedSF <- cts_groupedSF[, match(rownames(colData_k9_cts), colnames(cts_groupedSF))]

colData_k9_cts$cell_line <- as.character(colData_k9_cts$cell_line)
dds_groupedSF_batch <- DESeqDataSetFromMatrix(countData = cts_groupedSF,
                                              colData = colData_k9_cts,
                                              ~cell_line)

keep <- rowSums(counts(dds_groupedSF_batch)) >= 10
table(keep)
dds_groupedSF_batch <- dds_groupedSF_batch[keep,]

rlog_sum_cts <- rlog(dds_groupedSF_batch)
rlog_sum_cts_mat <- assay(rlog_sum_cts)

# remove batch effects, scale, and make heatmap
rlog_sum_cts_mat_batchRm <- limma::removeBatchEffect(rlog_sum_cts_mat, colData_k9_cts$replicate)
rlog_sum_cts_mat_batchRm_scale <- t(scale(t(rlog_sum_cts_mat_batchRm)))

# this heatmap is Figure 3E
pdf("TE_subfamily_heatmap_squire_sumCts_thenRlog_10T_batch.pdf", height = 10)
Heatmap(rlog_sum_cts_mat_batchRm_scale,
        cluster_columns = FALSE)
dev.off()



