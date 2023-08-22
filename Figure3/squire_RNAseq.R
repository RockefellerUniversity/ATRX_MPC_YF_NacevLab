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
library(rio)
library(GenomicRanges)
library(rtracklayer)
library(rGREAT)
library(org.Mm.eg.db)

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
# run squire 
################################################################################
################################################################################

### Squire Fetch funtion
system(paste("path_to_squire/Fetch.py",
             "-b mm10 -r -g -p 15 -k -f -c -x -v "))


### Squire Clean function
system(paste("path_to_squire/Clean.py",
             "-b mm10 -v"))

fastqs_R1 <- list.files("path_to_FQs", full.names = T, pattern = "_1.fq.gz")

### Squire Map function
squire_map <- function(fastqs_R1){
  r1_temp <- fastqs_R1
  r2_temp <- gsub("_1.fq.gz", "_2.fq.gz", r1_temp)
  
  fastqs_R2 <- list.files("path_to_FQs", full.names = T, pattern = "_2.fq.gz")
  
  if(r2_temp %in% fastqs_R2){
    system2(command  = "path_to_squire/Map.py",
            args = paste("--read1", r1_temp,
                         "--read2", r2_temp,
                         "--read_length 150",
                         "-f path_to_squire_output/squire_fetch",
                         "-b mm10",
                         "-p 20",
                         "-v"),
            stdout = paste0("path_to_map_logs/out_", 
                            basename(r1_temp), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
            stderr = paste0("path_to_map_logs/err_", 
                            basename(r1_temp), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
    )
  }
  
  
}

bplapply(fastqs_R1, squire_map)

### Squire Count function

bams <- list.files("squire_map", pattern = ".bam$")
bams2 <- gsub(".fastq.bam", "", bams)

count_squire <- function(bam){
  system2(command  = "path_to_squire/Count.py",
          args = paste("--read_length 150",
                       "-b mm10",
                       "-p 6",
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

system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 MSC_B10*",
                     "--group2 MSC_WT*",
                     "--condition1 B10",
                     "--condition2 WT",
                     "-s True", 
                     "-N MSC_B10_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_B10vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_B10vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_B10vWT_byFamily")


system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 MSC_G7*",
                     "--group2 MSC_WT*",
                     "--condition1 G7",
                     "--condition2 WT",
                     "-s True", 
                     "-N MSC_G7_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_G7vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_G7vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_G7vWT_byFamily")

system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 MSC_B10*",
                     "--group2 MSC_WT*",
                     "--condition1 B10",
                     "--condition2 WT",
                     # "-s True", # this would make it count only at the sub family level, we want locus so we don't specify this flag
                     "-N MSC_B10_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_B10vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_B10vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_B10vWT_byLocus")


system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 MSC_G7*",
                     "--group2 MSC_WT*",
                     "--condition1 G7",
                     "--condition2 WT",
                     # "-s True", # this would make it count only at the sub family level, we want locus so we don't specify this flag
                     "-N MSC_G7_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_G7vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_G7vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_G7vWT_byLocus")

################################################################################
################################################################################
# get TEs changing in both clones
################################################################################
################################################################################

rna_deseq_B10_wt <- read.table("squire_call_B10vWT_byLocus/DESeq2_TE_only.txt",
                               header = T)

rna_deseq_G7_wt <- read.table("squire_call_G7vWT_byLocus/DESeq2_TE_only.txt",
                              header = T)


rna_rep_merge <- merge(rna_deseq_B10_wt, 
                       rna_deseq_G7_wt, 
                       by = 0,
                       suffix = c("_B10", "_G7")) %>%
  dplyr::rename(Repeat_TE = Row.names) %>%
  arrange(pvalue_B10)

rna_rep_merge_bothP05 <- rna_rep_merge %>%
  dplyr::filter(padj_B10 < 0.05 & padj_G7 < 0.05) %>%
  dplyr::select(-contains("lfcSE"), -baseMean_G7) %>%
  dplyr::rename(baseMean = baseMean_B10)

# this table is Supplementary Table 13
rna_rep_merge_bothP05_UPonly <- rna_rep_merge_bothP05 %>%
  dplyr::filter(log2FoldChange_B10 > 0 & log2FoldChange_G7 > 0)
rio::export(rna_rep_merge_bothP05_UPonly, "bothClones_vsWT_squire_P05_UPonly_byLocus.xlsx")


################################################################################
################################################################################
# get genes near TEs that go up or down
################################################################################
################################################################################


### read in TEs from TE transcripts analysis, turn into GRanges, and annotate with genes

mm10_blacklist <- import.bed("mm10-blacklist.v2.bed") 

te_both <- rio::import("bothClones_vsWT_TElocal_P05_UPandDOWN_byLocus.xlsx")
te_both_up <- te_both %>%
  dplyr::filter(log2FoldChange_B10 > 0 & log2FoldChange_G7 > 0) %>%
  tidyr::separate(full_name, into = c("TE_id", "seqnames", "range", "strand"), sep = "\\:") %>%
  tidyr::separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(combined_ranges = paste(seqnames, start, end, sep = "_"))

te_both_up_gr <- GRanges(seqnames = te_both_up$seqnames, IRanges(start = as.numeric(te_both_up$start), end = as.numeric(te_both_up$end)))
te_both_up_gr <- te_both_up_gr[!te_both_up_gr %over% mm10_blacklist]
te_both_up_gr_anno <- submitGreatJob(te_both_up_gr, species = "mm10", request_interval = 0) %>%
  plotRegionGeneAssociationGraphs
te_both_up_df_anno <- as.data.frame(te_both_up_gr_anno) %>%
  mutate(combined_ranges = paste(seqnames, start, end, sep = "_")) %>%
  left_join(te_both_up %>% dplyr::select(combined_ranges, Repeat_TE), by = "combined_ranges")
te_both_up_genes <- te_both_up_df_anno$gene %>% unique
te_both_up_ids <- AnnotationDbi::select(org.Mm.eg.db, keys = te_both_up_genes, keytype = "SYMBOL", columns = "ENTREZID") %>%
  na.omit %>%
  pull("ENTREZID")

te_both_down <- te_both %>%
  dplyr::filter(log2FoldChange_B10 < 0 & log2FoldChange_G7 < 0) %>%
  tidyr::separate(full_name, into = c("TE_id", "seqnames", "range", "strand"), sep = "\\:") %>%
  tidyr::separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(combined_ranges = paste(seqnames, start, end, sep = "_"))

te_both_down_gr <- GRanges(seqnames = te_both_down$seqnames, IRanges(start = as.numeric(te_both_down$start), end = as.numeric(te_both_down$end)))
te_both_down_gr <- te_both_down_gr[!te_both_down_gr %over% mm10_blacklist]
te_both_down_gr_anno <- submitGreatJob(te_both_down_gr, species = "mm10", request_interval = 0) %>%
  plotRegionGeneAssociationGraphs
te_both_down_df_anno <- as.data.frame(te_both_down_gr_anno) %>%
  mutate(combined_ranges = paste(seqnames, start, end, sep = "_")) %>%
  left_join(te_both_down %>% dplyr::select(combined_ranges, Repeat_TE), by = "combined_ranges")
te_both_down_genes <- te_both_down_df_anno$gene %>% unique
te_both_down_ids <- AnnotationDbi::select(org.Mm.eg.db, keys = te_both_down_genes, keytype = "SYMBOL", columns = "ENTREZID") %>%
  na.omit %>%
  pull("ENTREZID")

# read in dds object from TE transcripts analysis
# includes genes and TEs
dds_te_local <- readRDS("B10_G7_TElocal_dds.rds")
rlog_te_local <- rlog(dds_te_local) %>% 
  assay %>%
  .[, c(7:9, 1:6)]  %>%
  data.frame 

####### get the TEs that go up and are also annotated by genes that go up in RNAseq

rlog_te_local_up <- rlog_te_local %>%
  dplyr::filter(rownames(rlog_te_local) %in% te_both_up_ids) %>%
  rownames_to_column(var = "gene_id") %>%
  mutate(gene_symbol = AnnotationDbi::select(org.Mm.eg.db, keys = gene_id, keytype = "ENTREZID", columns = "SYMBOL") %>% pull(SYMBOL)) %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames("gene_symbol")

# get the genes where expression goes up in the rRNA depleted RNAseq - from TE local calculations
deseq_res_B10 <- read.csv("B10_vs_WT_deseq_TEandGenes.csv", row.names = 1)
deseq_res_B10_up_gene <- deseq_res_B10 %>%
  dplyr::filter(log2FoldChange > 0 & padj < 0.05) %>%
  pull(SYMBOL)
deseq_res_G7 <- read.csv("G7_vs_WT_deseq_TEandGenes.csv", row.names = 1)
deseq_res_G7_up_gene <- deseq_res_G7 %>%
  dplyr::filter(log2FoldChange > 0 & padj < 0.05) %>%
  pull(SYMBOL)

genes_both_DE_up <- deseq_res_B10_up_gene[deseq_res_B10_up_gene %in% deseq_res_G7_up_gene]

# read in gene list of genes that go up in both clones from polyA RNAseq
DE_both_up_polya <- read.table("G7_B10_overlap_p05only_UPgenes_ID.txt") %>%
  .[,1] %>%
  as.character
DE_both_up_polya_symbol <- AnnotationDbi::select(org.Mm.eg.db, keytype = "ENTREZID", keys = DE_both_up_polya, columns = "SYMBOL") %>% pull(SYMBOL)

## get table with genes near TEs that go up and whether these genes go up in RNAseq experiments
sig_up_anno_both <- data.frame(row.names = rownames(rlog_te_local_up),
                               DE_UP_polyA = as.factor(rownames(rlog_te_local_up) %in% DE_both_up_polya_symbol),
                               DE_UP_rRNAdep = as.factor(rownames(rlog_te_local_up) %in% genes_both_DE_up))

# final table with TEs that go up and annotation of whether that TE is annotated by a gene that goes up in either RNAseq
# Supplementary Table 14
te_both_up_df_anno2 <- left_join(te_both_up_df_anno, sig_up_anno_both %>% rownames_to_column(var = "gene"), by = "gene") %>%
  dplyr::filter(!is.na(gene))
rio::export(te_both_up_df_anno2, "annoTable_genes_over_TElocal_UP_DEanno_polyA_and_rRNA.xlsx")

####### get the TEs that go up and are also annotated by genes that go up in RNAseq

rlog_te_local_down <- rlog_te_local %>%
  dplyr::filter(rownames(rlog_te_local) %in% te_both_down_ids) %>%
  rownames_to_column(var = "gene_id") %>%
  mutate(gene_symbol = AnnotationDbi::select(org.Mm.eg.db, keys = gene_id, keytype = "ENTREZID", columns = "SYMBOL") %>% pull(SYMBOL)) %>%
  dplyr::select(-gene_id) %>%
  column_to_rownames("gene_symbol")

# get the genes where expression goes down - from TE local calculations
deseq_res_B10 <- read.csv("B10_vs_WT_deseq_TEandGenes.csv", row.names = 1)
deseq_res_B10_down_gene <- deseq_res_B10 %>%
  dplyr::filter(log2FoldChange < 0 & padj < 0.05) %>%
  pull(SYMBOL)
deseq_res_G7 <- read.csv("G7_vs_WT_deseq/G7_vs_WT_deseq_TEandGenes.csv", row.names = 1)
deseq_res_G7_down_gene <- deseq_res_G7 %>%
  dplyr::filter(log2FoldChange < 0 & padj < 0.05) %>%
  pull(SYMBOL)

genes_both_DE_down <- deseq_res_B10_down_gene[deseq_res_B10_down_gene %in% deseq_res_G7_down_gene]


# read in gene list of genes that go down in both clones from polyA RNAseq
DE_both_down_polya <- read.table("G7_B10_overlap_p05only_DOWNgenes_ID.txt") %>%
  .[,1] %>%
  as.character
DE_both_down_polya_symbol <- AnnotationDbi::select(org.Mm.eg.db, keytype = "ENTREZID", keys = DE_both_down_polya, columns = "SYMBOL") %>% pull(SYMBOL)

## get table with genes near TEs that go down and whether these genes go down in RNAseq experiments
sig_down_anno_both <- data.frame(row.names = rownames(rlog_te_local_down),
                                 DE_down_polyA = as.factor(rownames(rlog_te_local_down) %in% DE_both_down_polya_symbol),
                                 DE_down_rRNAdep = as.factor(rownames(rlog_te_local_down) %in% genes_both_DE_down))

# final table with TEs that go down and annotation of whether that TE is annotated by a gene that goes down in either RNAseq
# Supplementary Table 15
te_both_down_df_anno2 <- left_join(te_both_down_df_anno, sig_down_anno_both %>% rownames_to_column(var = "gene"), by = "gene") %>%
  dplyr::filter(!is.na(gene))
rio::export(te_both_down_df_anno2, "annoTable_genes_over_TElocal_down_DEanno_polyA_and_rRNA.xlsx")

