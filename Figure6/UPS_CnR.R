
# load libraries

library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)
library(tibble)
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(rGREAT)
library(Herper)

library(vsn)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(hexbin)

library(BSgenome.Hsapiens.UCSC.hg38)
library(xlsx)

################################################################################
################################################################################
# Gene ontology for K27ac peaks that go up
################################################################################
################################################################################

# load in peaks 
sig_peaks_p01 <- import("k27ac_stringent_clone6andclone4_overlap_sigPeaks_p01_up.bed")

# annotate peaks with GREAT
great_all <- lapply(list(`p < 0.01` = sig_peaks_p01), 
                    submitGreatJob,
                    species = "hg38",
                    request_interval = 0)

great_all_regions <- lapply(great_all, rGREAT::plotRegionGeneAssociationGraphs, request_interval = 0)

great_all_geneid <- lapply(great_all_regions, function(x){
  sym <- x$gene
  ids <- AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "SYMBOL", columns = "ENTREZID") %>%
    pull(ENTREZID)
  ids <- ids[!is.na(ids)] %>%
    unique
  ids
})

# GO-BP using the GREAT gene lists
dir.create("GO_great_GOBP")
allgenes_great_gobp <- lapply(great_all_geneid, 
                              OrgDb = 'org.Hs.eg.db',
                              enrichGO,
                              ont = "BP",  
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5)

names(allgenes_great_gobp) <- gsub(" < ", "", names(allgenes_great_gobp))

# this table is Supplementary Table 33 and is used for Figure 6C
lapply(seq_along(allgenes_great_gobp), function(x){
  if(!is.null(allgenes_great_gobp[[x]])){
    rio::export(allgenes_great_gobp[[x]],
                paste0("GO_great_GOBP/GOBP_great_allgenes_result_", 
                       names(allgenes_great_gobp)[x],
                       "_k27ac_UPpeaks_overIgG.xlsx"))
  }
}) 

################################################################################
################################################################################
# Gene ontology for K9me3 peaks that go up or down
################################################################################
################################################################################

# DOWN peaks
# load in peaks 
sig_peaks_p05 <- import("k9_tiling90_k9_seacrP05_overIgG_reducedSet_DOWNonly.bed")

# annotate peaks with GREAT
great_all <- lapply(list(`p < 0.05` = sig_peaks_p05), 
                    submitGreatJob,
                    species = "hg38",
                    request_interval = 0)

great_all_regions <- lapply(great_all, rGREAT::plotRegionGeneAssociationGraphs, request_interval = 0)

great_all_geneid <- lapply(great_all_regions, function(x){
  sym <- x$gene
  ids <- AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "SYMBOL", columns = "ENTREZID") %>%
    pull(ENTREZID)
  ids <- ids[!is.na(ids)] %>%
    unique
  ids
})

# GO-BP using the GREAT gene lists
dir.create("GO_great_GOBP")

allgenes_great_gobp <- lapply(great_all_geneid, 
                              OrgDb = 'org.Hs.eg.db',
                              enrichGO,
                              ont = "BP",  
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5)

names(allgenes_great_gobp) <- gsub(" < ", "", names(allgenes_great_gobp))

# this table is Supplementary Table 31 and is used for Figure 6B
lapply(seq_along(allgenes_great_gobp), function(x) rio::export(allgenes_great_gobp[[x]], 
                                                               paste0("GO_great_GOBP/GOBP_great_allgenes_result_", 
                                                                      names(allgenes_great_gobp)[x],
                                                                      "_tile90_k9me3_DOWNpeaks_overIgG.xlsx")))

# UP peaks
# load in peaks 
sig_peaks_p05 <- import("k9_tiling90_k9_seacrP05_overIgG_reducedSet_UPonly.bed")

# annotate peaks with GREAT
great_all <- lapply(list(`p < 0.05` = sig_peaks_p05), 
                    submitGreatJob,
                    species = "hg38",
                    request_interval = 0)

great_all_regions <- lapply(great_all, rGREAT::plotRegionGeneAssociationGraphs, request_interval = 0)

great_all_geneid <- lapply(great_all_regions, function(x){
  sym <- x$gene
  ids <- AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "SYMBOL", columns = "ENTREZID") %>%
    pull(ENTREZID)
  ids <- ids[!is.na(ids)] %>%
    unique
  ids
})

# GO-BP using the GREAT gene lists
dir.create("GO_great_GOBP")

allgenes_great_gobp <- lapply(great_all_geneid, 
                              OrgDb = 'org.Hs.eg.db',
                              enrichGO,
                              ont = "BP",  
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5)

names(allgenes_great_gobp) <- gsub(" < ", "", names(allgenes_great_gobp))

# # this table is Supplementary Table 32 and is used for Supplementary Figure 11C
lapply(seq_along(allgenes_great_gobp), function(x) rio::export(allgenes_great_gobp[[x]], 
                                                               paste0("GO_great_GOBP/GOBP_great_allgenes_result_", 
                                                                      names(allgenes_great_gobp)[x],
                                                                      "_tile90_k9me3_UPpeaks_overIgG.xlsx")))


################################################################################
################################################################################
# Squire TE analysis
################################################################################
################################################################################

##########################################
# set up squire environment
##########################################

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

system("git clone https://github.com/wyang17/SQuIRE; cd SQuIRE; pip install -e .")


##########################################
# run squire
##########################################

### Squire fetch funtion
system(paste("path_to_squire/Fetch.py",
             "-b hg38 -r -g -p 15 -k -f -c -x -v "))


### Squire clean function
system(paste("path_to_squire/Clean.py",
             "-b hg38 -v"))


#### Squire map function
fastqs <- list.files("path_to_FQs", full.names = T, pattern = "pairs1.*fq.gz")

squire_map <- function(fastqs_R1){
  r1_temp <- fastqs_R1
  r2_temp <- gsub("pairs1", "pairs2", r1_temp)
  print(r1_temp)
  print(r2_temp)
  fastqs_R2 <- list.files("path_to_FQs", full.names = T, pattern = "pairs2")
  
  if(r2_temp %in% fastqs_R2){
    system2(command  = "path_to_squire/Map.py",
            args = paste("--read1", r1_temp,
                         "--read2", r2_temp,
                         "--read_length 100",
                         "-f path_to_squire_output/squire_fetch",
                         "-b hg38",
                         "-p 10",
                         "-v"),
            stdout = paste0("path_to_map_logs/out_", 
                            basename(r1_temp), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
            stderr = paste0("path_to_map_logs/err_", 
                            basename(r1_temp), "_", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
    )
  }
  
  
}

##### Squire count function

bams <- list.files("squire_map", pattern = ".bam$")
bams2 <- gsub("_pairs1.fq.bam", "", bams)

count_squire <- function(bam){
  system2(command  = "path_to_squire/Count.py",
          args = paste("--read_length 100",
                       "-b hg38",
                       "-p 9",
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


##### Squire call function 

system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 sg1_clone4*",
                     "--group2 WT*",
                     "--condition1 sg1_clone4",
                     "--condition2 WT",
                     "-s True", 
                     "-N sg1_clone4_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_clone4vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_clone4vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_clone4vWT_byFamily")


system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 sg1_clone6*",
                     "--group2 WT*",
                     "--condition1 sg1_clone6",
                     "--condition2 WT",
                     "-s True", 
                     "-N sg1_clone6_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_clone6vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_clone6vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_clone6vWT_byFamily")

system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 sg1_clone4*",
                     "--group2 WT*",
                     "--condition1 sg1_clone4",
                     "--condition2 WT",
                     "-N sg1_clone4_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_clone4vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_clone4vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_clone4vWT_byLocus")


system2(command  = "path_to_squire/Call.py",
        args = paste("--group1 sg1_clone6*",
                     "--group2 WT*",
                     "--condition1 sg1_clone6",
                     "--condition2 WT",
                     "-N sg1_clone6_v_WT",
                     "-f pdf",
                     "-p 20",
                     "-v"),
        stdout = paste0("path_to_call_logs/", 
                        "squire_call_clone6vWT_out", format(Sys.time(), "%Y%m%d_%H.%M.%S")),
        stderr = paste0("path_to_call_logs/", 
                        "squire_call_clone6vWT_err", format(Sys.time(), "%Y%m%d_%H.%M.%S"))
)

file.rename(from = "squire_call", to = "squire_call_clone6vWT_byLocus")

# find differential TEs using DESeq2

rna_deseq_clone4_wt <- read.table("squire_call_clone4vWT_byLocus/DESeq2_TE_only.txt",
                                  header = T)

rna_deseq_clone6_wt <- read.table("squire_call_clone6vWT_byLocus/DESeq2_TE_only.txt",
                                  header = T)


rna_rep_merge <- merge(rna_deseq_clone4_wt, 
                       rna_deseq_clone6_wt, 
                       by = 0,
                       suffix = c("_clone4", "_clone6")) %>%
  dplyr::rename(Repeat_TE = Row.names) %>%
  arrange(pvalue_clone4)


rna_rep_merge_bothP05 <- rna_rep_merge %>%
  dplyr::filter(padj_clone4 < 0.05 & padj_clone6 < 0.05) %>%
  dplyr::select(-contains("lfcSE"), -baseMean_clone6) %>%
  dplyr::rename(baseMean = baseMean_clone4)

# this table is Supplementary Table 36 and is used to make Figure 6F
rna_rep_merge_bothP05_UPonly <- rna_rep_merge_bothP05 %>%
  dplyr::filter(log2FoldChange_clone4 > 0 & log2FoldChange_clone6 > 0)
rio::export(rna_rep_merge_bothP05_UPonly, "UPS_bothClones_vsWT_squire_P05_UPonly_byLocus.xlsx")

################################################################################
################################################################################
# MEME motif analysis for ATAC-seq UP peaks
################################################################################
################################################################################

Herper::install_CondaTools(tools = "meme", 
                           env = "yan_meme", 
                           pathToMiniConda="path_to_miniconda")
Herper::local_CondaEnv(new = "yan_meme",
                       pathToMiniConda="path_to_miniconda")


# get consensus UP peaks for both clones

# make GRanges for clone 4 peaks that go up compared to WT
C4overWT <- read.table("../../DE_Genes/Antibody_ATAC_fromMACSisBlacklisted___Group_C4_ATAC_minus_WT_ATACDEG.xls",
                       header = TRUE) %>%
  na.omit
C4overWT_peaks <- strsplit(C4overWT$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
C4overWT_gr <- GRanges(seqnames = C4overWT_peaks[,2], 
                       IRanges(start = as.numeric(C4overWT_peaks[,3]),
                               end = as.numeric(C4overWT_peaks[,4])))

# get GRanges for clone 4 UP peaks
C4overWT_UP_P05 <- C4overWT[C4overWT$padj < 0.05 & C4overWT$log2FoldChange > 0, ]
head(C4overWT_UP_P05)
C4overWT_UP_P05_peaks <- strsplit(C4overWT_UP_P05$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
C4overWT_UP_P05_gr <- GRanges(seqnames = C4overWT_UP_P05_peaks[,2], 
                              IRanges(start = as.numeric(C4overWT_UP_P05_peaks[,3]),
                                      end = as.numeric(C4overWT_UP_P05_peaks[,4])))

# make GRanges for clone 6 peaks that go up compared to WT
C6overWT <- read.table("../../DE_Genes/Antibody_ATAC_fromMACSisBlacklisted___Group_C6_ATAC_minus_WT_ATACDEG.xls",
                       header = TRUE) %>%
  na.omit
C6overWT_peaks <- strsplit(C6overWT$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
C6overWT_gr <- GRanges(seqnames = C6overWT_peaks[,2], 
                       IRanges(start = as.numeric(C6overWT_peaks[,3]),
                               end = as.numeric(C6overWT_peaks[,4])))

# get GRanges for clone 6 UP peaks
C6overWT_UP_P05 <- C6overWT[C6overWT$padj < 0.05 & C6overWT$log2FoldChange > 0, ]
head(C6overWT_UP_P05)
C6overWT_UP_P05_peaks <- strsplit(C6overWT_UP_P05$IDs, split = ":")  %>%
  do.call(rbind.data.frame, .) 
C6overWT_UP_P05_gr <- GRanges(seqnames = C6overWT_UP_P05_peaks[,2], 
                              IRanges(start = as.numeric(C6overWT_UP_P05_peaks[,3]),
                                      end = as.numeric(C6overWT_UP_P05_peaks[,4])))

# make consensus set of peaks that go up
bothClones_consensus <- c(C4overWT_gr, C6overWT_gr) %>%
  reduce

bothClones_UP_P05_gr <- bothClones_consensus[bothClones_consensus %over% C4overWT_UP_P05_gr & bothClones_consensus %over% C6overWT_UP_P05_gr]
bothClones_UP_P05_gr
rtracklayer::export(bothClones_UP_P05_gr, "bothClones_C4andC6_UP_P05_ATAC.bed")


# get sequences under consensus UP peaks
atac_up_p05 <- rtracklayer::import("bothClones_C4andC6_UP_P05_ATAC.bed")
names(atac_up_p05) <- seq(length(atac_up_p05))
atac_up_p05_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, atac_up_p05)
writeXStringSet(atac_up_p05_seq, file = "atac_up_p05_wholePeak_seq.fa")

# run meme
# motif database from here: https://meme-suite.org/meme/doc/download.html on 20230106

# outputs from this were used for Figure 6E
system2(command = "meme-chip",
        args = c("-oc atac_up_p05_middle200bp",
                 "-time 240",
                 "-ccut 200", # this means it will only look at middle 200 bp (like homer)
                 "-dna -order 2 -minw 6 -maxw 15",
                 "-db ../motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
                 "-db ../motif_databases/MOUSE/uniprobe_mouse.meme",
                 "-db ../motif_databases/EUKARYOTE/jolma2013.meme",
                 "-meme-mod zoops -meme-nmotifs 3 -meme-searchsize 100000",
                 "-streme-pvt 0.05 -streme-totallength 4000000",
                 "-centrimo-score 5.0 -centrimo-ethresh 10.0",
                 "atac_up_p05_wholePeak_seq.fa"),
        stdout = paste0("meme_logs/motif_UP_P05_out_",format(Sys.time(), "%m.%d.%Y-%H.%M"),".txt"), 
        stderr = paste0("meme_logs/motif_UP_P05_err_",format(Sys.time(), "%m.%d.%Y-%H.%M"),".txt")
)



