source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
# setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401")
args <- commandArgs(TRUE)
fc_output_path <- args[1]
if(is.na(args[1])){
  fc_output_path <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S10/"
}
# 1. 计算表达量
add_rpkm <- function(counts,length){
  RPM <- (counts/sum(counts))*10^6
  RPKM <- (RPM/length)*10^3
  return(RPKM)
}
add_tpm <- function(counts, length) {
  # 1. 计算每百万转录本中每千碱基的reads数 (RPK)
  rpk <- counts / (length / 1000)  # 注意长度单位是kb
  # 2. 计算所有基因RPK的总和（按百万归一化）
  per_million_scale_factor <- sum(rpk) / 1e6
  # 3. 计算TPM
  tpm <- rpk / per_million_scale_factor
  return(tpm)
}

fread_c(paste0(fc_output_path,"rna_counts_1.txt"))  -> rna_counts_1
fread_c(paste0(fc_output_path,"rna_counts_2.txt")) -> rna_counts_2
fread_c(paste0(fc_output_path,"rna_counts_3.txt"))  -> rna_counts_3
fread_c(paste0(fc_output_path,"rna_counts_4.txt"))  -> rna_counts_4

merge(rna_counts_1,rna_counts_2,by="Geneid") %>% merge(rna_counts_3,by="Geneid") %>% 
  merge(rna_counts_4,by="Geneid") -> merged_rna_counts
merged_rna_counts %>% select(1,7,13,19,25) -> merged_rna_counts
merged_rna_counts$Mean_counts <- apply(merged_rna_counts[,2:5],1,mean)
add_rpkm(merged_rna_counts$Mean_counts,rna_counts_1$Length) -> rpkm
add_tpm(merged_rna_counts$Mean_counts,rna_counts_1$Length) -> tpm
merged_rna_counts$Mean_rpkm <- rpkm
merged_rna_counts$Mean_tpm <- tpm
fc_output_path -> output_path
fwrite_c(merged_rna_counts,o("merged_rna_counts.txt"))

