#!/usr/bin/env Rscript
# 功能：从 featureCounts 输出计算 RPKM
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
args <- commandArgs(T)
fc_output_file <- args[1]
libsize_file <- args[2]
bam_lst_file <- args[3]
if(is.na(args[1])){
  fc_output_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/fc_output_rna_seq/rna_seq_combined.txt"
  libsize_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/fc_output_rna_seq/libsize.txt"
  bam_lst_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/total_rna_bam_lst.txt"
}
fread_c(fc_output_file)  -> counts
counts[, colSums(!is.na(counts)) > 0] -> counts_1
fread_c(libsize_file) -> libsize
colnames(libsize) <- c("Sample","libsize")
# 如果是双端，文库大小减半
fread_c(bam_lst_file) -> bam_lst
merge(libsize,bam_lst,by.x="Sample",by.y="BamPath") %>% mutate(libsize_1=case_when(
  SingleEnd==0~libsize/2,
  TRUE~libsize
)) -> tmp
tmp %>% dplyr::select(Sample,libsize_1) %>% rename(libsize=libsize_1) -> libsize

counts_1[,7:ncol(counts_1)] %>% colnames() -> col_order
libsize[match(col_order,libsize$Sample),] -> libsize_reorder
all(libsize_reorder$Sample==col_order)
# 计算RPKM
get_rpkm <- function(counts,libsize){
  rkm <- counts[,7:ncol(counts)]/counts$Length*1000
  rpkm <- data.frame(t(t(rkm) / libsize) * 10^6 )
  cbind(GeneID = counts_1$Geneid, rpkm) -> rpkm_1
  return(rpkm_1)
}
get_rpkm(counts_1,libsize_reorder$libsize) -> rpkm
colnames(rpkm) <- c("GeneID",col_order)

# 输出结果
output_file <- sub(".txt$", "_RPKM.txt", fc_output_file)
write.table(rpkm, 
            file = output_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

message(paste("RPKM计算结果已保存到:", output_file))

