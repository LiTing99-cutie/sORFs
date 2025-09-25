args <- commandArgs(T)
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
# setwd("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/")
# read.table("./test_output_20250606/ribocode_offset_tab.txt") -> f_1
# read.table("./test_output_20250606/ribotish_offset_tab.txt") -> f_2
# Bam_name <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/01-output/call-orfs/p21_0523_1/output/alignment/p21_0523_1_Aligned.sortedByCoord.out.bam"
# output_file <- "./test_output_20250606/all_offset_tab.txt"
read.table(args[1]) -> f_1
read.table(args[2]) -> f_2
args[3] -> Bam_name
output_file <- args[4]
merge(f_1,f_2,by="V1",all = T) -> m
colnames(m) <- c("Length","Ribocode_offset","Ribotish_offset")
head(m)
# 如果一致，那么保留
# 如果不一致，保留被任意软件鉴定出来是12或者13的offset
# 如果都没有12或者13，那么就丢弃这个片段
mutate(m,offset_cura_1=case_when(
  Ribocode_offset == Ribotish_offset ~ Ribocode_offset,
  (Ribocode_offset!=12 | Ribocode_offset!=13 | is.na(Ribocode_offset)) & (Ribotish_offset==12|Ribotish_offset==13) ~ Ribotish_offset,
  (Ribotish_offset!=12 | Ribotish_offset!=13 | is.na(Ribotish_offset)) & (Ribocode_offset==12|Ribocode_offset==13) ~ Ribocode_offset,
)) -> m_1
mutate(m_1,Note=case_when(
  Ribocode_offset == Ribotish_offset ~ "Consistent",
  (Ribocode_offset!=12 | Ribocode_offset!=13 | is.na(Ribocode_offset)) & (Ribotish_offset==12|Ribotish_offset==13) ~ "Inconsistent",
  (Ribotish_offset!=12 | Ribotish_offset!=13 | is.na(Ribotish_offset)) & (Ribocode_offset==12|Ribocode_offset==13) ~ "Inconsistent",
  TRUE ~ "Inconsistent"
)) -> m_1
m_1$Bam_name <- Bam_name
fwrite(m_1,output_file,col.names = F,sep = '\t')
