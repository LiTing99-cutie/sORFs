source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528")
##### in house RNA-seq数据整理 #####
read.table("/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/featureCounts/bam.lst") -> in_house_bam_lst
in_house_bam_lst$libType <- 2
in_house_bam_lst$SingleEnd <- 0
colnames(in_house_bam_lst) <- c("BamPath","LibraryType","SingleEnd")
##### public RNA-seq 数据整理 #####
read.table("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S5/enhanced_results_with_path.tsv",header = T) -> public_bam_lst
public_bam_lst$SingleEnd <- 1
mutate(public_bam_lst,SingleEnd=case_when(
  grepl("SRR.*_1",BamPath) ~ 0,
  TRUE ~ SingleEnd
)) -> public_bam_lst
##### 合并 #####
rbind(in_house_bam_lst,public_bam_lst[,colnames(in_house_bam_lst)]) -> total_bam_lst
total_bam_lst -> total_rna_bam_lst

fwrite_c(total_rna_bam_lst,"./output/S2/total_rna_bam_lst.txt")

##### in house Ribo-seq数据整理 #####
read.table("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/merge_in_house_bam_20250530/in_house_bam_org.lst") -> in_house_bam_lst
in_house_bam_lst$libType <- 1
colnames(in_house_bam_lst) <- c("BamPath","LibraryType")
##### public Ribo-seq 数据整理 #####
read.table("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/Ribo-seq/enhanced_results_with_path.tsv",header = T) -> public_bam_lst
table(public_bam_lst$LibraryType)
# 手动修改
public_bam_lst$LibraryType <- 1
##### 合并 #####
rbind(in_house_bam_lst,public_bam_lst[,colnames(in_house_bam_lst)]) -> total_bam_lst
total_bam_lst -> total_ribo_bam_lst
fwrite_c(total_ribo_bam_lst,"./output/S2/total_ribo_bam_lst.txt")
total_ribo_bam_lst$SingleEnd <- 1
fwrite_c(total_ribo_bam_lst,"./output/S2/total_ribo_bam_lst.txt")