# 为了计算疏水性和等电点，一开始只导出了序列
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
# 读入之前注释好的最终的基于质谱的sep的list
fread_c("./output/sep_add_basic_ms_ribo_info_group_retained.txt") -> ms_sep
fread("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/anno.detected.seq",header = F) -> uniprot_detected
fread("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/anno.undetected.seq",header = F) -> uniprot_undetected
fread("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/trans_based_sorfs.txt") -> all_putative_sep
all_putative_sep %>% filter(!Seq %in% ms_sep$Seq) -> all_putative_sep_undetected
set.seed(123)
sample(1:nrow(all_putative_sep_undetected),nrow(ms_sep)) ->  sample_rows
all_putative_sep_undetected[sample_rows,] -> all_putative_sep_undetected_sampled
fwrite_c(all_putative_sep_undetected_sampled[,"ORF_id_trans"],path = "./output/S7/all_putative_sep_undetected_sampled.txt")
# 加上一些元数据信息
system(paste(
  "Rscript /home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/Uni.1_2.merge.v1.20250424.R",
  "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S7/all_putative_sep_undetected_sampled.txt",
  "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S7/all_putative_sep_undetected_sampled.tab_info.txt",
  "FALSE"
))
rbind(data.frame(Type="new_sep_detected",Seq=ms_sep$Seq),
      data.frame(Type="new_sep_undetected_sampled",Seq=all_putative_sep_undetected_sampled),
      data.frame(Type="uniprot_sep_detected",Seq=uniprot_detected$V1),
      data.frame(Type="uniprot_sep_undetected",Seq=uniprot_undetected$V1)) -> list_for_comparison

fwrite_c(list_for_comparison,path = "./output/S7/list_for_comparison.txt")
