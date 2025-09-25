# 设置工作目录
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401")
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
args <- commandArgs(T)
if(is.na(args[1])){
  sep_path = "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/new_sep_list.txt"
  output_file = "./output/sep_add_basic_ms_ribo_info_group_retained.txt"
  # 每个微蛋白被多少独特的psm以及肽段支持
  ms_info_path = "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/ms_info_updated.txt"
}else{
  sep_path <- args[1]
  output_file <- args[2]
  ms_info_path <- args[3]
}
output_path <- "output"
# 如果不存在就创建路径
create_path(output_path)
trans_based_sorfs_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/trans_based_sorfs.txt"
trans_mul_feature_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt"
gpe_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.10.gpe"
# 由于之前选择主要转录本的脚本存在错误，因此需要转换一下ORF_trans_id
correct_incorrect_map_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/correct_incorrect_map.rds"
sep_ribo_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_ribo_merge_call_orfs_20250338/merge/nonCano.sorf.filtered.add.orf_id_trans.txt"

# 【如果选择主要转录本的步骤正确，则此块代码忽略】【比较耗时】
correct_incorrect_map <- readRDS(correct_incorrect_map_path)
update_orf_id_trans <- function(df,correct_incorrect_map){
  merge(correct_incorrect_map,df,by.x="ORF_id_trans_incorrect",by.y="ORF_id_trans") -> df_new_orf_id_trans
  df_new_orf_id_trans %>% mutate(ORF_id_trans_incorrect=NULL) %>% rename(ORF_id_trans=ORF_id_trans_correct) -> df_new_orf_id_trans_1
  return(df_new_orf_id_trans_1)
}
fread(sep_path) -> all_orfs
update_orf_id_trans(all_orfs,correct_incorrect_map) -> all_orfs
# 提取染色体、起始位置、终止位置、链方向等信息
matches <- str_match(all_orfs$ORF_id_trans, "(ENS\\w+\\.\\d+)([+-])(chr\\w+):(\\d+)-(\\d+)")
all_orfs$Strand <- matches[, 3]
all_orfs$Chr <- matches[, 4]
all_orfs$Start <- matches[, 5] %>% as.numeric()
all_orfs$End <- matches[, 6] %>% as.numeric()

# 合并ORF_id_seq等列【相对比较耗时】
trans_based_sorfs <- fread_c(trans_based_sorfs_path)
setDT(all_orfs)
setDT(trans_based_sorfs)
setkey(all_orfs, ORF_id_trans)
setkey(trans_based_sorfs, ORF_id_trans)
all_orfs_1 <- trans_based_sorfs[, .(ORF_id_trans, ORF_id_seq, Seq, ENS_id, Scodon, Type)][all_orfs]
all_orfs_1$Length <- nchar(all_orfs_1$Seq)

# 加载基因名、转录本类型等信息并整合
trans_mul_feature <- fread(trans_mul_feature_path)
colnames(trans_mul_feature) <- c("Gene_ID", "Transcript_ID", "Transcript_type", "Gene_type", "Gene_name")
merge(all_orfs_1, trans_mul_feature[, c("Transcript_ID", "Gene_name", "Transcript_type", "Gene_type")],
    by.x = "ENS_id", by.y = "Transcript_ID") -> all_orfs_1
all_orfs_1 %>% as.data.frame() -> all_orfs_1

# 注释小肽类型
gpe <- fread(gpe_path, header = FALSE, data.table = FALSE)
gpe[, c("V1", "V6", "V7")] -> gpe_selected
colnames(gpe_selected) <- c("ENS_id", "CDS_start", "CDS_end")
merge(all_orfs_1, gpe_selected) -> all_orfs_2
get_orf_type_2(all_orfs_2) -> all_orfs_3
fwrite_c(all_orfs_3, o_f(output_path,"sep.add_anno.txt"))

# 添加质谱 (MS) 信息
fread_c(ms_info_path) -> ms_info
update_orf_id_trans(ms_info,correct_incorrect_map) -> ms_info
merge(all_orfs_3, ms_info) -> sep_add_ms_info

# 整合核糖体证据信息
fread_c(sep_ribo_path) -> sep_ribo
sep_add_ms_info$Ribo_evidence <- ifelse(sep_add_ms_info$ORF_id_trans %in% sep_ribo$ORF_id_trans, 1, 0)

# 输出最终结果
fwrite_c(sep_add_ms_info, output_file)



