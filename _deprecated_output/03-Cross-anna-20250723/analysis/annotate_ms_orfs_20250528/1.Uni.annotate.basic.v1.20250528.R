# 设置工作目录
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
args <- commandArgs(T)
if(is.na(args[1])){
  # 需要有列名ORF_id_trans
  sep_path = "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S1/noncano.sep.ms.specific.nonspecific.lst"
  orf_id_trans_map <- T
  ORF_id_seq_from_trans_based_sorfs <- T
  output_path <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S1"
  output_file <- "uncano_sep_tab_info.txt"
  # sep_path = "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/cano_sep_orf_id.txt"
  # orf_id_trans_map <- F
  # ORF_id_seq_from_trans_based_sorfs <- F
  # output_path <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S1"
  # output_file <- "test.txt"
}else{
  sep_path <- args[1]
  orf_id_trans_map <- args[2]
  ORF_id_seq_from_trans_based_sorfs <- args[3]
  output_path <- args[4]
  output_file <- args[5]
}

create_path(output_path)
trans_mul_feature_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt"
gpe_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.10.gpe"
trans_based_sorfs_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/trans_based_sorfs.txt"
anno_sep_info_path <- "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/cano_sep_orf_id_tab_info.txt"
fread(sep_path) -> all_orfs
colnames(all_orfs) <- "ORF_id_trans"
library(tictoc)
# 3.4min
if(orf_id_trans_map){
  # 由于之前选择主要转录本的脚本存在错误，因此需要转换一下ORF_trans_id
  # 【比较耗时】
  tic("更新ORF ID映射")
  update_orf_id_trans_dt(all_orfs) -> all_orfs
  toc()
}
# 提取染色体、起始位置、终止位置、链方向等信息
matches <- str_match(all_orfs$ORF_id_trans, "([A-Z]+\\d+\\.\\d+[a-zA-Z_]*)([+-])(chr\\w+):(\\d+)-(\\d+)")
all_orfs$ENS_id <-  matches[, 2]
all_orfs$Strand <- matches[, 3]
all_orfs$Chr <- matches[, 4]
all_orfs$Start <- matches[, 5] %>% as.numeric()
all_orfs$End <- matches[, 6] %>% as.numeric()

# 合并ORF_id_seq等列【相对比较耗时】
# 311.263s 5.2min
if(ORF_id_seq_from_trans_based_sorfs){
  tic("合并信息")
  trans_based_sorfs <- fread_c(trans_based_sorfs_path)
  setDT(all_orfs)
  setDT(trans_based_sorfs)
  setkey(all_orfs, ORF_id_trans)
  setkey(trans_based_sorfs, ORF_id_trans)
  all_orfs_1 <- trans_based_sorfs[, .(ORF_id_trans, ORF_id_seq, Seq, Scodon, Type)][all_orfs]
  all_orfs_1$Length <- nchar(all_orfs_1$Seq)
  toc()
}else{
  fread_c(anno_sep_info_path) -> anno_sep_info
  merge(all_orfs,anno_sep_info) -> all_orfs_1
}

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
fwrite_c(all_orfs_3, o(output_file))



