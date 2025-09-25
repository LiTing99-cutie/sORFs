# 加载需要的R包
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
# 设定工作路径
work_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/trans_based_database_20250322"
if(!dir.exists(work_path)){
  dir.create(work_path,recursive = T)
}
setwd(work_path)
output_path <- "./output/"
if(!dir.exists(output_path)){
  dir.create(output_path,recursive = T)
}
# ORF序列
trans_orf_fa_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/mouse_run_add_assemble_20250125/annotation/RibORF_annot/candidateORF.prot.tab"
# ORF ID，5'端位置和3'端位置
trans_orf_gpe_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/mouse_run_add_assemble_20250125/annotation/RibORF_annot/candidateORF.genepred.txt"
# 转录本信息，通过/home/user/data3/lit/project/sORFs/01-ribo-seq/S3.1b1.Run.Organize_transcript_meta.R生成，如果是其他物种，同样需要生成
trans_mul_feature_path <- "/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_Transcript_stable_ID_version_TSL_APPRIS_Transcript_length_Gene_stable_ID_version_Gene_name_CDS_length_Transcript_type.txt"
# uniprot序列
uniprot_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform.rmdup.seq"
# ncbi序列
ncbi_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup.seq"
# 1. 加载所需数据
trans_orf_gpe <- fread(trans_orf_gpe_path, select = c("V1", "V6", "V7"))
setnames(trans_orf_gpe, c("V1", "V6", "V7"), c("RibORF_id", "codon5", "codon3"))
trans_orf_fa <- fread(trans_orf_fa_path, select = c(1, 2),sep = '\t')
setnames(trans_orf_fa, c("RibORF_id", "Seq"))
all_orfs <- merge(trans_orf_fa, trans_orf_gpe, by = "RibORF_id")
saveRDS(all_orfs,"./output/all_orfs.rds")
# 2. 统一ID名称
tf_riborf_id <- function(df) {
  # 拆分RibORF_id生成统一ID
  df[, c("ENS_id", "Chr", "Strand", "Rank", "Span", "RelaStart", "RelaEnd", "Type", "Scodon") :=
       tstrsplit(RibORF_id, split = "[|:]", fixed = FALSE)]
  df[, ORF_id_trans := paste0(ENS_id, Strand, Chr, ":", codon5, "-", codon3)]
  df[, Location := str_extract(ORF_id_trans,"[+-]chr\\w+:\\d+[+-]\\d+")]
  df[, ORF_id_seq := paste0(Location, ":", Seq)]
  return(df)
}
all_orfs_tf <- tf_riborf_id(all_orfs)
saveRDS(all_orfs_tf,"./output/all_orfs_tf.rds")
rm(list = ls())
readRDS("./output/all_orfs_tf.rds") -> all_orfs_tf
# 3. 过滤相同ORF，assign转录本
Sys.time()
cat("Filtering&Assignment start\n")
# ORF_id可能对应多个ENS_id
ORF_id_ENS_id <- unique(all_orfs_tf[, .(ORF_id_seq, ENS_id)])
ORF_id_ENS_id %>% filter(grepl("^ENS",ENS_id)) -> ORF_id_ENS_id_annotated
ORF_id_ENS_id %>% filter(grepl("^STRG",ENS_id)) -> ORF_id_ENS_id_assembled
# 读取转录本特征信息并合并
trans_mul_feature_add <- fread(trans_mul_feature_path)
ORF_id_ENS_id_merge_feature <- merge(
  ORF_id_ENS_id_annotated,
  trans_mul_feature_add[, .(Transcript_ID, TSL, APPRIS, CDS_l)],
  by.x = "ENS_id",
  by.y = "Transcript_ID"
)
# 设置TSL和APPRIS的等级排序
tsl_levels <- c("tsl1", "tsl2", "tsl3", "tsl4", "tsl5", "tslNA")
appris_levels <- c("principal1", "principal2", "principal3", "principal4", "principal5", "alternative1", "alternative2", "")
ORF_id_ENS_id_merge_feature[, TSL := factor(TSL, levels = tsl_levels)]
ORF_id_ENS_id_merge_feature[, APPRIS := factor(APPRIS, levels = appris_levels)]
# 按条件排序并筛选出证据最强的转录本
ORF_id_ENS_id_top_trans <- ORF_id_ENS_id_merge_feature[
  order(ORF_id_seq, APPRIS, TSL, -CDS_l),
  .SD[1],
  by = ORF_id_seq
]
# Assign新的转录本类型，如果组装后的转录本和注释的转录本编码同一个小肽，那么优先注释的转录本
ORF_id_ENS_id_assembled$trans_type_1 <- "Assembled"
ORF_id_ENS_id_top_trans$trans_type_1 <- "Annotated"
rbind(ORF_id_ENS_id_assembled[,c("ENS_id","ORF_id_seq","trans_type_1")],
      ORF_id_ENS_id_top_trans[,c("ENS_id","ORF_id_seq","trans_type_1")]) -> tmp
tmp$trans_type_1 <- factor(tmp$trans_type_1,levels = c("Annotated","Assembled"))
tmp %>% 
  group_by(ORF_id_seq) %>%
  arrange(trans_type_1) %>%  # 按TSL、APPRIS等级以及cds_length降序排序
  slice(1) %>%  # 保留每组的第一行
  ungroup() -> ORF_id_ENS_id_top_trans_1
data.table(ORF_id_ENS_id_top_trans_1) -> ORF_id_ENS_id_top_trans_2
# 得到基于转录本去重的所有sORFs
all_orfs_top_trans <- merge(
  all_orfs_tf,
  ORF_id_ENS_id_top_trans_2[, .(ENS_id, ORF_id_seq)],
  by = c("ENS_id", "ORF_id_seq")
)
saveRDS(all_orfs_top_trans,file = "all_orfs_top_trans.rds")
Sys.time()
cat("Filtering&Assignment done\n")
rm(list = ls())
readRDS("all_orfs_top_trans.rds") -> all_orfs_top_trans
# 4.去掉和NCBI或Uniprot蛋白质序列相同的sORF
ncbi_seq <- fread(ncbi_seq_path, header = FALSE)
uniprot_seq <- fread(uniprot_seq_path, header = FALSE)
all_orfs_top_trans_filter <- all_orfs_top_trans[!all_orfs_top_trans$Seq %in% c(ncbi_seq$V1, uniprot_seq$V1)]
# 5.过滤得到小肽
## 计算长度
Sys.time()
cat("Get sORFs start\n")
all_orfs_top_trans_filter[,Length := str_length(Seq)]
all_orfs_top_trans_filter[,c("ORF_id_trans","ORF_id_seq","Seq","Length","Chr","codon5","codon3","Strand","ENS_id","Type","Scodon")] -> all_orfs_top_trans_filter
all_orfs_top_trans_filter[Length >= 6 & Length <= 150] -> trans_based_sorfs
fwrite_c(trans_based_sorfs,paste0(output_path,"/","trans_based_sorfs.txt"))
Sys.time()
cat("Get sORFs done\n")
# 5.统计每个步骤过滤掉的orf数量
data.frame(all_orfs_top_trans_N=nrow(all_orfs_top_trans),
           all_orfs_top_trans_filter_N=nrow(all_orfs_top_trans_filter),
           trans_based_sorfs_N=nrow(trans_based_sorfs)) -> trans_based_sorfs_stat
fwrite_c(trans_based_sorfs_stat,paste0(output_path,"/","trans_based_sorfs_stat.txt"))
