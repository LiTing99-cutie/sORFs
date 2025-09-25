source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
trans_mul_feature_path <- "./ts.meta.txt"
ncbi_seq_path <- "./GCF_000001405.40_GRCh38.p14_protein.rmdup.seq"
uniprot_seq_path <- "./uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.rmdup.seq"
output_path <- "./output_20250407"
readRDS("all_orfs_tf.rds") -> all_orfs_tf
# 3. 过滤相同ORF，assign转录本
Sys.time()
cat("Filtering&Assignment start\n")
# ORF_id可能对应多个ENS_id
ORF_id_ENS_id <- unique(all_orfs_tf[, .(ORF_id_seq, ENS_id)])
ORF_id_ENS_id %>% filter(grepl("^ENS",ENS_id)) -> ORF_id_ENS_id_annotated
# 读取转录本特征信息并合并
trans_mul_feature_add <- fread(trans_mul_feature_path)
colnames(trans_mul_feature_add) <- c("Transcript_ID","Gene_ID","Gene_name","TSL","APPRIS","CDS_l")
ORF_id_ENS_id_merge_feature <- merge(
  ORF_id_ENS_id_annotated,
  trans_mul_feature_add[, .(Transcript_ID, TSL, APPRIS, CDS_l)],
  by.x = "ENS_id",
  by.y = "Transcript_ID"
)
# 设置TSL和APPRIS的等级排序
tsl_levels <- c("tsl1", "tsl2", "tsl3", "tsl4", "tsl5", "tslNA")
appris_levels <- c("principal_1", "principal_2", "principal_3", "principal_4", "principal_5", "alternative_1", "alternative_2", "None")
ORF_id_ENS_id_merge_feature[, TSL := factor(TSL, levels = tsl_levels)]
ORF_id_ENS_id_merge_feature[, APPRIS := factor(APPRIS, levels = appris_levels)]
# 按条件排序并筛选出证据最强的转录本
ORF_id_ENS_id_top_trans <- ORF_id_ENS_id_merge_feature[
  order(ORF_id_seq, APPRIS, TSL, -CDS_l),
  .SD[1],
  by = ORF_id_seq
]
# 得到基于转录本去重的所有sORFs
all_orfs_top_trans <- merge(
  all_orfs_tf,
  ORF_id_ENS_id_top_trans[, .(ENS_id, ORF_id_seq)],
  by = c("ENS_id", "ORF_id_seq")
)
saveRDS(all_orfs_top_trans,file = o("all_orfs_top_trans.rds"))
Sys.time()
cat("Filtering&Assignment done\n")
rm(all_orfs_top_trans)
readRDS(o("all_orfs_top_trans.rds")) -> all_orfs_top_trans
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
fwrite_c(trans_based_sorfs,o("trans_based_sorfs.txt"))
Sys.time()

all_orfs_top_trans %>% .[,c("ORF_id_seq","ORF_id_trans")] -> all_orfs_ORF_id_seq_trans_map
saveRDS(all_orfs_ORF_id_seq_trans_map,o("all_orfs_ORF_id_seq_trans_map.rds"))
