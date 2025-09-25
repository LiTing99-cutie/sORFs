# 对于一个ORF_id，如果多个转录本编码这个ORF，那么选择一个转录证据最强的转录本

# 1.加载R包
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
trans_mul_feature_path <- "/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_Transcript_stable_ID_version_TSL_APPRIS_Transcript_length_Gene_stable_ID_version_Gene_name_CDS_length_Transcript_type.txt"
all_orfs_path <- "./mouse_brain_output_20241011/Ribo_ORFs_merge/nonCano.sorf.meta.merge.raw.3_ways.all_samples.txt"
script_path <- "./S3.1b1.Run.Organize_transcript_meta.R"
uniprot_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform.rmdup.seq"
ncbi_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup.seq"

# 2. 导入所有方法所有样本鉴定出的sORFs
all_orfs <- fread(all_orfs_path,header = F,data.table = F)
colnames(all_orfs) <- c("ORF_id_trans","Seq","Method","Sample")
matches <- str_match(all_orfs$ORF_id_trans, "(ENSMUST\\d+\\.\\d+)([+-])(chr\\w+):(\\d+)-(\\d+)")
# 添加三列，转录本ID，Location以及ORF_id
all_orfs$ENS_id <- matches[,2]
all_orfs$Location <- str_extract(all_orfs$ORF_id_trans,"[+-]chr\\w+:\\d+[+-]\\d+")
# ORF_id代表独特的小肽
all_orfs$ORF_id_seq <- paste0(all_orfs$Location,":",all_orfs$Seq)
# ORF_id可能对应多个ENS_id
distinct(all_orfs,ORF_id_seq,ENS_id) -> ORF_id_ENS_id

# 3.整合从biomart中下载的转录本的信息并与sORF合并
source(script_path)
trans_mul_feature_add <- fread(trans_mul_feature_path)
merge(ORF_id_ENS_id,trans_mul_feature_add[,c("Transcript_ID","TSL","APPRIS","CDS_l")],by.x = "ENS_id",by.y = "Transcript_ID") -> ORF_id_ENS_id_merge_feature

# 4.筛选出证据最强的转录本
# 设置TSL和APPRIS的等级排序
tsl_levels <- c("tsl1", "tsl2", "tsl3", "tsl4", "tsl5", "tslNA")
appris_levels <- c("principal1", "principal2", "principal3","principal4", "principal5", "alternative1", "alternative2", "")
# 转换TSL和APPRIS为因子并指定等级
ORF_id_ENS_id_merge_feature <- ORF_id_ENS_id_merge_feature %>%
  mutate(TSL = factor(TSL, levels = tsl_levels),
         APPRIS = factor(APPRIS, levels = appris_levels))
# 进行排序并筛选出证据最强的转录本
ORF_id_ENS_id_top_trans <- ORF_id_ENS_id_merge_feature %>%
  group_by(ORF_id_seq) %>%
  arrange(APPRIS,TSL,desc(CDS_l)) %>%  # 按TSL、APPRIS等级以及cds_length降序排序
  slice(1) %>%  # 保留每组的第一行
  ungroup()     # 解除分组

# 5.得到基于转录本去重的所有sORFs
merge(all_orfs,ORF_id_ENS_id_top_trans[,c("ENS_id","ORF_id_seq")],c("ENS_id","ORF_id_seq")) -> all_orfs_top_trans
fwrite(all_orfs_top_trans,file = "./mouse_brain_output_20241011/filter/nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.txt",sep = '\t')
all_orfs_top_trans[,c("ORF_id_trans","Seq","ENS_id","ORF_id_seq")] %>% distinct(ORF_id_trans,.keep_all = T) -> all_orfs_distinct
fwrite(all_orfs_distinct,file = "./mouse_brain_output_20241011/filter/nonCano.sorf.distinct.txt",sep = '\t')

# 6.去掉和NCBI或者uniprot蛋白质序列相同的sORF
ncbi_seq <- fread(ncbi_seq_path,header = F)
uniprot_seq <- fread(uniprot_seq_path,header = F)
all_orfs_distinct$Seq %in% c(ncbi_seq$V1,uniprot_seq$V1) -> idx
all_orfs_distinct[!idx,] -> all_orfs_filtered
## 导出unique的sorf
fwrite(all_orfs_filtered,file = "./mouse_brain_output_20241011/filter/nonCano.sorf.filtered.txt",sep = '\t')
## 导出原始的sorf
all_orfs_top_trans$ORF_id_trans %in% all_orfs_filtered$ORF_id_trans -> idx
all_orfs_top_trans[idx,] -> all_orfs_top_trans_filtered
fwrite(all_orfs_top_trans_filtered,file = "./mouse_brain_output_20241011/filter/nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.filtered.txt",sep = '\t')
