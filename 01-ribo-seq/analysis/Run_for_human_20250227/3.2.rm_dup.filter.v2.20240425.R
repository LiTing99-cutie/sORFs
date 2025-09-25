# 在/home/user/data3/lit/project/sORFs/01-ribo-seq/S3.1b.Run.Select_repre_trans.v1.1.R的基础上修改，去掉选择代表性的转录本这一步骤
# 1.加载R包
args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
# 输入文件所在路径nonCano.sorf.meta.merge.raw.3_ways.txt
all_orfs_path <- args[1]
# 输出文件夹
output_path <- args[2]
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
if(is.na(args[1])){
  all_orfs_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_ribo_merge_call_orfs_20250338/merge_include_cano/nonCano.sorf.meta.merge.raw.3_ways.txt"
  output_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_ribo_merge_call_orfs_20250338/merge_include_cano/"    
}
# uniprot序列
uniprot_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.rmdup.seq"
# ncbi序列
ncbi_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/hg38/GCF_000001405.40_GRCh38.p14_protein.rmdup.seq"

# 2. 导入所有方法鉴定出的sORFs
all_orfs <- fread(all_orfs_path,header = F,data.table = F)
colnames(all_orfs) <- c("ORF_id_trans","Seq","Method")
matches <- str_match(all_orfs$ORF_id_trans, "([A-Z]+\\d+\\.\\d+[a-zA-Z_]*)([+-])(chr\\w+):(\\d+)-(\\d+)")
# 添加三列，转录本ID，Location以及ORF_id
all_orfs$ENS_id <- matches[,2]
all_orfs$Location <- str_extract(all_orfs$ORF_id_trans,"[+-]chr\\w+:\\d+[+-]\\d+")
# ORF_id代表独特的小肽
all_orfs$ORF_id_seq <- paste0(all_orfs$Location,":",all_orfs$Seq)

# 3.得到去重的所有sORFs
all_orfs[,c("Seq","ORF_id_seq")] %>% distinct(ORF_id_seq,.keep_all = T) -> all_orfs_distinct
fwrite(all_orfs_distinct,file = o("nonCano.sorf.distinct.txt"),sep = '\t')

# 4.去掉和NCBI或者uniprot蛋白质序列相同的sORF
# ncbi_seq <- fread(ncbi_seq_path,header = F)
# uniprot_seq <- fread(uniprot_seq_path,header = F)
# all_orfs_distinct$Seq %in% c(ncbi_seq$V1,uniprot_seq$V1) -> idx
# all_orfs_distinct[!idx,] -> all_orfs_filtered
# ## 导出unique的sorf
# fwrite(all_orfs_filtered,file = o("nonCano.sorf.filtered.txt"),sep = '\t')
all_orfs_distinct -> all_orfs_filtered

# 5. 导入ORF_id_trans的对应关系
rds_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/all_orfs_ORF_id_seq_trans_map.rds"
readRDS(rds_path) -> all_orfs_ORF_id_seq_trans_map
merge(all_orfs_filtered,all_orfs_ORF_id_seq_trans_map,by="ORF_id_seq",all.x=T) -> all_orfs_filtered_add_orf_trans
all_orfs_filtered_add_orf_trans %>% na.omit() -> all_orfs_filtered_add_orf_trans
fwrite(all_orfs_filtered_add_orf_trans,file = o("nonCano.sorf.filtered.add.orf_id_trans.txt"),sep = '\t')
