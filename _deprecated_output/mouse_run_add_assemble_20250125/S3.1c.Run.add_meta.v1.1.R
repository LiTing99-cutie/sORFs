# 给过滤后的sORF增加元信息
# Usage: Rscript S3.1c.Run.add_meta.v1.1.R nonCano.sorf.filtered.txt nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.txt nonCano.sorf.filtered.add_meta.txt

args <- commandArgs(TRUE)
# 1.加载R包
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
# nonCano.sorf.filtered.txt
filtered_sorfs <- args[1]
# filtered_sorfs <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/Ribo_ORFs_add_assemble_20250125/filter/nonCano.sorf.filtered.txt"
# nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.txt
sorfs_all_m_s <- args[2]
# sorfs_all_m_s <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/Ribo_ORFs_add_assemble_20250125/filter/nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.filtered.txt"
# nonCano.sorf.filtered.add_meta.txt
output_file <- args[3]
# output_file <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/Ribo_ORFs_add_assemble_20250125/filter/nonCano.sorf.filtered.add_meta.txt"
trans_mul_feature_path <- "/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_Transcript_stable_ID_version_TSL_APPRIS_Transcript_length_Gene_stable_ID_version_Gene_name_CDS_length_Transcript_type.txt"

# 2.整合染色体，起始位置，终止位置，链，长度等信息
fread(filtered_sorfs,data.table = F) -> all_orfs
# 20250128修改
matches <- str_match(all_orfs$ORF_id_trans, "((?:[A-Z]+\\d+\\.\\d+)|(?:[A-Z]+\\.\\d+\\.\\d+))([+-])(chr\\w+):(\\d+)-(\\d+)")
all_orfs$Strand <- matches[,3]
all_orfs$Chr <- matches[,4]
all_orfs$Start <- matches[,5] %>% as.numeric()
all_orfs$End <- matches[,6] %>% as.numeric()
all_orfs$Length <- nchar(all_orfs$Seq)

# 3.整合鉴定到的样本数目和方法数目
fread(sorfs_all_m_s,data.table = F) %>% group_by(ORF_id_trans) %>% summarise(
  Method_n=n_distinct(Method),
  Sample_n=n_distinct(Sample)
) -> condition_n

merge(all_orfs,condition_n,by = "ORF_id_trans") -> all_orfs

## 4.整合基因名和转录本类型
# 20250128修改
all_orfs %>% filter(grepl("^ENS",ENS_id)) -> all_orfs_annotated
all_orfs %>% filter(grepl("^STRG",ENS_id)) -> all_orfs_assembled
trans_mul_feature <- fread(trans_mul_feature_path,data.table = F)
merge(all_orfs_annotated,trans_mul_feature[,c("Transcript_ID","Gene_name","Transcript_type")],by.x = "ENS_id",by.y="Transcript_ID") -> all_orfs_annotated_1
matches <- str_match(all_orfs_assembled$ORF_id_trans, "((?:[A-Z]+\\d+\\.\\d+)|(?:[A-Z]+\\.\\d+\\.\\d+))([+-])(chr\\w+):(\\d+)-(\\d+)")
all_orfs_assembled$ENS_id <- matches[,2]
all_orfs_assembled$Gene_name <- str_match(all_orfs_assembled$ENS_id,"[A-Z]+\\.\\d+")
all_orfs_assembled$Transcript_type <- "Newly_assembled"
col_names <- c("ORF_id_trans","ORF_id_seq","Seq","Length","Chr","Start","End","Strand","ENS_id","Gene_name","Transcript_type","Method_n","Sample_n")
rbind(all_orfs_annotated_1[,col_names],all_orfs_assembled[,col_names]) -> all_orfs

# 5.重新排列下顺序
fwrite(all_orfs,file = output_file,sep = '\t')