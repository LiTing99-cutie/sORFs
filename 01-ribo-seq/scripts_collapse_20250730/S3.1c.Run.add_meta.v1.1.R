# 给过滤后的sORF增加元信息
# Usage: Rscript S3.1c.Run.add_meta.v1.1.R nonCano.sorf.filtered.txt nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.txt nonCano.sorf.filtered.add_meta.txt

args <- commandArgs(TRUE)
# 1.加载R包
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
# nonCano.sorf.filtered.txt
filtered_sorfs <- args[1]
# nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.txt
sorfs_all_m_s <- args[2]
# nonCano.sorf.filtered.add_meta.txt
output_file <- args[3]
trans_mul_feature_path <- "/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_Transcript_stable_ID_version_TSL_APPRIS_Transcript_length_Gene_stable_ID_version_Gene_name_CDS_length_Transcript_type.txt"

# 2.整合染色体，起始位置，终止位置，链，长度等信息
fread(filtered_sorfs) -> all_orfs
matches <- str_match(all_orfs$ORF_id_trans, "(ENS\\w+\\.\\d+)([+-])(chr\\w+):(\\d+)-(\\d+)")
all_orfs$Strand <- matches[,3]
all_orfs$Chr <- matches[,4]
all_orfs$Start <- matches[,5] %>% as.numeric()
all_orfs$End <- matches[,6] %>% as.numeric()
all_orfs$Length <- nchar(all_orfs$Seq)

# 3.整合鉴定到的样本数目和方法数目
fread(sorfs_all_m_s) %>% group_by(ORF_id_trans) %>% summarise(
  Method_n=n_distinct(Method),
  Sample_n=n_distinct(Sample)
) -> condition_n

merge(all_orfs,condition_n,by = "ORF_id_trans") -> all_orfs

## 4.整合基因名和转录本类型
trans_mul_feature <- fread(trans_mul_feature_path)
merge(all_orfs,trans_mul_feature[,c("Transcript_ID","Gene_name","Transcript_type")],by.x = "ENS_id",by.y="Transcript_ID") -> all_orfs

# 5.重新排列下顺序
all_orfs[,c(2,4,3,9,6,7,8,5,1,12,13,10,11)] -> all_orfs
fwrite(all_orfs,file = output_file,sep = '\t')