# 整理合并关于每个转录本的注释信息，包括了每个转录本的ID，TSL，APPRIS，转录本的长度，基因ID，基因名，CDS的长度以及转录本的类型
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
fread("/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_ranscript_stable_ID_version_CDS_length.txt",data.table = F) -> trans_cds_l
fread("/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_Transcript_stable_ID_version_TSL_APPRIS_Transcript_length_Gene_stable_ID_version_Gene_name.txt",data.table = F) -> trans_mul_feature
trans_cds_l[is.na(trans_cds_l)] <- 0
distinct(trans_cds_l) -> trans_cds_l
merge(trans_mul_feature,trans_cds_l) -> trans_mul_feature_add
colnames(trans_mul_feature_add) <- c("Transcript_ID","TSL","APPRIS","Transcript_l","Gene_ID","Gene_name","CDS_l")
trans_mul_feature_add$TSL %<>% str_extract("tsl\\w+")
fread("/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.trans_id.trans_type.txt",
      header = F,data.table = F) -> trans_id.trans_type
colnames(trans_id.trans_type) <- c("Transcript_ID","Transcript_type")
merge(trans_mul_feature_add,trans_id.trans_type,all.x = T) -> trans_mul_feature_add_1
trans_mul_feature_add_1[is.na(trans_mul_feature_add_1)] <- "None"
fwrite(trans_mul_feature_add_1,file = "/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_Transcript_stable_ID_version_TSL_APPRIS_Transcript_length_Gene_stable_ID_version_Gene_name_CDS_length_Transcript_type.txt",sep = '\t')