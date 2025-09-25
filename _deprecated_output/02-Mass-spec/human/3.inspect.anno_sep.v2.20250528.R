# 将uniprot protein id转换为ORF_id_trans
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3")
# 读入文件
trans_mul_feature_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt"
fread_c("./uniprot.human.sep.15.gpe") -> uniprot.human.sep.15.gpe
# uniprot id和seq
fread("./uniprot.human.sep.tab",sep='\t',header = F,data.table = F) -> uniprot_id_seq
colnames(uniprot_id_seq) <- c("Uniprot_id","Seq")
# uniprot id和转录本id
fread("./perfect_match_to_ensId.uniq.txt",sep='\t',header = F,data.table = F) -> uniprot_id_ensembl_t_id_l
colnames(uniprot_id_ensembl_t_id_l) <- c("Uniprot_id","ENS_id","Length",".",".")
merge(uniprot_id_seq,uniprot_id_ensembl_t_id_l[,c("Uniprot_id","ENS_id","Length")],by="Uniprot_id") %>% 
  merge(uniprot.human.sep.15.gpe,by.x="ENS_id",by.y="V1") -> m_1
m_1 %>% dplyr::rename(Chr=V2,Strand=V3,Start=V6,End=V7) -> m_2
m_2[,!grepl("^V",colnames(m_2))] -> m_2
mutate(m_2,
       ORF_id_trans=paste0(ENS_id,Strand,Chr,":",Start,"-",End),
       ORF_id_seq=paste0(Strand,Chr,":",Start,"-",End,":",Seq),
       Scodon="ATG",
       Type="Canonical") -> m_3
m_3 %>% select(ORF_id_trans,ORF_id_seq,Seq,Scodon,Type) -> m_4
# 导出必要的列
fwrite_c(m_4,"./cano_sep_orf_id_tab_info.txt")
fwrite_c(m_4[,"ORF_id_trans",drop=FALSE],"./cano_sep_orf_id.txt")
fwrite_c(m_3[,c("Uniprot_id","ORF_id_trans"),drop=FALSE],"./cano_sep_orf_id_uniprot_id.txt")


# 20250605获取uniprot蛋白的注释等级
fread_c("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/cano_sep_orf_id_uniprot_id.txt") -> map
read.table("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/uniprot.human.sep.id.type.txt") -> uniprot.human.sep.id.type
merge(map,uniprot.human.sep.id.type,by.x="Uniprot_id",by.y="V1") -> sep_uniprot_id_orf_id_type
colnames(sep_uniprot_id_orf_id_type) <- c(colnames(sep_uniprot_id_orf_id_type)[1:2],"Type") 
saveRDS(sep_uniprot_id_orf_id_type,"/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/sep_uniprot_id_orf_id_type.rds")