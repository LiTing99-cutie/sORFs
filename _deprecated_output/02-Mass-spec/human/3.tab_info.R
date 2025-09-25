# 需要给已经注释的uniprot蛋白添加以下的信息
# [1] "ENS_id"          "ORF_id_trans"    "ORF_id_seq"      "Seq"             "Scodon"          "Type"            "Strand"         
# [8] "Chr"             "Start"           "End"             "Length"          "Gene_name"       "Transcript_type" "Gene_type"      
# [15] "CDS_start"       "CDS_end"  "Unique_psm_n"	"Unique_peptide_n" "Ribo_evidence" "ENS_p_id"
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3")
trans_mul_feature_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt"
sep_ribo_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_ribo_merge_call_orfs_20250338/merge_include_cano/nonCano.sorf.filtered.add.orf_id_trans.txt"
ms_info_path <- "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/anno_sep_ms_info.txt"
fread_c("./uniprot.human.sep.15.gpe") -> uniprot.human.sep.15.gpe
head(uniprot.human.sep.15.gpe)
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
       Type="Canonical",
       ORF_type_1="Canonical",
       ORF_type_2="Canonical",
       ORF_type_3="Canonical") -> m_2
# 加载基因名、转录本类型等信息并整合
trans_mul_feature <- fread(trans_mul_feature_path)
colnames(trans_mul_feature) <- c("Gene_ID", "Transcript_ID", "Transcript_type", "Gene_type", "Gene_name")
add_info_1 <- function(df){
  merge(df, trans_mul_feature[, c("Transcript_ID", "Gene_name", "Transcript_type", "Gene_type")],
        by.x = "ENS_id", by.y = "Transcript_ID") -> df
  df %>% as.data.frame() -> df
  return(df)
}
add_info_1(m_2) -> m_3
# 添加质谱 (MS) 信息
fread_c(ms_info_path) -> ms_info
# 168/188在可以转换为ensembl的list里面
sum(ms_info$ORF_id_trans %in% m_3$Uniprot_id)
merge(m_3, ms_info,by.x = "Uniprot_id",by.y="ORF_id_trans",all.x = T) -> m_4
# 整合核糖体证据信息
fread_c(sep_ribo_path) -> sep_ribo
m_4$Ribo_evidence <- ifelse(m_4$Seq %in% sep_ribo$Seq, 1, 0)
m_4[is.na(m_4)] <- 0
colnames(m_4)
fwrite_c(m_4,"./uniprot.human.sep.tab_info.txt")
