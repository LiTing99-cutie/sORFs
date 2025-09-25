source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528")

fread_c("./output/S3/sep_info_20250605.txt") -> sep_info
mutate(sep_info,MS_detection=case_when(
  Unique_psm_n!=0  ~ "MS_unique_peptide",
  Unique_psm_n==0 & All_psm_n!=0 ~ "MS_nonUnique_peptide",
  TRUE ~ "MS_no_detection"
)) -> sep_info

mutate(sep_info,Cano_type=case_when(
  ORF_type_1=="Canonical"  ~ "Cano",
  TRUE ~ "Uncano"
)) -> sep_info

table(sep_info$MS_detection,sep_info$Cano_type)
table(sep_info$Cano_type)
##### 经典小肽卡不同的基因表达阈值，被质谱鉴定出来的比例 #####
sep_info %>% filter(ORF_type_1=="Canonical") -> cano_sep
cus <- function(df,cutoff){
  if(cutoff=="None"){
    df -> tmp_1
  }else{
    df %>% filter(A>cutoff) -> tmp_1
  }
  nrow(tmp_1) -> N
  tmp_1 %>% count(MS_detection) -> tmp_2
  tmp_2 %>% .[3,2] -> MS_unique_peptide_n
  tmp_2[2,2]+tmp_2[3,2] -> MS_all_peptide_n
  return(list(N,MS_unique_peptide_n/N,MS_all_peptide_n/N))
}
cus(cano_sep,"None")
cus(cano_sep,0)
cus(cano_sep,1)

# 查看MS检测到的等级对应到的uniprot蛋白的等级
readRDS("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/sep_uniprot_id_orf_id_type.rds") -> sep_uniprot_id_orf_id_type
merge(cano_sep,sep_uniprot_id_orf_id_type,by.x = "ORF_id_trans",by.y = "ORF_id_trans") -> cano_sep_add_uniprot_pe_type
table(cano_sep_add_uniprot_pe_type$MS_detection,cano_sep_add_uniprot_pe_type$Type.y)

# 只选择PE=1的蛋白质，一共1572个
cano_sep_add_uniprot_pe_type %>% filter(Type.y=="PE_1") -> cano_sep_pe_1

# 0.41;0.56
cus(cano_sep_pe_1,"None")
cus(cano_sep_pe_1,0)
# 0.65;0.83
# 表达量大于1且PE等级为1，也就是Experimental evidence at protein level的蛋白质中，有65%有特异性的肽段支持，有83%有肽段支持
cus(cano_sep_pe_1,1)
cus(cano_sep_pe_1,5)
cus(cano_sep_pe_1,10)
cus(cano_sep_pe_1,20)

##### 经典小肽卡不同的基因表达阈值和Ribo-seq表达阈值，被质谱鉴定出来的比例 #####
cus <- function(df,cutoff_1,cutoff_2){
  if(cutoff_1=="None"){
    df -> tmp_1
  }else{
    df %>% filter(A>cutoff_1) -> tmp_1
  }
  if(cutoff_2=="None"){
    tmp_1 -> tmp_2
  }else{
    tmp_1 %>% filter(Ribo_rpkm>cutoff_2) -> tmp_2
  }
  nrow(tmp_2) -> N
  tmp_2 %>% count(MS_detection) -> tmp_3
  tmp_3 %>% .[3,2] -> MS_unique_peptide_n
  tmp_3[3,2]+tmp_3[2,2] -> MS_all_peptide_n
  return(list(N,MS_unique_peptide_n/N,MS_all_peptide_n/N))
}

cus(cano_sep_pe_1,1,0)
cus(cano_sep_pe_1,0,1)
# 66%,84%
cus(cano_sep_pe_1,1,1)
cus(cano_sep_pe_1,1,5)
