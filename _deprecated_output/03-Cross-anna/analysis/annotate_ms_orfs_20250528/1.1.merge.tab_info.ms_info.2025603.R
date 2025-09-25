source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528")
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S1/cano_sep_tab_info.txt") -> cano_sep_tab_info
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S1/uncano_sep_tab_info.txt") -> uncano_sep_tab_info
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S1/psm_sep_all.txt") -> psm_sep_all
fread_c("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/cano_sep_orf_id_uniprot_id.txt") -> map

cano_sep_tab_info$Length <- nchar(cano_sep_tab_info$Seq)
uncano_sep_tab_info %>% select(-i.ORF_id_seq) -> uncano_sep_tab_info
rbind(cano_sep_tab_info,uncano_sep_tab_info) -> sep_tab_info
table(sep_tab_info$ORF_type_1)

# 合并cano以及uncano的ms info
psm_sep_all %>% filter(grepl("^ENS",V35)) -> tmp_1
psm_sep_all %>% filter(grepl("^sp",V35)) -> tmp_2
tmp_2[,c(1:42,44,43)] -> tmp_3
colnames(tmp_3) <- colnames(tmp_1)
rbind(tmp_1,tmp_3) -> psm_sep_all_1
# 计算psm_n以及peptide_n
cus <- function(df){
  count(df,V35) -> pro_psm_n
  count(df,V35,V3) %>% count(V35) -> pro_peptide_n
  merge(pro_psm_n,pro_peptide_n,by="V35") -> m
  return(m)
}
cus(psm_sep_all_1) %>% rename(ORF_id_trans=V35,All_psm_n=n.x,All_peptide_n=n.y) -> df_1
cus(psm_sep_all_1 %>% filter(V43=="Unique")) %>% rename(ORF_id_trans=V35,Unique_psm_n=n.x,Unique_peptide_n=n.y) -> df_2
merge(df_1,df_2,by='ORF_id_trans',all.x = T) -> ms_info
ms_info[is.na(ms_info)] <- 0
# 更改ms info中的uniprot id为orf_id_trans
merge(ms_info,map,by.x="ORF_id_trans",by.y = "Uniprot_id",all.y = T) %>% mutate(ORF_id_trans=NULL) %>% 
  rename(ORF_id_trans=ORF_id_trans.y) -> cano_ms_info
cano_ms_info[is.na(cano_ms_info)] <- 0
ms_info %>% filter(grepl("^ENS",ORF_id_trans)) -> uncano_ms_info
update_orf_id_trans_dt(uncano_ms_info) -> uncano_ms_info_up
uncano_ms_info_up %>% as.data.frame() -> uncano_ms_info_up
rbind(uncano_ms_info_up[,colnames(uncano_ms_info)],cano_ms_info[,colnames(uncano_ms_info)]) -> ms_info_1

# 与sep_tab_info合并
merge(sep_tab_info,ms_info_1,by='ORF_id_trans') -> sep_info
sep_info %>% count(ORF_type_1)
sep_info %>% filter(Unique_psm_n>0) %>% count(ORF_type_1)
sep_info %>% filter(All_psm_n>0) %>% count(ORF_type_1)

# 导出
fwrite_c(sep_info,"./output/S1/sep_info_20250603.txt")
