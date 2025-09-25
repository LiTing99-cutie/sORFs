source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
total_psm_path <- "./total_psm.txt"
read_psm(total_psm_path) -> total_psm
total_psm %>% filter(grepl("ENS",Protein)) %>% filter(Is.Unique=="true") -> sep_unique_psm
fwrite_c(sep_unique_psm,"./sep_unique_psm.txt")
distinct(sep_unique_psm,Protein) -> sep
colnames(sep) <- "ORF_id_trans"
fwrite_c(sep,"./sep.txt")
sep_unique_psm %>% count(Protein) -> sep_unique_psm_n
sep_unique_psm %>% distinct(Peptide,.keep_all = T) %>% count(Protein) -> sep_unique_pep_n
merge(sep_unique_psm_n,sep_unique_pep_n,by="Protein") %>% rename(Unique_psm_n=n.x,
                                                                 Unique_peptide_n=n.y,
                                                                 ORF_id_trans=Protein) -> ms_info
update_orf_id_trans <- function(df){
  correct_incorrect_map_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/correct_incorrect_map.rds"
  correct_incorrect_map <- readRDS(correct_incorrect_map_path)
  merge(correct_incorrect_map,df,by.x="ORF_id_trans_incorrect",by.y="ORF_id_trans") -> df_new_orf_id_trans
  df_new_orf_id_trans %>% mutate(ORF_id_trans_incorrect=NULL) %>% rename(ORF_id_trans=ORF_id_trans_correct) -> df_new_orf_id_trans_1
  return(df_new_orf_id_trans_1)
}
update_orf_id_trans(ms_info) -> ms_info
fwrite_c(ms_info,"./ms_info.txt")
