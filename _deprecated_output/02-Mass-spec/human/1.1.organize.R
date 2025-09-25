setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/")
read_psm("./sep_unique_psm.txt") -> old_psm
read_psm("./psm_retained_group_based_on_transcript.txt") -> new_psm
rbind(old_psm,new_psm[,colnames(old_psm)]) -> updated_psm

get_ms_info <- function(psm){
  psm %>% count(Protein) -> sep_unique_psm_n
  psm %>% distinct(Peptide,.keep_all = T) %>% count(Protein) -> sep_unique_pep_n
  merge(sep_unique_psm_n,sep_unique_pep_n,by="Protein") %>% rename(Unique_psm_n=n.x,
                                                                   Unique_peptide_n=n.y,
                                                                   ORF_id_trans=Protein) -> ms_info
  return(ms_info)
}

get_ms_info(updated_psm) -> ms_info

fwrite_c(ms_info,"./ms_info_updated.txt")
