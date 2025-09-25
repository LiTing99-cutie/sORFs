
setwd("/rd1/user/lit/project/sORFs/analysis/20250404_tmp_from_nanjing")
source("~/bin/lit_utils.R")
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
lib_text()
fread("tmp_match/out_list.txt",fill = TRUE,skip=2) -> peptide_match
fread("tmp_match/match_n.txt",header = FALSE) -> match_info
get_type <- function(peptide_match,match_info){
  # 匹配到已注释蛋白的肽段
  filter(peptide_match,grepl("^sp",V2)) -> match_to_annotated
  match_to_annotated %>% distinct(V1) %>% mutate(V2="has match to annotated protein") -> match_to_annotated_1
  # 没有匹配，以及特异性匹配到小蛋白的肽段
  match_info %>% filter(V2=="has no match") -> no_match
  match_info %>% filter(V2=="has 1 match") -> unique_match
  unique_match %>% filter(!V1 %in% match_to_annotated_1$V1) -> unique_match_to_sorfs
  rbind(no_match,unique_match_to_sorfs) %>% rbind(match_to_annotated_1) -> res
  merge(match_info[,"V1",drop = FALSE],res,by="V1",all.x = T) -> res_1
  res_1[is.na(res_1)] <- "has multiple matches to SEPs"
  return(res_1)
}
get_type(peptide_match,match_info) -> df_in_house
df_in_house$study_pmid <- "In_house"
fread("out_list.txt",fill = TRUE,skip=2) -> peptide_match_1
fread("match_n.txt",header = FALSE) -> match_info_1
get_type(peptide_match_1,match_info_1) -> df_public
fread_c("/rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data/public_processed_data/2024_Wacholder_bioarxiv/SupplementalTable1_all_unannotated_peptides.csv") -> pub
filter(pub,is_HLA!="TRUE") -> pub_non_HLA
pub_non_HLA %>% distinct(study_pmid,unmodified_peptide) -> pub_peptide
merge(pub_peptide,df_public,by.x = "unmodified_peptide",by.y = "V1") -> df_public_1

df_in_house %>% select(1,3,2) -> df_in_house_1
colnames(df_in_house_1) <- colnames(df_public_1)

rbind(df_in_house_1,df_public_1) %>% count(study_pmid,V2) %>% as.data.frame()-> plot_input
colnames(plot_input) <- c("Study_pmid","Type","N")
mutate(plot_input,Type=case_when(
  Type=="has 1 match"~ "has unique match to SEPs",
  Type=="match to annotated protein"~"has match to annotated protein",
  TRUE~Type
)) -> plot_input
bar_plot_basic(plot_input,x = "Study_pmid",y = "N",fil_col = "Type",label = T)
bar_plot_basic_stack(plot_input,x = "Study_pmid",y = "N",fil_col = "Type")->p
p+theme(legend.position = "right")->p_1
ggsave(p_1,filename="output/public_peptide_map.pdf",height = 5,width = 8)

