setwd("/rd1/user/lit/project/sORFs/analysis/20250404_tmp_from_nanjing")
source("~/bin/lit_utils.R")
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
lib_text()
fread("out_list.txt",fill = TRUE,skip=2) -> peptide_match
fread_c("sep_add_overlap_pro_filter.txt") -> in_house_res

# 匹配到已注释蛋白的肽段
filter(peptide_match,grepl("^sp",V2)) -> match_to_annotated
match_to_annotated %>% distinct(V1) %>% mutate(V2="match to annotated protein") -> match_to_annotated_1
match_to_annotated_1
# 没有匹配，以及特异性匹配到小蛋白的肽段
fread("match_n.txt",header = FALSE) -> match_info
match_info %>% filter(V2=="has no match") -> no_match
match_info %>% filter(V2=="has 1 match") -> unique_match
unique_match %>% filter(!V1 %in% match_to_annotated_1$V1) -> unique_match_to_sorfs
rbind(no_match,unique_match_to_sorfs) %>% rbind(match_to_annotated_1) -> res
merge(match_info[,"V1",drop = FALSE],res,by="V1",all.x = T) -> res_1
res_1[is.na(res_1)] <- "has multiple match to sorfs"
# 整理所有的结果
table(res_1$V2) %>% as.data.frame() %>% fwrite_c("output/public_peptide_match_statistics.txt")
# 比较和in house的小蛋白的交集
intersect(in_house_res$ORF_id_trans,unique_match_to_sorfs$V2)
# 查看不同的研究中不同类型的数量
fread_c("/rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data/public_processed_data/2024_Wacholder_bioarxiv/SupplementalTable1_all_unannotated_peptides.csv") -> pub
filter(pub,is_HLA!="TRUE") -> pub_non_HLA
pub_non_HLA %>% distinct(study_pmid,unmodified_peptide) -> pub_peptide
merge(pub_peptide,res_1,by.x = "unmodified_peptide",by.y = "V1") %>% 
  count(study_pmid,V2) -> plot_input
plot_input$study_pmid %<>% as.character(plot_input$study_pmid)
colnames(plot_input) <- c("Study_pmid","Type","N")
table(plot_input$Type)
mutate(plot_input,Type=case_when(
  Type=="has 1 match"~ "has unique match to SEPs",
  Type=="has multiple match to SEPs"~ "has multiple match to sorfs",
  Type=="match to annotated protein"~"has match to annotated protein",
  TRUE~Type
)) -> plot_input
bar_plot_basic(plot_input,x = "Study_pmid",y = "N",fil_col = "Type",label = T)
bar_plot_basic_stack(plot_input,x = "Study_pmid",y = "N",fil_col = "Type")->p
p+theme(legend.position = "right")->p_1
ggsave(p_1,filename="output/public_peptide_map.pdf",height = 5,width = 8)