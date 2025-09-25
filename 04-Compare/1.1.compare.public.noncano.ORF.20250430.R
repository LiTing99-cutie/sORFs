setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human")
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()

#数字和之前文章搜集到的不太一样，因此这部分搜集结果可用性比较低
readxl::read_xlsx("./Public_SEP/human/organized/SEPs_with_MS_evidence.xlsx") -> public_noncano_sorf
table(public_noncano_sorf$paper)
table(public_noncano_sorf$`MS evidence support`)

# 比较id的时候需要注意是更新前的ID还是更新后的ID
# grep -f /home/user/data3/lit/project/sORFs/02-Mass-spec/human/new_sep_list.txt /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/group_specific_db/uniprot.contam.sorfs.trans.modiHeader.fa |wc -l
# 把肽段比对到我们的数据库中，看看是否有比对上的ORF
fread_c("/home/user/data3/lit/project/sORFs/04-Compare/output/S1/pub_in_house_db_map_perfect.txt") -> pub_in_house_map_perfect
pub_in_house_map_perfect %>% distinct(V1,V2) %>% group_by(V1) %>% 
  mutate(n=n()) %>% 
  filter(n==1) -> unique_map
# 426
n_distinct(unique_map$V1)
# 374
n_distinct(unique_map$V2)
filter(unique_map,!grepl("^sp",V2)) -> unique_map_to_sep
# 368
n_distinct(unique_map_to_sep$V1)
# 344
n_distinct(unique_map_to_sep$V2)
fread_c("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/new_sep_list.txt") -> ms_sep
# 和我们鉴定出来的ORF有2个overlap
intersect(ms_sep$ORF_id_trans,unique_map$V2)
filter(unique_map_to_sep,V2 %in% intersect(ms_sep$ORF_id_trans,unique_map$V2))
# # 即使是非unique的比对
# pub_in_house_map_perfect %>% distinct(V1,V2) %>% group_by(V1) %>% 
#   mutate(n=n()) %>% 
#   filter(n>1) -> nonunique_map
# intersect(ms_sep$ORF_id_trans,nonunique_map$V2)
