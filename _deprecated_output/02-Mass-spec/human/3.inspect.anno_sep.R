setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human")
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()

fread_c("./total_psm.txt") -> total_psm
read.table("S3/uniprot.human.sep.id.txt") -> uniprot.human.sep.id
total_psm %>% filter(grepl("^sp",Protein)) %>% 
  filter(!grepl("sp",Mapped.Proteins)) -> anno_unique_psms_1
##### 更新了获取唯一比对的uniprot蛋白的计算方式，应该还需要包括哪些protein里面没有sp，但是mapped protein里面只有一个sp的谱图 #####
total_psm %>% filter(!grepl("^sp",Protein)) %>% 
  filter(str_count(Mapped.Proteins,fixed("sp"))==1) -> anno_unique_psms_2
pattern <- "sp[^,]+"
# 将Protein替换为Mapped.Proteins中的Protein
anno_unique_psms_2$Protein <- str_extract_all(anno_unique_psms_2$Mapped.Proteins,pattern) %>% unlist()
rbind(anno_unique_psms_1,anno_unique_psms_2) -> anno_unique_psms
anno_unique_psms %>% filter(Protein %in% uniprot.human.sep.id$V1) -> anno_sep_unique_psm
n_distinct(anno_sep_unique_psm$Protein)
anno_sep_unique_psm %>% fwrite_c("./S3/anno_sep.psm.txt")
anno_sep_unique_psm %>% distinct(Protein) %>% fwrite("./S3/anno_sep.id.txt",col.names = F)
get_ms_info(anno_sep_unique_psm) -> anno_sep_ms_info
mean(anno_sep_ms_info$Unique_psm_n)
mean(anno_sep_ms_info$Unique_peptide_n)
fwrite_c(anno_sep_ms_info,"S3/anno_sep_ms_info.txt")

###### 唯一比对加上非唯一比对 #####
pattern <- "sp[^,]+"
str_extract_all(total_psm$Protein,pattern) %>% unlist() %>% unique() -> c_1
str_extract_all(total_psm$Mapped.Proteins,pattern) %>% unlist() %>% unique() -> c_2
c(c_1,c_2) %>% unique -> uniprot_id
# 1592
sum(uniprot.human.sep.id$V1 %in% uniprot_id)

intersect(unique(anno_sep_unique_psm$Protein),uniprot.human.sep.id$V1) %>% length()

# 计算非特异性比对
total_psm %>% filter(grepl("^sp",Protein)) %>% 
  filter(grepl("sp",Mapped.Proteins)) -> non_unique_psm_1
total_psm %>% filter(!grepl("^sp",Protein)) %>% 
  filter(str_count(Mapped.Proteins,fixed("sp"))>1) -> non_unique_psm_2
get_sp_id <- function(psm){
  pattern <- "sp[^,]+"
  str_extract_all(psm$Protein,pattern) %>% unlist()  -> c_1
  str_extract_all(psm$Mapped.Proteins,pattern) %>% unlist() -> c_2
  c(c_1,c_2) -> uniprot_id
  return(uniprot_id)
}
get_sp_id(non_unique_psm_1) -> sp_id_1
get_sp_id(non_unique_psm_2) -> sp_id_2
data.frame(c(sp_id_1,sp_id_2)) -> sp_id_non_unique
colnames(sp_id_non_unique) <- "Protein"
sp_id_non_unique %>% count(Protein) -> sp_id_non_unique_psm_n 
sp_id_non_unique_psm_n %>% filter(Protein %in% uniprot.human.sep.id$V1) -> anno_sep_non_unique_psm

unique(c(anno_sep_non_unique_psm$Protein,anno_sep_unique_psm$Protein)) %>% length()
