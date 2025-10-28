setwd("/rd1/user/lit/project/sORFs/analysis/20250404_tmp_from_nanjing")
source("~/bin/lit_utils.R")
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
lib_text()

fread_c("../20250331_stat_human_ms_data/output/total_psm.txt") -> total_psm
read.table("./uniprot.human.sep.id.txt") -> uniprot.human.sep.id
total_psm %>% filter(grepl("^sp",Protein)) %>% 
  filter(!grepl("sp",Mapped.Proteins)) -> anno_psms
anno_psms %>% filter(Protein %in% uniprot.human.sep.id$V1) -> anno_sep_unique_psm
n_distinct(anno_sep_unique_psm$Protein)
anno_sep_unique_psm %>% fwrite_c("./output/anno_sep.psm.txt")
anno_sep_unique_psm %>% distinct(Protein) %>% fwrite("./output/anno_sep.id.txt",col.names = F)
# total_psm %>% filter(grepl("ENS",Protein)) %>% 
#   filter(grepl("sp",Mapped.Proteins)) -> tmp
# tmp %>% filter(grepl(uniprot.human.sep.id$V1,Mapped.Proteins,fixed = T))
