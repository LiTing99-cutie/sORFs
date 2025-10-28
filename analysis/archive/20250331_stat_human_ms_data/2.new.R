source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
lib_plot()

setwd("/rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data")
output_path <- "output/"
##### 1. Statistics #####
# 统计下总的二级谱图的数量
readRDS("./output/sample_metadata_ordered.rds") -> sample_metadata
## 3645648 一共是365万的总谱图数量
sum(sample_metadata$ms2_count)
res_path <- "/rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/output/db_search_20250326"
get_sample_ms_stat_from_path(res_path) -> sample_ms_stat
merge(sample_ms_stat,sample_metadata,by="Sample") %>% mutate(ID_rate=PSM_N/ms2_count) -> sample_ms_stat_add_id_rate
## 平均的谱图鉴定率为20.9%
mean(sample_ms_stat_add_id_rate$ID_rate)
fwrite_c(sample_ms_stat_add_id_rate,"output/sample_ms_stat_add_id_rate.txt")

get_total_psm(res_path) -> total_psm
total_psm %>% filter(grepl("ENS",Protein)) %>% filter(Is.Unique=="true") -> all_sample_sep_unique_psm_1

# 1070个unique谱图支持的蛋白质
length(unique(all_sample_sep_unique_psm_1$Protein))
# 1296个unique谱图支持的肽段
length(unique(all_sample_sep_unique_psm_1$Peptide))

all_sample_sep_unique_psm_1 %>% count(Protein) -> sep_unique_psm_count
table(sep_unique_psm_count$n)

all_sample_sep_unique_psm_1 %>% count(Protein,Peptide) %>% count(Protein) -> sep_unique_peptide_count
table(sep_unique_peptide_count$n)

merge(sep_unique_psm_count,sep_unique_peptide_count,by="Protein") %>% rename(Psm_count=n.x,Peptide_count=n.y) -> sep_unique_psm_peptide_count

total_psm %>% filter(grepl("ENS",Protein)) %>% filter(Is.Unique=="false") %>% filter(!grepl("sp",Mapped.Proteins)) -> all_sample_sep_group_psm
length(unique(all_sample_sep_group_psm$Protein))

fwrite_c(total_psm,path = "output/total_psm.txt")
fwrite_c(all_sample_sep_unique_psm_1,path = "output/all_sample_sep_unique_psm.txt")
##### 2. Plot #####
##### 2.1 Saturation Plot #####
set.seed(123)
sample_order <- sample(sample_metadata$Sample)
cumu_sep_plot_v1(all_sample_sep_unique_psm_1,sample_order) -> cumu_res_1
ggsave(cumu_res_1[[2]],filename="plot/cumu_sample_n.pdf",height = 5,width = 10)
ggsave(cumu_res_1[[3]],filename="plot/cumu_sample_name.pdf",height = 5,width = 10)
cumu_sep_peptide_plot(all_sample_sep_unique_psm_1,sample_order) -> cumu_res_2
ggsave(cumu_res_2[[2]],filename="plot/cumu_pep_sample_n.pdf",height = 5,width = 10)
ggsave(cumu_res_2[[3]],filename="plot/cumu_pep_sample_name.pdf",height = 5,width = 10)

##### 2.2 Overlap between different methods #####
output_path <- "plot/"
all_sample_sep_unique_psm_1 %>% distinct(Protein,Sample) -> all_sample_sep
all_sample_sep %>% merge(sample_metadata,by = "Sample") -> all_sample_sep_m_meta
##### 2.2.1 Compare enzyme #####
compare_enzyme <- function(method){
  all_sample_sep_m_meta %>% filter(Enrichment==method) %>% filter(Eyzyme!="Null") -> df
  split(df[, "Protein"], df$Eyzyme) -> split_df
  venn_plot_n(split_df) -> p
  path <- paste0("eyzyme_compare_",method,".pdf")
  ggsave(p,filename=o(path),width = 8,height = 5)
  return(p)
}
lapply(c("MWCO","PAGE"), compare_enzyme) -> l
##### 2.2.2 Compare enrichment method #####
compare_method <- function(enzyme){
  all_sample_sep_m_meta %>% filter(Eyzyme==enzyme) -> df
  split(df[, "Protein"], df$Enrichment) -> split_df
  venn_plot_n(split_df) -> p
  path <- paste0("method_compare_",enzyme,".pdf")
  ggsave(p,filename=o(path),width = 8,height = 5)
  return(p)
}
lapply(c("Trypsin","Trypsin_LysN","Trypsin_Chymotrypsin"), compare_method) -> l

##### 2.2.3 Compare different methods #####
split(all_sample_sep_m_meta[, "Protein"], all_sample_sep_m_meta$Replicate) -> split_df
venn_plot_n(split_df)->p
ggsave(p,filename=o("replicate_compare.pdf"),width = 8,height = 5)
