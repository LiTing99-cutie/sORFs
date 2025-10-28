args <- commandArgs(TRUE)
source("~/bin/lit_utils.R")
lib_text()
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")

stat_path <- "output/stat/default_trans/stat.rds"
sample_order_path <- "output/sample_order.rds"
output_path_1 <- "output/stat/default_trans/"
output_path_2 <- "output/visual/default_trans/"


readRDS(stat_path) -> l
readRDS(sample_order_path) -> sample_order
l$all_sample_all_level_sep %>% filter(Confidence_Level==1) -> df_level_1
df_level_1 %>% count(Protein) -> Protein_MS_sample_n

# 得到level 1的小肽的长度
l$all_sample_all_level_sep %>% filter(Confidence_Level==1) %>%
  distinct(Protein,Length) -> final_pro_l

# 大于一个样本支持
filter(Protein_MS_sample_n,n>1) -> Protein_MS_sample_over_1
# 大于一个谱图支持
merge(filter(Protein_MS_sample_n,n==1),
      filter(l$all_sample_all_level_sep,Confidence_Level==1),
      by="Protein") -> one_sample_supp_sep
filter(one_sample_supp_sep,`Total Spectral Count`>1) %>% .[,"Protein",drop=FALSE] -> sep_over_1_spec
mutate(sep_over_1_spec,n=1) -> sep_over_1_spec
# 得到最终的list
rbind(Protein_MS_sample_over_1,sep_over_1_spec) -> final_list

final_list %>% merge(final_pro_l) -> final_list_1
output_path <- output_path_1
fwrite_c(final_list_1,o("final_list.txt"))

# 绘制累积曲线
df_level_1[df_level_1$Protein %in% final_list$Protein,] -> df_level_1_filtered
cumu_sep_plot(df_level_1_filtered) -> l
l[["cumulative_data"]] -> df
output_path <- output_path_2
ggsave(l[["p_1"]],filename=o("sample_count_cumu_filtered.pdf"),width = 10,height = 5)
ggsave(l[["p_2"]],filename=o("sample_name_cumu_filtered.pdf"),width = 10,height = 5)
l[["p_1"]]
l[["p_2"]]