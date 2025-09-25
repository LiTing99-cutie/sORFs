source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human")
total_psm_path <- "./total_psm.txt"
psm_group_retain_path <- "./psm_retained_group_based_on_transcript.txt"
sample_order_path <- "./sample_order.rds"
sample_metadata_path <- "./sample_metadata_ordered.txt"
output_path <- "S5/plot/"
create_path(output_path)
fread_c(total_psm_path) -> total_psm
total_psm %>% filter(grepl("ENS",Protein)) %>% filter(Is.Unique%in%c("true","TRUE","True")) -> psm_1
fread_c(psm_group_retain_path) -> psm
readRDS(sample_order_path) -> sample_order
fread_c(sample_metadata_path) -> sample_metadata
intersect(colnames(psm),colnames(psm_1)) -> com_cols
rbind(psm[,com_cols],psm_1[,com_cols]) -> psm_new_list
df <- psm_new_list
fwrite_c(psm_new_list,path = "./psm_new_list.txt")
##### 1. Saturation Plot #####
# 统计累积 unique 蛋白质数量
# 输入一个蛋白质和样本对应的数据框
set.seed(123)
seeds <- sample(1:33,33)
sample_order[seeds] -> sample_order_sampled
df$Sample <- factor(df$Sample,levels = sample_order_sampled)
arrange(df,Sample) -> df
df_1 <- df %>% distinct(Protein,.keep_all = T)
cumulative_data <- df_1 %>%
  group_by(Sample) %>%
  summarise(Proteins = list(unique(Protein))) %>%
  mutate(CumulativeUniqueProteins = cumsum(sapply(Proteins, length)),
         SampleCount = row_number())
df %>% distinct(Sample,Protein) %>% count(Sample) -> sample_unique_protein
merge(cumulative_data,sample_unique_protein,by="Sample") -> cumu_data_m_sample_unique
# 创建主图，柱状图和折线图共存
library(colorspace)
hcl_palettes(plot = TRUE)
ggplot(cumu_data_m_sample_unique, aes(x = Sample)) +
  geom_line(aes(y = CumulativeUniqueProteins, group = 1), size = 1,color=hcl.colors(3,"Dark 3")[1]) +
  # geom_area(aes(y = CumulativeUniqueProteins, group = 1),fill = "darkblue", alpha = 0.5)+
  geom_bar(aes(y = n), stat = "identity", fill = hcl.colors(3,"Dark 3")[3], width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_3(rotate = T) -> p
ggsave(p,filename=o("cumu_pep_sample_name_add_bar.pdf"),height = 5,width = 10)

ggplot(cumu_data_m_sample_unique, aes(x = SampleCount)) +
  geom_line(aes(y = CumulativeUniqueProteins, group = 1), size = 1,color=hcl.colors(3,"Dark 3")[1]) +
  # geom_area(aes(y = CumulativeUniqueProteins, group = 1),fill = "darkblue", alpha = 0.5)+
  geom_bar(aes(y = n), stat = "identity", fill = hcl.colors(3,"Dark 3")[3], width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_3(rotate = T) -> p
ggsave(p,filename=o("cumu_pep_sample_n_add_bar.pdf"),height = 5,width = 6)
##### 2. Overlap between different methods #####
psm_new_list %>% distinct(Protein,Sample) -> all_sample_sep
all_sample_sep %>% merge(sample_metadata,by = "Sample") -> all_sample_sep_m_meta
##### 2.1 Compare enzyme #####
compare_enzyme <- function(method){
  all_sample_sep_m_meta %>% filter(Enrichment==method) %>% filter(Eyzyme!="Null") -> df
  split(df[, "Protein"], df$Eyzyme) -> split_df
  venn_plot_n(split_df) -> p
  path <- paste0("eyzyme_compare_",method,".pdf")
  ggsave(p,filename=o(path),width = 8,height = 5)
  return(p)
}
lapply(c("MWCO","PAGE"), compare_enzyme) -> l
##### 2.2 Compare enrichment method #####
compare_method <- function(enzyme){
  all_sample_sep_m_meta %>% filter(Eyzyme==enzyme) -> df
  split(df[, "Protein"], df$Enrichment) -> split_df
  venn_plot_n(split_df) -> p
  path <- paste0("method_compare_",enzyme,".pdf")
  ggsave(p,filename=o(path),width = 8,height = 5)
  return(p)
}
lapply(c("Trypsin","Trypsin_LysN","Trypsin_Chymotrypsin"), compare_method) -> l

##### 2.3 Compare different replicate #####
split(all_sample_sep_m_meta[, "Protein"], all_sample_sep_m_meta$Replicate) -> split_df
venn_plot_n(split_df)->p
ggsave(p,filename=o("replicate_compare.pdf"),width = 8,height = 5)