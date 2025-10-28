args <- commandArgs(TRUE)

source("~/bin/lit_utils.R")
lib_text()
lib_plot()
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
stat_res_path <- args[1]
output_path <- args[2]
sample_order_path <- args[3]
sample_metadata_path <- args[4]
if(is.na(args[1])){
  stat_res_path <- "./output/stat/"
  output_path <- "./output/visual/"
  sample_order_path <- "output/sample_order.rds"
  sample_metadata_path <- "output/sample_metadata.rds"
}
if(!dir.exists(output_path)){
  dir.create(output_path,recursive = T)
}
i <- function(file_name){
  paste0(stat_res_path,"/",file_name)
}
o <- function(file_name){
  paste0(output_path,"/",file_name)
}

# 1. 基本的统计数据
## 1.1 总计
fread_c(i("all_sep_sta.txt")) -> all_sep_sta
bar_plot_basic(all_sep_sta,"Feature","N",label = T)+ylim(0,max(all_sep_sta$N)) -> p
ggsave(p,filename=o("sep_type_n.pdf"),width = 10,height = 5)

## 1.2.1 分样本 小肽
# 作图时样本的排列顺序
readRDS(file = sample_order_path) -> sample_order
order_sample <- function(df){
  df$Sample <- factor(df$Sample,levels = sample_order)
  return(df)
}
fread_c(i("sample_sep_ms_sta.txt")) -> sample_sep_ms_sta
# 去掉两个来自于buffer的样本
sample_sep_ms_sta %>% filter(.$Sample %in% sample_order) -> sample_sep_ms_sta
sample_sep_ms_sta %>% .[,1:6] %>% melt(id.vars = "Sample",variable.name = "Type",value.name = "Freq") -> df
order_sample(df) %>% bar_plot_basic(x = "Sample",y = "Freq",fil_col = "Type") -> p
ggsave(p,filename=o("sep_type_n_per_sample.pdf"),width = 10,height = 5)

## 1.2.2 分样本 质谱
sample_sep_ms_sta %>%  .[,c(1,7:9)] %>% melt(id.vars = "Sample",variable.name = "Type",value.name = "Freq") -> df
order_sample(df) %>% bar_plot_basic(x = "Sample",y = "Freq",fil_col = "Type") -> p
ggsave(p,filename=o("ms_stat_n_per_sample.pdf"),width = 10,height = 5)
order_sample(sample_sep_ms_sta) %>% bar_plot_basic("Sample","ID_rate") -> p
ggsave(p,filename=o("ID_rate.pdf"),width = 10,height = 5)

# 2. 每个小肽被支持的程度
fread_c(i("all_sample_all_level_sep.txt")) -> df
df %>% filter(Confidence_Level==1) -> df_level_1
df_level_1 %>% count(Protein) -> Protein_MS_sample_n
bar_plot_v1(Protein_MS_sample_n,"n",label=T)+ylim(0,max(table(Protein_MS_sample_n$n))) -> p
ggsave(p,filename=o("sep_sample_n.pdf"),width = 10,height = 5)

cumu_sep_plot <- function(df){
  # 统计累积 unique 蛋白质数量
  # 输入一个蛋白质和样本对应的数据框
  df$Sample <- factor(df$Sample,levels = sample_order)
  arrange(df,Sample) -> df
  df_1 <- df %>% distinct(Protein,.keep_all = T)
  cumulative_data <- df_1 %>%
    group_by(Sample) %>%
    summarise(Proteins = list(unique(Protein))) %>%
    mutate(CumulativeUniqueProteins = cumsum(sapply(Proteins, length)),
           SampleCount = row_number())
  
  # 绘制累积频数曲线
  ggplot(cumulative_data, aes(x = SampleCount, y = CumulativeUniqueProteins)) +
    geom_line() +
    geom_point() +
    theme_3() -> p_1
  
  # 确保 Sample 顺序
  cumulative_data$Sample <- factor(cumulative_data$Sample, levels = unique(cumulative_data$Sample))
  # cumulative_data$Sample <- factor(cumulative_data$Sample, levels = sample_order)
  ggplot(cumulative_data, aes(x = Sample, y = CumulativeUniqueProteins,group=1)) +
    geom_line() +
    geom_point() +
    theme_3(rotate = T) -> p_2
  
  return(list(cumulative_data=cumulative_data,p_1=p_1,p_2=p_2))
}
cumu_sep_plot(df_level_1) -> l
ggsave(l[["p_1"]],filename=o("sample_count_cumu.pdf"),width = 10,height = 5)
ggsave(l[["p_2"]],filename=o("sample_name_cumu.pdf"),width = 10,height = 5)

# 3. 样本之间的重叠情况
sample_split_protein <- split(df_level_1[, "Protein"], df_level_1$Sample)
readRDS(sample_metadata_path) -> sample_metadata
# comparison enzyme
cus_1 <- function(enzyme){
  sample_lst <- filter(sample_metadata,Enrichment=="MWCO") %>% filter(Eyzyme==enzyme) %>% .$Experiment
  unlist(sample_split_protein[sample_lst])
}
eyzymes <- c("Trypsin", "Trypsin_ArgC", "Trypsin_AspN", "Trypsin_GluC", "Trypsin_LysN", "Trypsin_LysC", "Trypsin_Chymotrypsin")
lapply(eyzymes, cus_1) -> l
names(l) <- eyzymes
venn_plot_n(l) -> p
ggsave(p,filename=o("eyzyme_compare_mwco.pdf"),width = 8,height = 5)

# comparison enrichment method
cus_2 <- function(method){
  sample_lst <- filter(sample_metadata,Eyzyme=="Trypsin") %>% filter(Enrichment==method) %>% .$Experiment
  unlist(sample_split_protein[sample_lst])
}
methods <- c("MWCO","PAGE")
lapply(methods, cus_2) -> l
names(l) <- methods
venn_plot_n(l) -> p
ggsave(p,filename=o("method_compare_trypsin.pdf"),width = 8,height = 5)

# 比较LC_T这个酶切方法的是三个不同重复
sample_metadata$Experiment[grep("LC_T",sample_metadata$Experiment)] -> sample_lst
venn_plot_n(sample_split_protein[sample_lst])
