args <- commandArgs(TRUE)

source("~/bin/lit_utils.R")
lib_text()
lib_plot()
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")

# stat_res_path <- "./output/MS/stat/2024_12_13_default_loose_para"
stat_res_path <- args[1]
# output_path <- "./output/MS/visual/2024_12_13_default_loose_para"
output_path <- args[2]
sample_order_path <- "./output/MS/visual/sample_order.rds"
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
ggsave(p,filename=o("sep_type_n.pdf"),width = 5,height = 5)

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
ggsave(p,filename=o("sep_type_n_per_sample.pdf"),width = 5,height = 5)

## 1.2.2 分样本 质谱
sample_sep_ms_sta %>%  .[,c(1,7:9)] %>% melt(id.vars = "Sample",variable.name = "Type",value.name = "Freq") -> df
order_sample(df) %>% bar_plot_basic(x = "Sample",y = "Freq",fil_col = "Type") -> p
ggsave(p,filename=o("ms_stat_n_per_sample.pdf"),width = 5,height = 5)
order_sample(sample_sep_ms_sta) %>% bar_plot_basic("Sample","ID_rate") -> p
ggsave(p,filename=o("ID_rate.pdf"),width = 5,height = 5)

# 2. 样本之间的重叠情况
fread_c(i("all_sample_all_level_sep.txt")) -> df
df %>% filter(Confidence_Level==1) -> df_level_1
sample_split_protein <- split(df_level_1[, "Protein"], df_level_1$Sample)

# comparison 1
compare_lst <- c("m14_less3K","m14_3_10K","m14_10_30K")
l <- list(MWCO=unlist(sample_split_protein[compare_lst]),
            C8=unlist(sample_split_protein["m14_C8"]),
            PCP=unlist(sample_split_protein["m14_PCP"]))
venn_plot_n(l) -> p
ggsave(p,filename=o("m14_method_compare.pdf"),width = 5,height = 5)

# comparison 1-1
compare_lst <- c("m14_less3K","m14_3_10K","m14_10_30K")
venn_plot_n(sample_split_protein[compare_lst]) -> p
ggsave(p,filename=o("m14_MWCO.pdf"),width = 5,height = 5)

# comparison 2
l <- list(Native=unlist(sample_split_protein[c("E16_N_T1","E16_N_T2")]),
            Tristricine=unlist(sample_split_protein["E16_T_T"]))
venn_plot_n(l) -> p
ggsave(p,filename=o("E16_TN_compare.pdf"),width = 5,height = 5)
compare_lst <- c("E16_N_T1","E16_N_T2")
venn_plot_n(sample_split_protein[compare_lst]) -> p
ggsave(p,filename=o("E16_T_compare_replicate.pdf"),width = 5,height = 5)

# comparison 3
compare_lst <- sample_order[grepl("E16_T",sample_order)]
venn_plot_n(sample_split_protein[compare_lst]) -> p
ggsave(p,filename=o("E16_TN_compare.pdf"),width = 5,height = 5)

venn_plot_n(sample_split_protein[c("F2","F4")])+venn_plot_n(sample_split_protein[c("S2","S4")])+venn_plot_n(sample_split_protein[c("P2","P4")]) -> p
ggsave(p,filename=o("MWCO_compare_per_replicate.pdf"),width = 15,height = 5)
venn_plot_n(sample_split_protein[c("F2","F4","S2","S4","P2","P4")]) -> p
ggsave(p,filename=o("MWCO_compare.pdf"),width = 5,height = 5)

# 3. 每个小肽被支持的程度

df_level_1 %>% count(Protein) -> Protein_MS_sample_n
bar_plot_v1(Protein_MS_sample_n,"n",label=T)+ylim(0,max(table(Protein_MS_sample_n$n))) -> p
ggsave(p,filename=o("sep_sample_n.pdf"),width = 5,height = 5)

cumu_sep_plot <- function(df){
  # 统计累积 unique 蛋白质数量
  # 输入一个蛋白质和样本对应的数据框
  df <- df %>% distinct(Protein,.keep_all = T)
  cumulative_data <- df %>%
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
  ggplot(cumulative_data, aes(x = Sample, y = CumulativeUniqueProteins,group=1)) +
    geom_line() +
    geom_point() +
    theme_3(rotate = T) -> p_2
  
  return(list(p_1,p_2))
}
cumu_sep_plot(df_level_1) -> l
ggsave(l[[1]],filename=o("sample_count_cumu.pdf"),width = 5,height = 5)
ggsave(l[[2]],filename=o("sample_name_cumu.pdf"),width = 5,height = 5)
