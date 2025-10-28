source("~/bin/lit_utils.R")
# 整合所有的元数据，加上ms2 count的数据
output_path <- "/rd1/user/lit/project/sORFs/analysis/20250523_human_ms_run/stat_output/S1"
create_path(output_path)
# 函数
## 替换df的sample列中的-为_
reformat_s <- function(df){
  df$Sample %>% gsub("-","_",.) -> df$Sample
  return(df)
}
# 导入数据
fread("/rd1/user/lit/project/sORFs/raw_data/MS_metadata_20250326_new_batch_human.csv",data.table = F) -> sample_metadata
fread("/rd1/user/lit/project/sORFs/raw_data/MS_metadata_20250523_human.txt") -> sample_metadata_20250525
## 这里的ms2 count包括了小鼠的
fread("/rd1/user/lit/project/sORFs/output/MS/stat/sample_ms2_count_formal_hm.add_h.txt",data.table = F) -> sample_ms2_count
reformat_s(sample_ms2_count) -> sample_ms2_count
fread("/rd1/user/lit/project/sORFs/analysis/20250523_human_ms_run/stat_output/S0/ms2_count.txt",data.table = F) -> sample_ms2_count_250525
# sample_metadata
sample_metadata <- mutate(sample_metadata,V5=NULL,V6=NULL)
colnames(sample_metadata) <- c("Sample","Eyzyme","Period","Enrichment")
reformat_s(sample_metadata) -> sample_metadata
sample_metadata_20250525 <- mutate(sample_metadata_20250525,V5=NULL)
colnames(sample_metadata_20250525) <- c("Sample","Eyzyme","Period","Enrichment")
rbind(sample_metadata,sample_metadata_20250525) -> sample_metadata_human_all
table(sample_metadata_human_all$Eyzyme,sample_metadata_human_all$Enrichment)
# 添加replicate 列
sample_metadata_human_all$Sample %>% str_split("_") %>% do.call(rbind,.) %>% .[,2] -> replicates
sample_metadata_human_all$Replicate <- replicates
# ms2 count
colnames(sample_ms2_count_250525) <- colnames(sample_ms2_count)
rbind(sample_ms2_count,sample_ms2_count_250525) -> sample_ms2_count_human_all
# 合并
## 20250525目前有部分的样本还没有跑完，因此有部分样本的ms2count是缺失的
merge(sample_metadata_human_all,sample_ms2_count_human_all) -> sample_metadata_human_all_1
# 导出
sample_metadata_human_all_1 %>% arrange(Enrichment) -> sample_metadata_ordered
sample_metadata_ordered$Sample -> sample_order
fwrite_c(sample_metadata_ordered[,"Sample",drop=FALSE],o("sample_order.txt"))
fwrite_c(sample_metadata_ordered,o("sample_metadata_ordered.txt"))

##### plot #####
##### 样本分布 #####
library(ggplot2)
library(dplyr)
library(scales)

# 计算每个组合的频数
count_data <- sample_metadata_human_all %>% 
  group_by(Enrichment, Eyzyme) %>% 
  summarise(count = n(), .groups = 'drop')

# 创建堆叠柱状图
library(ggplot2)
library(dplyr)

# 计算频数数据（包含每个组合的计数和每列总数）
plot_data <- sample_metadata_human_all %>% 
  group_by(Enrichment, Eyzyme) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  group_by(Enrichment) %>% 
  mutate(total = sum(count))  # 计算每个Enrichment的总数

# 创建堆叠柱状图
p <- ggplot(plot_data, aes(x = Enrichment, y = count, fill = Eyzyme)) +
  geom_col(position = "stack", width = 0.7, alpha = 0.8) +
  
  # 添加各部分的频数标签（显示在每个色块中间）
  geom_text(aes(label = ifelse(count > 0, count, "")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  
  # 添加总频数标签（显示在柱子顶部）
  geom_text(aes(x = Enrichment, y = total, label = total),
            inherit.aes = FALSE,  # 不继承前面的美学映射
            vjust = -0.5,         # 将标签放在柱子上方
            size = 3.5) +
  
  # 设置颜色和主题
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(x = "Enrichment",
       y = "Sample number",
       fill = "Eyzyme") +
  
  # 调整y轴范围，给顶部标签留出空间
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 

# 显示图形
p+theme_3() -> p_1
p_1

# 保存图形（可选）
ggsave(p_1,filename=o("enzyme_enrichment_distribution.pdf"),
       width = 8, height = 6, dpi = 300, bg = "white")

##### ms2count #####
# 10552762 接近1kw的原始谱图数量
sum(sample_metadata_human_all_1$ms2_count)

library(ggplot2)
ggplot(sample_metadata_human_all_1, aes(x = Enrichment, y = ms2_count, fill = Eyzyme)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.7,outlier.shape = NA) +
  stat_summary(
    fun = mean, geom = "point", shape = 18, size = 3, color = "red",
    position = position_dodge(0.8), show.legend = FALSE
  ) +
  labs(x = "Enrichment",
       y = "MS2 count",
       fill = "Eyzyme") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") -> p
p+theme_3() -> p_2
ggsave(p_2,filename=o("enzyme_enrichment_ms2_count.pdf"),
       width = 8, height = 6, dpi = 300, bg = "white")