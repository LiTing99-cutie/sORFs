process_1 <- function(df){
  high_quality_sample <- fread("./human_brain_output_20250227/high_quality_run.add.rm.20250328.txt",header = F) %>% .$V1
  high_quality_sample <- gsub("_GSM.*RNA-Seq", "", high_quality_sample)
  df$Sample <- gsub("_GSM.*RNA-Seq", "", df$Sample)
  df %>% arrange(Sample) -> df
  df$SampleGroup <- ifelse(df$Sample %in% high_quality_sample, "High_q", "Low_q")
  return(df)
}

##### uniq比对的reads比例和数目 #####

fread("./human_brain_output_20250227/stat/Merged_Uniquely_mapped_reads.ribo.txt") -> df
process_1(df) -> df  
# 将 uniquely_mapped_reads_rate 转换为数值形式
df$uniquely_mapped_reads_rate <- as.numeric(gsub("%", "", df$uniquely_mapped_reads_rate)) / 100
# 绘制 uniquely_mapped_reads_rate 的柱状图
p1 <- ggplot(df, aes(x = uniquely_mapped_reads_rate, y = Sample)) +
  geom_col() +
  labs(x = "Uniquely Mapped Reads Rate", y = "Sample") +
  theme_3() +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous( expand = c(0,0))
#  facet_wrap(~ SampleGroup, scales = "free_y")

# 绘制 Uniquely_mapped_reads_number 的柱状图
p2 <- ggplot(df, aes(x = Uniquely_mapped_reads_number, y = Sample)) +
  geom_col() +
  labs(x = "Uniquely Mapped Reads Number", y = "Sample") +
  theme_3() +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous( expand = c(0,0))
#  facet_wrap(~ SampleGroup, scales = "free_y")

# 显示两个图
p1
ggsave("./plot/mapped_rate.pdf",height = 12,width = 8)
p2
ggsave("./plot/mapped_number.pdf",height = 12,width = 8)


##### 不同区域的reads比例 #####
# 读取数据
df <- fread_c("human_brain_output_20250227/reads_distri_20250328/summary_tag_count.tsv")
process_1(df) -> df 
# 过滤只保留 CDS 和 UTR
df_filtered <- df %>%
  filter(Type %in% c("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons"))
# 计算每个样本总Tag_count
df_fraction <- df_filtered %>%
  group_by(Sample) %>%
  mutate(Total = sum(Tag_count),
         Fraction = Tag_count / Total)
# 绘制堆叠柱状图（水平）
ggplot(df_fraction, aes(y = Sample, x = Fraction, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 0.90, 1),
    labels = c("0%", "25%", "50%", "75%", "90%", "100%"),
    expand = c(0,0)
  ) +
  scale_fill_brewer(palette="Set3")+
  labs(x = "Fraction", y = "Sample", fill = "Type") +
  theme_3()
#  facet_wrap(~ SampleGroup, scales = "free_y")
# theme(
#   axis.text.y = element_blank(),
#   axis.ticks.y = element_blank()
# )
ggsave("./plot/reads_distri.pdf",height = 12,width = 8)
##### 三碱基周期性 #####
df <- fread_c("human_brain_output_20250227/stat/periodicity.txt")
process_1(df) -> df
ggplot(df, aes(x = Periodicity, y = Sample)) +
  geom_col() +
  labs(x = "Periodicity", y = "Sample") +
  theme_3() +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous( expand = c(0,0))+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red")+
  geom_vline(xintercept = 0.65, linetype = "dashed", color = "red")
#  facet_wrap(~ SampleGroup, scales = "free_y")
ggsave("./plot/periodicity.pdf",height = 12,width = 8)

##### Reads不同长度的比例 #####
fread_c("./human_brain_output_20250227/stat/reads_length_distribution.txt") -> df
colnames(df) <- c("Sample","Length","Count")
process_1(df) -> df
filter(df,Length %in% seq(26,34)) -> df
df_norm <- df %>%
  group_by(Sample) %>%
  mutate(proportion = Count / sum(Count)) %>%
  ungroup()
df_matrix <- df_norm %>%
  select(Sample, Length, proportion) %>%
  pivot_wider(names_from = Length, values_from = proportion, values_fill = 0)
rownames_df <- df_matrix$Sample
df_matrix <- df_matrix[,-1]
rownames(df_matrix) <- rownames_df
heat_data <- as.matrix(df_matrix)
pheatmap(
  heat_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "deepskyblue4"))(100),
  main = "Read Length Distribution (Normalized)",
  fontsize_row = 10,
  fontsize_col = 10
) -> p
ggsave(p,filename = "./plot/reads_l_distri.pdf",
       width = 5,height = 10)

heat_data[!grepl("SRR15906",rownames(heat_data)),] -> heat_data_1
pheatmap(
  heat_data_1,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "deepskyblue4"))(100),
  main = "Read Length Distribution (Normalized)",
  fontsize_row = 10,
  fontsize_col = 10
) -> p
ggsave(p,filename = "./plot/reads_l_distri_1.pdf",
       width = 5,height = 10)
          