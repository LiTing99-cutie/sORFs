source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250509")
output_path <- "01-output/plot/"
create_path(output_path)
file_1 <- "01-output/stat/Uniquely_mapped_reads_rate_number.txt"
file_2 <- "01-output/stat/summary_tag_count.tsv"
file_3 <- "./01-output/stat/reads_length_distribution.txt"
file_4 <- "./01-output/stat/ribotish/total.frame_distr_file.txt"
file_5 <- "./01-output/stat/mapping_statistics.txt"
file_6 <- "./01-output/stat/raw_reads.txt"
file_7 <- "./01-output/stat/p_sites_number.txt"
##### 不同区域的reads比例 #####
# 读取数据
fread_c(file_2) -> df
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
ggsave(o("reads_distri.pdf"),height = 12,width = 8)
##### 三碱基周期性 #####
df <- fread_c(file_4)
df %>% .[c(2,4,6,8),] -> df
samples <- c("p21_40_1_0425","p21_40_0422","p21_0321","p21_40_0425")
df$Sample <- samples
df %>% mutate(Periodicity=V2/(V2+V3+V4)) %>% mutate(CDS_P_site_number=V2+V3+V4) -> df
df -> sample_periodicity
ggplot(df, aes(x = Periodicity, y = Sample)) +
  geom_col() +
  labs(x = "Periodicity", y = "Sample") +
  theme_3() +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous( expand = c(0,0))+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red")+
  geom_vline(xintercept = 0.65, linetype = "dashed", color = "red")
#  facet_wrap(~ SampleGroup, scales = "free_y")
ggsave(o("periodicity.pdf"),height = 12,width = 8)

##### Reads不同长度的比例 #####
fread_c(file_3) -> df
colnames(df) <- c("Sample","Length","Count")
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
ggsave(p,filename = o("reads_l_distri.pdf"),
       width = 5,height = 10)
          
##### uniq比对的reads比例和数目 #####
fread_c(file_1) -> df_1
fread_c(file_5) -> df_5
df_1$Sample <- gsub(".raw","",df_1$Sample)
merge(df_1,df_5) -> tmp
tmp %>% merge(sample_periodicity[,c(5,6,7)]) -> tmp_1
# 需要修改输出路径
# fwrite_c(tmp_1,"./01-output/stat/mapping_statistics_uniq_stat.txt")

fread_c(file_6) %>% select(1,4) -> total_reads
colnames(total_reads) <- c("file","Raw_reads_number")
total_reads$Sample <- total_reads$file %>% str_extract(.,"p21\\w+")
merge(total_reads[,c(2,3)],tmp_1) -> tmp_2
tmp_2$Raw_reads_number <- as.numeric(gsub(pattern = ",",replacement = "",tmp_2$Raw_reads_number))
tmp_2$CDS_P_site_rate <- tmp_2$CDS_P_site_number/tmp_2$Raw_reads_number

fread_c(file_7) -> p_sites_number
colnames(p_sites_number) <- c("file","P_site_number")
p_sites_number$Sample <- p_sites_number$file %>% str_extract(.,"p21\\w+")
merge(p_sites_number[,c(2,3)],tmp_2) -> tmp_3
tmp_3$P_site_rate <- tmp_3$P_site_number/tmp_3$Raw_reads_number
fwrite_c(tmp_3,"./01-output/stat/stat.all.20250516.txt")
