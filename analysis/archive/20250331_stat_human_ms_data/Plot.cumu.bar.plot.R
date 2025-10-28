df <- all_sample_sep_unique_psm_1
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
ggsave(p,filename="plot/cumu_pep_sample_name_add_bar.pdf",height = 5,width = 10)

ggplot(cumu_data_m_sample_unique, aes(x = SampleCount)) +
  geom_line(aes(y = CumulativeUniqueProteins, group = 1), size = 1,color=hcl.colors(3,"Dark 3")[1]) +
  # geom_area(aes(y = CumulativeUniqueProteins, group = 1),fill = "darkblue", alpha = 0.5)+
  geom_bar(aes(y = n), stat = "identity", fill = hcl.colors(3,"Dark 3")[3], width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_3(rotate = T) -> p
ggsave(p,filename="plot/cumu_pep_sample_n_add_bar.pdf",height = 5,width = 6)



