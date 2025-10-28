##### 展示小肽的unique peptides中FDR控制前后来自decoy的肽段的比例 ######
cus_1 <- function(df){
  Decoy_n <- df %>% filter(grepl("^rev",Protein)) %>% nrow()
  Target_n <-  nrow(df) - Decoy_n
  data.frame(Target_n=Target_n,Decoy_n=Decoy_n)
}
cus_1(sorf_unique_psm_add_decoy) -> df_1
cus_1(sorf_unique_psm_group_specific) -> df_2

library(ggplot2)

# 合并两个数据框并添加组别信息
df_1$group <- "Default FDR"
df_2$group <- "Group_specific FDR"

# 合并数据框
data <- rbind(
  df_1,
  df_2
)

# 将数据转化为长格式，适合 ggplot2 作图
data_long <- tidyr::pivot_longer(data, cols = c("Target_n", "Decoy_n"), 
                                 names_to = "Strategy", values_to = "count")

# 绘图
## 调整柱子之间的距离
n=n_distinct(data_long$group)
alternative_x <- c(rep(0.4,n),rep(0.52,n))
ggplot(data_long, aes(x = alternative_x, y = count, fill = Strategy)) +
  geom_bar(stat = "identity", position = position_stack(0.5),width = 0.1) +
  labs(x = "Strategy", y = "Count",fill="") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_3(rotate = T)+
  scale_fill_brewer(palette = "Accent")+
  scale_x_continuous(
    breaks = unique(alternative_x),
    labels = unique(data_long$group)
  ) 
ggsave(filename = "plot/Stra_decoy_target_n.pdf",width = 3,height = 5)


##### 三次数据手动检查的谱图质量 ######
data.frame(
  Default_FDR=c(13,8,19),
  Group_specific_FDR=c(9,5,6),
  Group_specific_FDR_new_db=c(5,10,5)
) -> df
rownames(df) <- c("Rate_3","Rate_2","Rate_1")
# 转换为长格式
df_long <- df %>%
  mutate(Rate = rownames(df)) %>%
  pivot_longer(cols = -Rate, names_to = "Strategy", values_to = "Count")
df_long %<>% group_by(Strategy) %>% mutate(Proportion = Count / sum(Count))
# 绘图
ggplot(df_long, aes(x = Strategy, y = Proportion, fill = Rate)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Strategy", y = "Proportion", fill = "Rate") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_3(rotate = T)+
  scale_fill_brewer(palette = "Set3")
ggsave(filename = "plot/Manual_inspection_rating.pdf",width = 3,height = 6)
