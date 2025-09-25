fread("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/output/Ribo-ORFs/RiboCode/meta_3_pre_config.txt") -> meta
sum(meta$`proportion(per total mapped reads)`)
sum(meta$f0_sum)+sum(meta$f1_sum)+sum(meta$f2_sum)

fread("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/p21_0321.r_len_distri.txt") -> human_pcw_21_r_len
filter(human_pcw_21_r_len,V2 %in% c(29,30,31,32,33,34,35)) %>% .$V3 %>% sum()

fread("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/stat/reads_length_distribution.txt") %>% 
  filter(V1=="human_brain_ribo_1") %>% 
  filter(.,V2 %in% c(25,28,29)) %>% .$V3 %>% sum()

fread("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/stat/reads_length_distribution.txt") %>% 
  filter(V1=="human_brain_ribo_1") %>% 
  filter(.,V2 %in% c(24,25,26,27,28,29,30)) %>% .$V3 %>% sum()

fread("./p21_0321.r_len_distri.txt") -> df
colnames(df) <- c("Sample","Length","Count")
filter(df,Length %in% seq(26,34)) -> df
df_norm <- df %>%
  mutate(proportion = Count / sum(Count))
