source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
lib_plot()
setwd("/rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data")

fread_c("./public_processed_data/2024_Wacholder_bioarxiv/SupplementalTable1_all_unannotated_peptides.csv") -> pub
fread_c("./output/all_sample_sep_unique_psm.txt") -> all_sample_sep_unique_psm

# I和L等同
all_sample_sep_unique_psm$Peptide %<>% sub("I","L",.) 
n_distinct(pub$unmodified_peptide)
pub$unmodified_peptide %>% sub("I","L",.) -> pub$unmodified_peptide_n
# 选择非HLA的研究
filter(pub,is_HLA!="TRUE") -> pub_non_HLA
n_distinct(pub_non_HLA$unmodified_peptide)
n_distinct(pub_non_HLA$unmodified_peptide_n)
pub_non_HLA %>% distinct(unmodified_peptide_n,.keep_all = T) -> pub_non_HLA_unique
# 每个研究的肽段数量
bar_plot(pub_non_HLA_unique,"study_pmid")+ggtitle("Number of identified noncanonical peptides")

##### venn plot #####
venn_plot_n(
  list("Pub"=pub_non_HLA_unique$unmodified_peptide_n,
       "In_house"=all_sample_sep_unique_psm$Peptide)
)
##### upset plot #####
df <- pub_non_HLA %>% distinct(study_pmid,unmodified_peptide_n)
# 转换为 UpSet 格式
df_wide <- df %>% 
  mutate(value = 1) %>%
  pivot_wider(names_from = study_pmid, values_from = value, values_fill = 0)
# 去掉肽段列，保留 研究(PMID) 作为列
df_upset <- df_wide %>% select(-unmodified_peptide_n) %>% as.data.frame()
# 绘制 UpSet Plot
library(UpSetR)
upset(df_upset,
      nsets = 10,
      nintersects = NA)
rbind(
  data.frame(study_pmid="In_house",
             unmodified_peptide_n=all_sample_sep_unique_psm$Peptide)%>% distinct(study_pmid,unmodified_peptide_n),
  df
)  -> add_in_house_df

cus_upset <- function(df){
  # 转换为 UpSet 格式
  df_wide <- df %>% 
    mutate(value = 1) %>%
    pivot_wider(names_from = study_pmid, values_from = value, values_fill = 0)
  # 去掉肽段列，保留 研究(PMID) 作为列
  df_upset <- df_wide %>% select(-unmodified_peptide_n) %>% as.data.frame()
  # 绘制 UpSet Plot
  library(UpSetR)
  upset(df_upset,
        nsets = 10,
        nintersects = NA) -> p
  return(p)
}
cus_upset(add_in_house_df)

##### 查看允许一个错配的情况下交集的个数 ##### 
library(stringdist)
peptides_1 <- unique(pub_non_HLA_unique$unmodified_peptide_n)
peptides_2 <- unique(all_sample_sep_unique_psm$Peptide)
find_one_mismatch_pairs <- function(peptides_1, peptides_2) {
  # 去重
  peptides_1 <- unique(peptides_1)
  peptides_2 <- unique(peptides_2)
  
  # 仅保留长度匹配的肽段
  common_lengths <- intersect(nchar(peptides_1), nchar(peptides_2))
  peptides_1 <- peptides_1[nchar(peptides_1) %in% common_lengths]
  peptides_2 <- peptides_2[nchar(peptides_2) %in% common_lengths]
  
  # 构建所有可能组合（可优化为分块执行）
  combinations <- expand.grid(pep1 = peptides_1, pep2 = peptides_2, stringsAsFactors = FALSE)
  combinations <- combinations %>% filter(nchar(pep1) == nchar(pep2))
  
  # 计算 Hamming 距离
  combinations$hamming <- stringdist::stringdist(combinations$pep1, combinations$pep2, method = "hamming")
  
  # 返回 Hamming 距离为1的组合
  mismatch1_df <- combinations %>% filter(hamming == 1)
  
  return(mismatch1_df)
}
result <- find_one_mismatch_pairs(peptides_1, peptides_2)
print(result)

