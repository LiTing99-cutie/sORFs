source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
lib_plot()
setwd("/rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search")
res_path <- "/rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search/output/brain_HLA_20250501/split_20_all_files_20250506/"
# 文件名不同，生成Sample的规则不同，需要修改
read_psm <- function(psm_path) {
  # 检查文件是否只有标题行
  if (length(readLines(psm_path)) == 1) {
    return(NULL)  # 返回空数据框以便后续合并
  }else{
    df <- read.table(psm_path, header = TRUE, sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE)
    # 分割路径并提取psm.tsv的上一级目录名作为样本名
    split_result <- str_split(psm_path, "/")[[1]]
    non_empty_parts <- split_result[split_result != ""]
    target_string <- non_empty_parts[length(non_empty_parts) - 1]
    df$Sample <- target_string
    return(df)
  }
}
get_total_psm(res_path) -> total_psm
files <- list.files(res_path, pattern = "^psm.tsv", full.names = TRUE,recursive = T)
files[1]
# 总鉴定到的谱图数目,440614
nrow(total_psm)
# 总鉴定到的肽段数目
distinct(total_psm,Peptide) %>% nrow()
colnames(total_psm)
head(total_psm[,c("Sample")])
# 根据样本名生成HLA类型以及组织类型列
mutate(total_psm,HLA_type=case_when(
  grepl("W6_32",Sample) ~ "HLA-I",
  grepl("Tue39L243",Sample) ~ "HLA-II"
),Tissue=case_when(
  grepl("Brain",Sample) ~ "Brain",
  grepl("Cerebellum",Sample) ~ "Cerebellum"
)) -> total_psm
head(total_psm[,c("HLA_type","Tissue")])
# 分别查看不同类型的谱图数量
table(total_psm$HLA_type,total_psm$Tissue)

# 加起来一共51376条肽段
total_psm %>% distinct(Peptide) -> unique_peptide
nrow(unique_peptide)
unique_peptide %>% head()

# 单独看看大脑中是否和已经报道的一样，是2w左右的肽段
## 35846
total_psm %>% filter(Tissue=="Brain") %>% distinct(Peptide) %>% nrow()

###### 肽段来自于什么分型 ######
peptide_summary <- total_psm %>%
  group_by(Peptide) %>%
  summarise(
    HLA_I_count = sum(HLA_type == "HLA-I"),
    HLA_II_count = sum(HLA_type == "HLA-II"),
    .groups = 'drop'
  )

peptide_summary <- peptide_summary %>%
  mutate(
    Type = case_when(
      HLA_I_count > 0 & HLA_II_count == 0 ~ "HLA-I only",
      HLA_I_count == 0 & HLA_II_count > 0 ~ "HLA-II only",
      HLA_I_count > 0 & HLA_II_count > 0 ~ "Both HLA-I and HLA-II",
      TRUE ~ "Neither HLA-I nor HLA-II"
    )
  )

table(peptide_summary$Type)

###### 蛋白质来自于什么分型 ######
protein_summary <- total_psm %>% filter(Is.Unique %in% c("true","TRUE")) %>% 
  group_by(Protein) %>%
  summarise(
    HLA_I_count = sum(HLA_type == "HLA-I"),
    HLA_II_count = sum(HLA_type == "HLA-II"),
    .groups = 'drop'
  )

protein_summary <- protein_summary %>%
  mutate(
    Type = case_when(
      HLA_I_count > 0 & HLA_II_count == 0 ~ "HLA-I only",
      HLA_I_count == 0 & HLA_II_count > 0 ~ "HLA-II only",
      HLA_I_count > 0 & HLA_II_count > 0 ~ "Both HLA-I and HLA-II",
      TRUE ~ "Neither HLA-I nor HLA-II"
    )
  )

table(protein_summary$Type)

# 特异性的比对
total_psm %>% filter(Is.Unique %in% c("true","TRUE")) -> unique_psm
# 1174个谱图
unique_psm %>% filter(grepl("ENST",Protein)) -> unique_uncano_sep_psm
# 500个肽段
nrow(distinct(unique_uncano_sep_psm,Peptide))
## 500个肽段来自于什么分型
distinct(unique_uncano_sep_psm,Peptide) %>% merge(peptide_summary) %>% count(Type)
# 495个uncano_sep
nrow(distinct(unique_uncano_sep_psm,Protein))
fwrite(distinct(unique_uncano_sep_psm,Protein),"./stat_output/HLA_peptide_supported_uncano_sep.txt")

##### 肽段分型来源的差异 #####
obs <- c(`Both HLA-I and HLA-II` = 6255, `HLA-I only` = 21137, `HLA-II only` = 23984)
ref <- c(`Both HLA-I and HLA-II` = 17, `HLA-I only` = 292, `HLA-II only` = 191)
# 将参考分布转换为比例
ref_prop <- ref / sum(ref)
# 计算期望频数（基于参考比例）
expected <- ref_prop * sum(obs)
# 执行卡方检验
chisq_test <- chisq.test(x = obs, p = ref_prop)
# 输出结果
print(chisq_test)

library(ggplot2)
# 合并数据
df <- data.frame(
  Type = rep(names(obs), 2),
  Count = c(obs, ref),
  Group = rep(c("Total", "Uncano_SEP"), each = 3)
)
# 绘制堆叠条形图，显示相对比例
ggplot(df, aes(x = Group, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Distribution Comparison", y = "Proportion", x = "Group") +
  scale_y_continuous(labels = scales::percent) +  # 将 y 轴标签转换为百分比格式
  theme_minimal()

###### 从protein.tsv中再查看有多少个非经典蛋白质被鉴定出来 ######
fil_contam <- function(df=df){
  # 读取污染蛋白ID
  fread("/rd1/user/lit/project/sORFs/custom_database/crap.Entry.Name.txt") -> contam_id
  # 将所有 contam_id 合并为一个正则表达式（用"|"分隔）
  combined_pattern <- paste(contam_id$contam_id, collapse = "|")
  # df$`Entry Name` %in% contam_id$contam_id -> idx
  str_detect(string = df$Protein,pattern = combined_pattern) -> idx
  df[!idx,] -> df_fil_contam
  return(df_fil_contam)
}
read_file_fil_contam <- function(path=path){
  # 读取蛋白质结果，并从蛋白质的结果中去掉去掉污染蛋白，并从路径中提取样本名
  df <- fread(path,data.table = F)
  # 分割路径并提取protein.tsv的上一级目录名作为样本名
  split_result <- str_split(path, "/")[[1]]
  non_empty_parts <- split_result[split_result != ""]
  target_string <- non_empty_parts[length(non_empty_parts) - 1]
  if(nrow(df)>1){
    df <- fil_contam(df)
    df$Sample <- target_string
  }
  return(df)
}
get_combined_protein_df_f_path <- function(path){
  protein_files <- list.files(path, pattern = "^protein.tsv", full.names = TRUE,recursive = T)
  lapply(protein_files,read_file_fil_contam) -> l
  merged_protein <- do.call(rbind, l) 
  return(merged_protein)
}
get_combined_protein_df_f_path(res_path) -> total_protein
n_distinct(total_protein$Protein)
total_protein[,"Sample"] %>% head()
total_protein[,"Protein"] %>% head()
colnames(total_protein)
total_protein %>% filter(`Indistinguishable Proteins`=="") -> unique_protein
# 2686
unique_protein %>% filter(grepl("^sp",Protein)) %>% .$Protein %>% unique() %>% length()
unique_protein %>% filter(grepl("^ENST",Protein)) -> unique_uncano_sep_protein
# 560
nrow(distinct(unique_uncano_sep_protein,Protein))
# 495
nrow(distinct(unique_uncano_sep_protein %>% filter(`Unique Peptides`>0),Protein))

##### 看有多少经典的小肽被鉴定出来 #####
read.table("stat_output/Cano.sep.id.txt") -> Cano.sep.id
total_psm %>% filter(grepl("^sp",Protein)) %>% 
  filter(!grepl("sp",Mapped.Proteins)) -> anno_psms
# 2476
length(unique(anno_psms$Protein))
anno_psms %>% filter(Protein %in% Cano.sep.id$V1) -> cano_sep_unique_psm
# 22
n_distinct(cano_sep_unique_psm$Protein)

