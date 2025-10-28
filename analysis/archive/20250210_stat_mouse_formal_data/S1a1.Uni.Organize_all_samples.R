# 适配了新的数据
args <- commandArgs(TRUE)

source("~/bin/lit_utils.R")
lib_text()
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
res_path <- args[1]
output_path <- args[2]
sample_metadata_path <- args[3]
if(is.na(args[1])){
  res_path <- "/rd1/user/lit/project/sORFs/analysis/20250224_test_second_run/output/search/open/"
  output_path <- "/rd1/user/lit/project/sORFs/analysis/20250224_test_second_run/output/stat/open/"
  sample_metadata_path <- "/rd1/user/lit/project/sORFs/analysis/20250210_stat_mouse_formal_data/output/sample_metadata.rds"
}
if(!dir.exists(output_path)){
  dir.create(output_path,recursive = T)
}
# 定义函数
get_grouped_protein <- function(df){
# 展开Indistinguishable Proteins中的蛋白
# 输入：protein.tsv数据框，其中Indistinguishable Proteins这一列非空
# 输出：展开了indistinguishable Proteins的数据框，其中蛋白质名字被替换，且subset new_sep出来，但是Length长度需要进一步修改
df %>% separate_rows(`Indistinguishable Proteins`, convert = TRUE,sep=", ") -> tmp
tmp %>% mutate(Protein=`Indistinguishable Proteins`,
              `Protein ID`=`Indistinguishable Proteins`,
              `Entry Name`=`Indistinguishable Proteins`,
              `Protein Description`=`Indistinguishable Proteins`) -> tmp_1
return(tmp_1[grepl("^ENS|^STRG",tmp_1$Protein),])
}

get_leveled_new_sep_df <- function(protein_file){
# 对于每个蛋白质结果文件，得到一个new_sep的df
# 输入：protein.tsv
# 输出：样本中新鉴定到的所有小肽，展开了Indistinguishable Proteins，并且Confidence_Level中标注了可信度
read_file(protein_file) -> protein
# 得到新鉴定的小肽
get_new_sep(protein) -> s_pep
# 对新鉴定的小肽进行分组，分为被unique peptide支持，protein group中的leader protein，protein group中只有new_sep，protein group中两者兼有或者只有注释蛋白
mutate(
  s_pep,
  sp_present = grepl("sp", s_pep$`Indistinguishable Proteins`),
  starts_with_ENS = grepl("^ENS|^STRG", s_pep$`Indistinguishable Proteins`),
  contains_ENS = grepl(" ENS| STRG", s_pep$`Indistinguishable Proteins`),
  Group_type = case_when(
    `Unique Peptides`>0 ~ "Unique",
    `Unique Peptides`==0 & `Indistinguishable Proteins`=="" ~ "Group_none",
    !sp_present & (starts_with_ENS | contains_ENS) ~ "Group_only_new_sep",
    sp_present ~ "Group_mixed"
  )
)  %>% 
select(-sp_present, -starts_with_ENS, -contains_ENS) -> s_pep_add_type
# 得到不同证据支持程度的小肽
s_pep_add_type %>% filter(Group_type=="Unique") %>% mutate(Confidence_Level=1) -> Level_1_new_sep
s_pep_add_type %>% filter(Group_type=="Group_none") %>% mutate(Confidence_Level=2) -> Level_2_new_sep
## 对于protein列和Indistinguishable Proteins同时含有新鉴定小肽的情况，需要把原数据和展开后的数据进行合并
s_pep_add_type %>% filter(Group_type=="Group_only_new_sep") -> df 
rbind(df,get_grouped_protein(df)) %>% mutate(Confidence_Level=3) -> Level_3_new_sep
s_pep_add_type %>% filter(Group_type=="Group_mixed") -> df_1
## 对于protein列为注释蛋白，Indistinguishable Proteins中可能含有新鉴定小肽，展开
protein %>% filter(!grepl("^ENS|^STRG", protein$Protein) & `Indistinguishable Proteins`!="") -> df_2
rbind(df_1,get_grouped_protein(df_1),get_grouped_protein(df_2)) %>% 
  mutate(Confidence_Level=4) -> Level_4_new_sep
rbind(Level_1_new_sep,Level_2_new_sep,Level_3_new_sep,Level_4_new_sep) -> All_level_sep  
All_level_sep %>% select(-Group_type) -> All_level_sep
return(All_level_sep)
}

get_sample_ms_stat_from_path <- function(PATH){
  # 输入：结果路径，子路径中包含了protein.tsv、peptide.tsv、ion.tsv以及psm.tsv，可以输入多个路径构成的向量
  # 输出：每个样本中的质谱统计数据：蛋白总数，肽段总数，以及肽段谱图匹配PSM总数
  protein_files <- c()
  peptide_files <- c()
  ion_files <- c()
  psm_files <- c()
  for (path in PATH){
    protein_files <- c(list.files(path, pattern = "^protein.tsv", full.names = TRUE,recursive = T),protein_files)
    peptide_files <- c(list.files(path, pattern = "^peptide.tsv", full.names = TRUE,recursive = T),peptide_files)
    psm_files <- c(list.files(path, pattern = "^psm.tsv", full.names = TRUE,recursive = T),psm_files)
  }
  df <- data.frame(
    Sample=str_extract(protein_files, "(?<=min_)\\w+(?=_Slot)"),
    Protein_Group_N=sapply(protein_files, function(x){read_file_fil_contam(x) %>% nrow()}),
    Peptdide_N=sapply(peptide_files, function(x){read_file(x) %>% nrow()}),
    PSM_N=sapply(psm_files, function(x){read.table(x, header = TRUE, sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE) %>% nrow()})
  )
  # 原始的行名为文件名，去除
  rownames(df) <- NULL
  return(df)
}

get_stat_from_sample <- function(PATH){
# 输入：结果路径，子路径中包含了protein.tsv，可以输入多个路径构成的向量
# 输出：小肽总体统计；分样本小肽以及质谱数量统计；小肽构成的数据框
protein_files <- c()  
 for (path in PATH){
    protein_files <- c(list.files(path, pattern = "^protein.tsv", full.names = TRUE,recursive = T),protein_files)}
lapply(protein_files, get_leveled_new_sep_df) %>% do.call(rbind,.) -> all_sample_all_level_sep
# 统计小肽结果
## 统计总数
# 统计不同置信度的小肽时，如果在某个样本中置信度是1，但是在其他样本中置信度是2，那往高了选取，在任意一个样本中置信度是1，那么就是1
all_sample_all_level_sep %>% group_by(Protein) %>% filter(Confidence_Level==min(Confidence_Level)) %>% 
  distinct(Protein,.keep_all = T) %>% ungroup() %>% count(Confidence_Level) -> tmp
all_sep_sta <- data.frame(Feature=c("New_sep_level_1_N","New_sep_level_2_N",
                              "New_sep_level_3_N","New_sep_level_4_N","New_sep_total_N"),
                      N=c(tmp$n,n_distinct(all_sample_all_level_sep$Protein)))
all_sep_sta
## 分样本进行统计
all_sample_all_level_sep %>% count(Sample,Confidence_Level) %>% dcast(Sample~Confidence_Level) -> sample_sep_sta
sample_sep_sta[is.na(sample_sep_sta)] <- 0
mutate(sample_sep_sta,New_sep_total_N=rowSums(sample_sep_sta[2:5])) -> sample_sep_sta
colnames(sample_sep_sta) <- c("Sample","New_sep_level_1_N","New_sep_level_2_N",
                              "New_sep_level_3_N","New_sep_level_4_N","New_sep_total_N")
sample_sep_sta
# 统计质谱搜库结果
## 分样本进行统计
get_sample_ms_stat_from_path(PATH) -> sample_ms_sta
merge(sample_sep_sta,sample_ms_sta,by="Sample") -> sample_sep_ms_sta  

return(list(sample_sep_ms_sta=sample_sep_ms_sta,all_sep_sta=all_sep_sta,all_sample_all_level_sep=all_sample_all_level_sep))
}

get_combined_protein_df_f_path <- function(path){
protein_files <- list.files(path, pattern = "^protein.tsv", full.names = TRUE,recursive = T)
lapply(protein_files,read_file_fil_contam,species="mouse") -> l
names(l) <- str_extract(protein_files, "(?<=min_)\\w+(?=_Slot)")
merged_protein <- do.call(rbind, l) 
return(merged_protein)
}

# 在新样本上进行统计，使用函数
get_stat_from_sample(res_path) -> l
# 保存到文件
saveRDS(l,o("stat.rds"))

## 导出
readRDS(sample_metadata_path) -> sample_metadata
merge(l$sample_sep_ms_sta,sample_metadata,by.x="Sample",by.y = "Experiment") %>% mutate(ID_rate=PSM_N/ms2_count) -> sample_sep_ms_sta
fwrite_c(sample_sep_ms_sta,path = o("sample_sep_ms_sta.txt"))
fwrite_c(l$all_sep_sta,path = o("all_sep_sta.txt"))
fwrite_c(l$all_sample_all_level_sep,path = o("all_sample_all_level_sep.txt"))

# 得到所有的protein
get_combined_protein_df_f_path(res_path) -> all_sample_all_protein
fwrite_c(all_sample_all_protein,path = o("all_sample_all_protein.txt"))

# 得到level_1的小肽以及被支持的样本个数
l$all_sample_all_level_sep %>% filter(Confidence_Level==1) -> df_level_1
df_level_1 %>% count(Protein) -> Protein_MS_sample_n
fwrite_c(Protein_MS_sample_n,o("level_1_protein_ms_sample_n.txt"))

