# 设置工作目录
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()

# 手动填写
readxl::read_excel("./output/S2/Public-RNA-seq样本信息-20250530.xlsx") -> sample_info
colnames(sample_info)
fread("/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/2022-NN/phs002489.v1.pht011444.v1.p1.c1.Human_Brain_Development_Subject_Phenotypes.GRU.txt",skip=10,data.table = F) -> tmp_1
fread_c("/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/2022-NN/SraRunTable.txt") -> tmp_2
colnames(tmp_1)
colnames(tmp_2)
merge(tmp_1,tmp_2,by.x="SUBJECT_ID",by.y="biospecimen_repository_sample_id") -> NN_sample_info
NN_sample_info[,c("Run","AGE","GESTATIONAL_AGE","SUBJECT_ID")] -> NN_sample_info_less
NN_sample_info_less %>% mutate(Developmental_stage=case_when(
  grepl("Fetal",SUBJECT_ID)~ paste0(GESTATIONAL_AGE,"pcw"),
  grepl("Adult",SUBJECT_ID)~ gsub(" ","",AGE)
)) -> NN_sample_info_less
is.na(sample_info$Developmental_stage) -> inx
match(sample_info$Sample[inx],NN_sample_info_less$Run) -> inx_1
sample_info$Developmental_stage[inx] <- NN_sample_info_less$Developmental_stage[inx_1]
sample_info %>%  mutate(Ananomical_region=case_when(
  grepl("pcw",Developmental_stage) & grepl("SRR",Sample)~ "Prenatal cortex",
  grepl("years",Developmental_stage) & grepl("SRR",Sample) ~ "Adult dorsolateral prefrontal cortex",
  TRUE ~ Ananomical_region
)) -> sample_info

sample_info %>% mutate(Developmental_stage_1=case_when(
  grepl("pcw",Developmental_stage)~ "Fetal",
  grepl("years",Developmental_stage)~ "Adult",
  TRUE ~ Developmental_stage
)) -> sample_info

colnames(sample_info)

table(sample_info$Developmental_stage)
table(sample_info$Developmental_stage_1)
table(sample_info$Ananomical_region)

# sample_info %>% fwrite_c("./output/S2/sample_info_rna_seq.txt")

#### 合并文库类型，单双端以及文库大小 ####
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/total_rna_bam_lst.txt") -> libtype
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/fc_output_rna_seq/libsize.txt") -> libsize
merge(libtype,libsize,by.x="BamPath",by.y="V1") -> m
colnames(m) <- c("BamPath","LibraryType","SingleEnd","LibSize")
rename_cus <- function(str){
  # 获取最后一个'/'后的内容
  filename <- basename(str)
  
  # 移除指定后缀
  clean_name <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", filename)
  
  return(clean_name)
}
rename_cus(m$BamPath) -> m$Sample
merge(select(m,-BamPath),sample_info) -> sample_info_1
sample_info_1$LibraryType <- factor(sample_info_1$LibraryType, levels = c("0","1","2"))
sample_info_1$SingleEnd <- factor(sample_info_1$SingleEnd, levels = c("0","1"))
sample_info_1 %>% fwrite_c("./output/S2/sample_info_rna_seq.txt")
sample_info_1 %>% saveRDS("./output/S2/sample_info_rna_seq.rds")
