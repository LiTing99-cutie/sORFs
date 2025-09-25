source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528")

fread_c("./output/S3/sep_info_20250605.txt") -> sep_info
mutate(sep_info,MS_detection=case_when(
  Unique_psm_n!=0  ~ "MS_unique_peptide",
  Unique_psm_n==0 & All_psm_n!=0 ~ "MS_nonUnique_peptide",
  TRUE ~ "MS_no_detection"
)) -> sep_info

mutate(sep_info,Cano_type=case_when(
  ORF_type_1=="Canonical"  ~ "Cano",
  TRUE ~ "Uncano"
)) -> sep_info

table(sep_info$MS_detection,sep_info$Cano_type)

filter(sep_info,MS_detection=="MS_unique_peptide") -> ms_detection_sep

# MS检测到的SEP是否表达都大于0
## 324
filter(ms_detection_sep,A==0) %>% nrow()
## 293
filter(ms_detection_sep,Ribo_rpkm==0) %>% nrow()

# 1981
filter(ms_detection_sep,A>=1 & Ribo_rpkm>=1) %>% nrow()
# 1462
filter(ms_detection_sep,A>=5 & Ribo_rpkm>=5) %>% nrow()

ms_detection_sep %>% filter(ORF_type_1!="Canonical") -> ms_detection_uncano_sep
# 1493
filter(ms_detection_uncano_sep,A>=1 & Ribo_rpkm>=1) %>% nrow()
# 1147
filter(ms_detection_uncano_sep,A>=5 & Ribo_rpkm>=5) %>% nrow()
