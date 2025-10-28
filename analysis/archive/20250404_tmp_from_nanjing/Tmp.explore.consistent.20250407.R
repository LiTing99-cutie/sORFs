# 确认基于psm筛选特异性肽段支持的蛋白质与基于protein.tsv筛选的蛋白质结果是否一致
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
lib_plot()
read_psm("/rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/output/db_search_20250326/Trypsin/CAD20241224licq_BSEP_DDA_60min_21pcw_1_PAGE_T_T_Slot1_46_1_4747_d/psm.tsv")-> psm
psm %>% filter(grepl("ENS",Protein)) %>% filter(Is.Unique=="true") -> sep_unique_psm
sep_unique_psm %>% count(Protein) -> sep_unique_psm_count
fread_c("/rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/output/db_search_20250326/Trypsin/CAD20241224licq_BSEP_DDA_60min_21pcw_1_PAGE_T_T_Slot1_46_1_4747_d/protein.tsv")-> proteins
filter(proteins,`Unique Spectral Count`>0) %>% filter(grepl("ENS",Protein)) -> sep_unique_psm_res

sep_unique_psm_count %>% arrange(Protein) -> tmp_1
sep_unique_psm_res %>% arrange(Protein) -> tmp_2
all(tmp_1$n==tmp_2$`Unique Spectral Count`)
