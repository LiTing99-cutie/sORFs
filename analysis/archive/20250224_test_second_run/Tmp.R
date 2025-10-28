source("~/bin/lit_utils.R")
lib_text()
setwd("/rd1/user/lit/project/sORFs/analysis/20250224_test_second_run")
##### 1、得到ID进行second run #####
file_list <- list.files("/rd1/user/lit/project/sORFs/analysis/20250210_stat_mouse_formal_data/output/stat/",pattern = "all_sample_all_level_sep.txt",
           recursive = T,full.names = T)
names <- c("defaul_ribo","default_trans","open_ribo","open_trans")
lapply(1:4, function(x){
  fread_c(file_list[x]) -> df
  df$Method <- names[x]
  return(df)
  }) %>% do.call(rbind,.) -> all_sample_all_level_all_method_sep
# filter(all_sample_all_level_all_method_sep,Sample=="E16_1_PAGE_T_T") -> all_level_all_method_sep
filter(all_sample_all_level_all_method_sep,grepl("T_T",Sample)) -> all_level_all_method_sep

# /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_ribo_database/Trypsin/CAD20241224licq_BSEP_DDA_60min_E16_1_PAGE_T_T_Slot2_49_1_4736_d
# 首先过滤掉第4级别的
filter(all_level_all_method_sep,Confidence_Level!=4) -> all_level_all_method_sep_fil
fwrite(all_level_all_method_sep_fil,col.names = T,file = "./output/all_level_all_method_sep_fil.txt")
# 如果肽段1同时比对到sep1和sep2，肽段2同时比对到sep1和sep3，此时sep1没有unique peptide
# 是否可以只找到那些比对到小肽上的肽段，进行后续蛋白质的组装
all_level_all_method_sep_fil %>% distinct(Protein) -> all_level_all_method_sep_fil_unique
fwrite(all_level_all_method_sep_fil_unique,col.names = F,file = "./output/trypsin_all_method_sep_unique.txt")



