source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528")
output_path <- "./output/tmp"
fread_c("./output/S3/sep_info_20250605.txt") -> sep_info
filter(sep_info,Unique_psm_n>0) %>% select(ORF_id_trans) %>% fwrite_c(o("sorf_id.20250606.txt"))
