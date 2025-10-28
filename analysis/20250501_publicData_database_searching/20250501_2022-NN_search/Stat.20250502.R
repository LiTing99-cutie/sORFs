source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
lib_plot()
setwd("/rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search")
res_path <- "/rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search/output/Trypsin/"
get_total_psm(res_path) -> total_psm
create_path("stat_output")
fwrite_c(total_psm,path = "stat_output/total_psm.txt")

res_path <- "/rd1/user/lit/project/sORFs/analysis/20250501_publicData_database_searching/20250501_2022-NN_search/output/uniprot_db/"
get_total_psm(res_path) -> total_psm
create_path("stat_output")
fwrite_c(total_psm,path = "stat_output/total_psm_uniprot_db.txt")