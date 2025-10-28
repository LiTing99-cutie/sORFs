source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
lib_plot()
output_path <- "../results/S2"
create_path(output_path)
fread_c("../results/S1/sample_metadata_ordered.txt") -> sample_metadata
res_path <- "/rd1/user/lit/project/sORFs/analysis/20250523_human_ms_run/output/"
get_sample_ms_stat_from_path(res_path) -> sample_ms_stat
merge(sample_ms_stat,sample_metadata,by="Sample") %>% mutate(ID_rate=PSM_N/ms2_count) -> sample_ms_stat_add_id_rate
# 平均的谱图鉴定率为33.0%（比上次的20.9%要好）
mean(sample_ms_stat_add_id_rate$ID_rate)
fwrite_c(sample_ms_stat_add_id_rate,o("sample_ms_stat_add_id_rate.txt"))
# 初步查看有多少比对到小肽上的特异性肽段，以及多少被支持的小肽
get_total_psm(res_path) -> total_psm
total_psm %>% filter(grepl("^ENS",Protein)) %>% filter(Is.Unique=="true") -> all_sample_sep_unique_psm_1
# 2141个unique谱图支持的蛋白质
length(unique(all_sample_sep_unique_psm_1$Protein))
# 2570个unique谱图支持的肽段
length(unique(all_sample_sep_unique_psm_1$Peptide))
fwrite_c(total_psm,path = o("total_psm.txt"))
fwrite_c(all_sample_sep_unique_psm_1,path = o("all_sample_sep_unique_psm.txt"))