source("~/bin/lit_utils.R")
lib_text()
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")

setwd("/rd1/user/lit/project/sORFs/analysis/20250303_check_sublist_spectrum/")

get_total_psm("/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_default_trans_database") -> total_psm

fread_c("ms_plus_sep_expr_0_uORF_default_trans.txt") -> ms_plus_sep_expr_0_uORF_default_trans

filter(total_psm,`Is.Unique`=="true") %>% filter(Protein %in% ms_plus_sep_expr_0_uORF_default_trans$ORF_id_trans) -> ms_plus_sep_expr_0_uORF_default_trans_psm

select(ms_plus_sep_expr_0_uORF_default_trans_psm,Spectrum,Sample) %>% reformat_spectrum() %>% 
  select(Spectrum,Spectrum_1,Sample) -> ms_plus_sep_expr_0_uORF_default_trans_psm_output

fwrite_c(ms_plus_sep_expr_0_uORF_default_trans_psm_output,"./ms_plus_sep_expr_0_uORF_default_trans_psm_output.txt")
