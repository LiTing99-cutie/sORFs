source("~/bin/lit_utils.R")
lib_text()
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")

fread_c("/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_ribo_database/Trypsin/CAD20241224licq_BSEP_DDA_60min_E16_1_3_30K_T_T_Slot2_4_1_4673_d/psm.tsv") -> psm
filter(psm,`Is Unique`=="TRUE") %>% filter(grepl("^ENS",.$Protein)) -> sep_unique_psm
sep_unique_psm$ID <- paste0(sep_unique_psm$Protein,":",sep_unique_psm$Peptide)
fwrite_c(sep_unique_psm[,c("ID","Peptide")],"./sep_unique_psm_peptide.txt",col.names = F)

fread_c("./exact_same_pep.txt") -> blastp_res
distinct(blastp_res,V1) %>% nrow()
count(blastp_res,V1) %>% filter(n==1) %>% nrow()
count(blastp_res,V1) %>% filter(n>1) %>% nrow()

fread_c("sep_unique_psm_peptide.redup.txt") -> dedup_sep
distinct(blastp_res,V1) -> tmp_1

dedup_sep[!dedup_sep$V1 %in% tmp_1$V1,]

fread("out.match_n.txt",header=F) -> out.match_n
as.numeric(str_extract(out.match_n$V2, "\\d+")) -> out.match_n$match_n
