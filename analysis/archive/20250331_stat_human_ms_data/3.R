source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
setwd("/rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data")

fread("./public_processed_data/2024_Wacholder_bioarxiv/SupplementalTable1_all_unannotated_peptides.csv") -> pub
n_distinct(pub$unmodified_peptide)

# I和L等同
pub$unmodified_peptide %>% sub("I","L",.) -> pub$unmodified_peptide_normalized
res_path <- "/rd1/user/lit/project/sORFs/analysis/20250326_human_ms_run/output/db_search_20250326"
get_total_psm(res_path) -> total_psm

######
filter(pub,is_HLA!="TRUE") -> pub_non_HLA
venn_plot_n(
  list("Pub"=pub_non_HLA$unmodified_peptide_normalized,
       "In_house"=all_sample_sep_unique_psm_1$Peptide %>% sub("I","L",.))
)
venn_plot_n(
  list("Pub"=pub_non_HLA$unmodified_peptide_normalized,
       "In_house"=total_psm$Peptide%>% sub("I","L",.))
)
rbind(all_sample_sep_unique_psm_1,all_sample_sep_group_psm) -> psm_1
venn_plot_n(
  list("Pub"=pub_non_HLA$unmodified_peptide_normalized,
       "In_house"=psm_1$Peptide %>% sub("I","L",.))
)
######
pub_non_HLA %>% filter(study_pmid!=36171426) -> pub_non_HLA_non_Duffy
fwrite_c(data.frame(unique(pub_non_HLA_non_Duffy$unmodified_peptide_normalized)),
         "./check_overlap_20250401/pub_non_HLA_non_Duffy.txt",col.names = F)
filter(pub,is_HLA!="TRUE") -> pub_non_HLA
# 和in-house unique psm的交集
venn_plot_n(
  list("Pub"=pub_non_HLA_non_Duffy$unmodified_peptide_normalized,
       "In_house"=all_sample_sep_unique_psm_1$Peptide %>% sub("I","L",.))
) -> p
ggsave(p,filename=o("pub_non_HLA_non_Duffy.in_house_unique_peptide.pdf"),width = 8,height = 5)

# 和in-house total psm的交集
venn_plot_n(
  list("Pub"=pub_non_HLA_non_Duffy$unmodified_peptide_normalized,
       "In_house"=total_psm$Peptide %>% sub("I","L",.))
)
# 和in-house unique_group psm的交集
rbind(all_sample_sep_unique_psm_1,all_sample_sep_group_psm) -> psm_1
venn_plot_n(
  list("Pub"=pub_non_HLA_non_Duffy$unmodified_peptide_normalized,
       "In_house"=psm_1$Peptide %>% sub("I","L",.))
)

######
# 31155234 The Translational Landscape of the Human Heart
# 31340039 A hidden human proteome encoded by 'non-coding' genes
# 33510483 Noncanonical open reading frames encode functional proteins essential for cancer cell survival
# 34663921 Unannotated proteins expand the MHC-I-restricted immunopeptidome in cancer 
# 35841888 A high-resolution map of human RNA translation 
# 非HLA非Duffy的研究的交集
output_path <- "plot/"
## 由于最多只能展示7个研究，因此不展示最少的两个研究
pub_non_HLA_non_Duffy %>% filter(study_pmid!=35788065 & study_pmid!=34193551) -> input
split(input$unmodified_peptide_normalized,input$study_pmid) %>% 
  venn_plot_n() -> p
ggsave(p,filename=o("pub_non_HLA_non_Duffy_compare.pdf"),width = 8,height = 5)
# HLA研究的交集 
filter(pub,is_HLA=="TRUE") -> pub_HLA
split(pub_HLA$unmodified_peptide_normalized,pub_HLA$study_pmid) %>% 
  venn_plot_n()

