source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
sep_ms_detected_sampled_tab_info_path <- "./output/sep_add_basic_ms_ribo_info_group_retained.txt"
fread_c(sep_ms_detected_sampled_tab_info_path) -> sep_ms_detected_sampled_tab_info
sep_ms_detected_sampled_tab_info -> df
df %>% filter(ORF_type_1!="Canonical") -> df
output_path <- "output/S2/plot/"
create_path(output_path)
## 不同的ORF类型
bar_plot_v1(df,col = "ORF_type_1",fil_col = "ORF_type_1",label = T,ylim = 750)+
  scale_fill_manual(values = hcl.colors(7,"Dark 3"))  -> p
ggsave(p,filename = o("ORF_type.pdf"),height = 5,width = 5)

## start codon
pie_plot(df,"Scodon") -> p
ggsave(p,filename = o("scodon_type.pdf"),height = 5,width = 5)

## ms evidence
df %>% mutate(Unique_peptide_n_1=
                case_when(Unique_peptide_n > 2 ~ ">2",
                          TRUE ~ as.character(Unique_peptide_n))) -> df
pie_plot(df,"Unique_peptide_n_1")+
  guides(fill = guide_legend(title = "Unique peptide number")) -> p
ggsave(p,filename = o("unique_peptide_n.pdf"),height = 5,width = 5)

