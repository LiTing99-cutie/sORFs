setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/public_psm/")
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()

fread_c("./2022_NN_PXD035950/total_psm.txt") -> total_psm
read.table("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/uniprot.human.sep.id.txt") -> uniprot.human.sep.id
total_psm %>% filter(grepl("^sp",Protein)) %>% 
  filter(!grepl("sp",Mapped.Proteins)) -> anno_psms
anno_psms %>% filter(Protein %in% uniprot.human.sep.id$V1) -> anno_sep_unique_psm
n_distinct(anno_sep_unique_psm$Protein)
unique(anno_sep_unique_psm$Protein) -> custom_db_protein
# 我们只鉴定出89个
get_ms_info(anno_sep_unique_psm) -> anno_sep_ms_info
mean(anno_sep_ms_info$Unique_psm_n)
mean(anno_sep_ms_info$Unique_peptide_n)
uniprot.human.sep.id$V1[uniprot.human.sep.id$V1 %in% total_psm$Protein] -> p_1
uniprot.human.sep.id$V1[uniprot.human.sep.id$V1 %in% total_psm$Mapped.Proteins] -> p_2
# 126
c(p_1,p_2) %>% unique() %>% length()

total_psm %>% filter(Is.Unique=="TRUE") -> unique_psm
unique_psm %>% filter(grepl("ENST",Protein)) -> res
n_distinct(res$Protein)
fread("../human/new_sep_list.txt") -> sep_list
venn_plot_n(list = list(
  in_house_identified_uncano_sep=sep_list$ORF_id_trans,
  public=unique(res$Protein)))
  
# 和我们鉴定出的188个注释的小肽进行比较
read.table("../human/S3/anno_sep.id.txt") -> in_house_identified_sep
venn_plot_n(list = list(
  in_house_identified_sep=in_house_identified_sep$V1,
  public=unique(anno_sep_unique_psm$Protein)
))

# processed data鉴定出359个uniprot蛋白质
readxl::read_xlsx("../../04-Compare/Public_SEP/human/Duffy et al. 2022_Nat Neurosci.xlsx",sheet = 2) -> prenatal_brain_sep
readxl::read_xlsx("../../04-Compare/Public_SEP/human/Duffy et al. 2022_Nat Neurosci.xlsx",sheet = 3) -> adult_brain_sep
readxl::read_xlsx("../../04-Compare/Public_SEP/human/Duffy et al. 2022_Nat Neurosci.xlsx",sheet = 4) -> neuro_brain_sep
intersect(colnames(prenatal_brain_sep),
          colnames(adult_brain_sep)
          ) %>% intersect(colnames(neuro_brain_sep)) -> comm_cols
rbind(prenatal_brain_sep[,comm_cols],adult_brain_sep[,comm_cols],neuro_brain_sep[,comm_cols]) -> comm_df
comm_df %>% filter(type=="A-canonical") %>% filter(Peptide_Length<=150) %>% distinct(Peptide_Sequence,.keep_all = T) -> processed_data

read.table("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/uniprot.human.sep.rmdup.tab") -> uniprot.human.sep.rmdup.tab

venn_plot_n(list = list(
  processed_data=processed_data$Peptide_Sequence,
  uniprot.human.sep.rmdup.tab=uniprot.human.sep.rmdup.tab$V2
))

# 读入db为uniprot的结果【待修改，可能是因为设置了group】
fread("./2022_NN_PXD035950/total_psm_uniprot_db.txt") -> total_psm
read.table("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/uniprot.human.sep.id.txt") -> uniprot.human.sep.id
total_psm %>% filter(Is.Unique=="TRUE")-> unique_psms
anno_psms %>% filter(Protein %in% uniprot.human.sep.id$V1) -> anno_sep_unique_psm
# 89
n_distinct(anno_sep_unique_psm$Protein)
unique(anno_sep_unique_psm$Protein) -> uniprot_db_protein
get_ms_info(anno_sep_unique_psm) -> anno_sep_ms_info
mean(anno_sep_ms_info$Unique_psm_n)
mean(anno_sep_ms_info$Unique_peptide_n)
uniprot.human.sep.id$V1[uniprot.human.sep.id$V1 %in% total_psm$Protein] -> p_1
uniprot.human.sep.id$V1[uniprot.human.sep.id$V1 %in% total_psm$Mapped.Proteins] -> p_2
# 694
c(p_1,p_2) %>% unique() %>% length()