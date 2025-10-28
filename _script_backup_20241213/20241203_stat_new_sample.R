get_sample_stat_from_path <- function(path=path,stat_feature=stat_feature){
  Sample=str_extract(path, "(?<=min_)\\w+(?=_Slot)")
  if(stat_feature=="Protein"){
    read_file_fil_contam(path) -> df
    df %>% get_new_s_pep -> df_new
    stat <- data.frame(Sample,nrow(df),nrow(df_new))
    colnames(stat) <- c("Sample","Protein_N","New_s_pep_N")
  }else if(stat_feature=="Peptide"){
    read_file(path) -> df
    stat <- data.frame(Sample,nrow(df))
    colnames(stat) <- c("Sample","Peptide_N")
  }
  return(stat)
}

search_path <- "/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024_11_14_batch_run"
protein_files_1 <- list.files(search_path, pattern = "^protein.tsv", full.names = TRUE,recursive = T)
peptide_files_1 <- list.files(search_path, pattern = "^peptide.tsv", full.names = TRUE,recursive = T)

search_path <- "/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/C8_acid/"
protein_files_2 <- list.files(search_path, pattern = "^protein.tsv", full.names = TRUE,recursive = T)
peptide_files_2 <- list.files(search_path, pattern = "^peptide.tsv", full.names = TRUE,recursive = T)

protein_files <- c(protein_files_1,protein_files_2)
peptide_files <- c(peptide_files_1,peptide_files_2)

lapply(protein_files, get_sample_stat_from_path,stat_feature="Protein") %>% do.call(rbind,.) -> stat_protein
lapply(peptide_files, get_sample_stat_from_path,stat_feature="Peptide") %>% do.call(rbind,.) -> stat_peptide
merge(stat_protein,stat_peptide) -> stat

get_combined_protein_df_f_path(search_path) -> merged_df
library(ggpubr)
ggplot(filter(merged_df,Sample=="m14_C8"), aes(x=Length)) +
  geom_density(fill="blue", alpha=0.5) +
  xlab("Protein Length") +
  ylab("Density")+theme_pubr()+
  theme(axis.text.x = element_text(angle=35,vjust = 0.5,hjust = 0.5))+
  xlim(0,1000)
ggplot(filter(merged_df,Sample=="m14_PCP"), aes(x=Length)) +
  geom_density(fill="blue", alpha=0.5) +
  xlab("Protein Length") +
  ylab("Density")+theme_pubr()+
  theme(axis.text.x = element_text(angle=35,vjust = 0.5,hjust = 0.5))+
  xlim(0,1000)