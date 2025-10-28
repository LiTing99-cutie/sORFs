##### 2、整合结果 #####
# 对于每个有unique peptide的小肽，提取出对应的谱图
path="/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_06_default_ribo_database/Trypsin/CAD20241224licq_BSEP_DDA_60min_E16_1_PAGE_T_T_Slot2_49_1_4736_d"
fread_c(paste0(path,"/protein.tsv")) -> protein
fread_c(paste0(path,"/psm.tsv")) -> psm
fread_c(paste0(path,"/peptide.tsv")) -> peptide
psm %>% filter(`Is Unique`=="TRUE") -> unique_psm
filter(protein,`Unique Peptides`>0 & grepl("^ENS|^STRG", protein$Protein)) -> sep_with_unique_pep
sep_with_unique_pep %>% merge(unique_psm,by="Protein") -> m_1
m_1 %>% count(Protein) -> protein_spectra_n
# check是否一致
merge(protein_spectra_n,sep_with_unique_pep,by="Protein") -> tmp
all(tmp$`Unique Spectral Count`==tmp$n)


reformat_spectrum <- function(psm){
  psm$Spectrum %>% stringr::str_split_fixed(.,"\\.",4) %>% .[,2] %>% gsub("^0+", "", .) -> scan_number
  paste0(stringr::str_split_fixed(psm$Spectrum,"\\.",4) %>% .[,1],".",
         scan_number,".",scan_number,".",psm$Charge) -> spectrum
  return(spectrum)
}
m_1$Spectrum_1 <- reformat_spectrum(m_1)
fwrite_c(m_1[,c("Spectrum","Spectrum_1","Protein")],"output/protein_spectrum/E16_1_PAGE_T_T.txt")


