source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/total_ribo_bam_lst.txt") -> libtype
rename_cus <- function(str){
  # 获取最后一个'/'后的内容
  filename <- basename(str)
  
  # 移除指定后缀
  clean_name <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", filename)
  
  return(clean_name)
}
sample_info_ribo <- data.frame(Sample=rename_cus(libtype$BamPath))
# fwrite_c(sample_info_ribo,"./output/S2/sample_info_ribo_seq.txt")
head(sample_info_ribo)
fread_c("/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/2022-NN/SraRunTable.txt") -> tmp_2
filter(tmp_2,`Assay Type`=="RNA-Seq") -> nn_rna
filter(tmp_2,`Assay Type`=="OTHER") -> nn_ribo
merge(nn_rna,nn_ribo,by="biospecimen_repository_sample_id") %>% select(Run.x,Run.y) %>% fwrite_c("./output/S2/ribo_rna_corre.nn.20250602.txt")
fread_c("./output/S2/sample_info_rna_seq.txt") -> sample_info_rna
# 手动制作样本的对应信息
readxl::read_excel("./output/S2/ribo_rna_corre.20250602.xlsx") -> ribo_rna_corre
merge(sample_info_rna,ribo_rna_corre,by.x="Sample",by.y="RNA_for_comparison") -> ribo_rna_sample_info
ribo_rna_sample_info %>% select(-Sample,-RNA,-LibraryType,-SingleEnd,-LibSize) %>% rename(Sample=Ribo) -> sample_info_ribo
sample_info_ribo %>% fwrite_c("./output/S2/sample_info_ribo_seq.txt")

ribo_rna_sample_info %>% select(-LibraryType,-SingleEnd,-LibSize) %>% rename(RNA_1 = Sample) -> ribo_rna_sample_info
ribo_rna_sample_info %>% fwrite_c("./output/S2/sample_info_rna_ribo.txt")
