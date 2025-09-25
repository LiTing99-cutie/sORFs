source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()

path="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/"
output_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/stat/"

# 获取所有 fq.stat.txt 文件的路径
files <- list.files(path=path,pattern = "metaplots_pre_config.txt", full.names = TRUE, recursive = TRUE)

get_peri <- function(file){
  file %>% stringr::str_split("/") %>% unlist() %>% .[14] -> sample_name
  fread(file,skip = 1) -> peri
  sum(peri$f0_sum)/(sum(peri$f1_sum)+sum(peri$f2_sum)+sum(peri$f0_sum)) -> periodicity
  return(data.frame(Sample=sample_name,Periodicity=periodicity))
}
lapply(files, get_peri) %>% do.call(rbind,.) -> peri_all_s
peri_all_s[is.na(peri_all_s)] <- 0
fwrite(peri_all_s,"./human_brain_output_20250227/stat/periodicity.txt")

get_rl_offset <- function(file){
  file %>% stringr::str_split("/") %>% unlist() %>% .[14] -> sample_name
  fread(file,skip = 1,data.table = F) -> peri
  peri %>% dplyr::select(1,3) -> rl_offset
  colnames(rl_offset) <- c("Read_length","Offset")
  rl_offset$Read_length <-  sub("# ","",rl_offset$Read_length)
  rl_offset$Sample <- sample_name
  rl_offset
} 
lapply(files, get_rl_offset) %>% do.call(rbind,.) -> rl_offset_all_s
rl_offset_all_s %>% filter(Offset==12) -> rl_offset_all_s
fwrite(rl_offset_all_s,"./human_brain_output_20250227/stat/read_length_offset.txt")
