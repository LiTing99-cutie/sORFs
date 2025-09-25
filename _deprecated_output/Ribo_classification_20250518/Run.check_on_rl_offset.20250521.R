source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()

path="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Ribo_classification_20250518/output/bam/bam_3"
output_path <- path

files <- list.files(path=path,pattern = "metaplots_pre_config.txt", full.names = TRUE, recursive = TRUE)

get_rl_offset <- function(file){
  file %>% stringr::str_split("/") %>% unlist() %>% .[14] -> sample_name
  fread(file,skip = 1,data.table = F) -> peri
  peri %>% dplyr::select(1,3) -> rl_offset
  colnames(rl_offset) <- c("Read_length","Offset")
  rl_offset$Read_length <-  sub("# ","",rl_offset$Read_length)
  rl_offset$Sample <- sample_name
  rl_offset
} 
lapply(files, get_rl_offset) %>% do.call(rbind,.) -> rl_offset_all_s_1
rl_offset_all_s_1$Offset <- as.factor(rl_offset_all_s_1$Offset)
rl_offset_all_s_1 %>% count(Read_length,Offset) %>% bar_plot_basic_stack(.,"Read_length","n",fil_col = "Offset")+scale_fill_discrete_divergingx()
table(rl_offset_all_s_1$Offset)
rl_offset_all_s_1 %>% filter(Offset==12) -> rl_offset_all_s_2
# fwrite(rl_offset_all_s,"./human_brain_output_20250227/stat/read_length_offset.txt")

path="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227"
files <- list.files(path=path,pattern = "metaplots_pre_config.txt", full.names = TRUE, recursive = TRUE)
files[!grepl("cutoff",files)] -> files_1
get_rl_offset <- function(file){
  file %>% stringr::str_split("/") %>% unlist() %>% .[13] -> sample_name
  fread(file,skip = 1,data.table = F) -> peri
  peri %>% dplyr::select(1,3) -> rl_offset
  colnames(rl_offset) <- c("Read_length","Offset")
  rl_offset$Read_length <-  sub("# ","",rl_offset$Read_length)
  rl_offset$Sample <- sample_name
  rl_offset
} 
lapply(files_1, get_rl_offset) %>% do.call(rbind,.) -> rl_offset_all_s_1
rl_offset_all_s_1 %>% count(Read_length,Offset) %>% bar_plot_basic_stack(.,"Read_length","n",fil_col = "Offset")+scale_fill_discrete_divergingx()
table(rl_offset_all_s_1$Offset)


path="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Ribo_classification_20250518/output/bam/bam_4"
files <- list.files(path=path,pattern = "metaplots_pre_config.txt", full.names = TRUE, recursive = TRUE)
lapply(files, get_rl_offset) %>% do.call(rbind,.) -> rl_offset_all_s_1
rl_offset_all_s_1$Offset <- as.factor(rl_offset_all_s_1$Offset)
rl_offset_all_s_1 %>% count(Read_length,Offset) %>% bar_plot_basic_stack(.,"Read_length","n",fil_col = "Offset")+scale_fill_discrete_divergingx()
table(rl_offset_all_s_1$Offset)