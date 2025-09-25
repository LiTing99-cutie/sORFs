source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
args <- commandArgs(TRUE)
# 指定路径，读入文件
p_site_count_path <- args[1]
rna_seq_counts_path <- args[2]
rna_seq_lib_path <- args[3]
ribo_seq_lib_path <- args[4]
output_path <- args[5]
if(is.na(args[1])){
  p_site_count_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/output/anno_orfs//ribo_counts/p_site.counts.txt"
  rna_seq_counts_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/output/anno_orfs//rna_counts/"
  rna_seq_lib_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/libsize/rna/libsize.txt"
  ribo_seq_lib_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/libsize/ribo/libsize.txt"
  output_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/output/anno_orfs//te"
}
if(!dir.exists(output_path)){
  dir.create(output_path,recursive = T)
}
fread_c(p_site_count_path) %>% rename(ORF_id_trans=Geneid) -> p_site.counts
p_site.counts %>% .[,7:ncol(p_site.counts)] -> p_site.counts_1
## 读取特定文件夹下以.txt结尾的featureCounts结果文件，并且合并
file_lst <- list.files(rna_seq_counts_path,pattern = "*txt$",full.names = T)
rna_seq_counts_list <- lapply(file_lst, function(path) {
  fread_c(path) %>% rename(ORF_id_trans = Geneid)
})

# 修改列名为样本名
basename(colnames(p_site.counts_1)) -> tmp_str
sapply(strsplit(tmp_str, "\\."), function(x) x[2]) -> sample_name
colnames(p_site.counts_1) <- sample_name

# 修改列名为样本名
lapply(rna_seq_counts_list, function(x)x[,7:ncol(x)]) %>% do.call(cbind,.) -> rna_seq_counts
basename(colnames(rna_seq_counts)) %>% sub("^(.*)_.*$", "\\1", .) -> sample_name
colnames(rna_seq_counts) <- sample_name

# 1、首先计算特定样本的翻译效率（对长度的归一化是一样的）
fread_c(ribo_seq_lib_path) -> libsize
get_rpkm <- function(counts,length,libsize){
  # 输入counts：行为基因，列为样本，基因长度，以及libsize
  counts/length*1000 -> counts_norm_by_l
  col_n <-  ncol(counts_norm_by_l)
  row_n <-  nrow(counts_norm_by_l)
  counts_norm_by_l/matrix(rep(libsize$V2, each = row_n), 
                            nrow = row_n, ncol = col_n)*10^6 -> rpkm
  return(rpkm)
}
get_rpkm(p_site.counts_1,length=p_site.counts$Length,libsize=libsize) -> ribo_rpkm

# 2、其次计算ORF的表达量
## 计算TPM
rna_seq_counts/p_site.counts$Length*1000 -> rna_seq_norm_by_length
rna_seq_norm_by_length %>% apply(., 2, function(x) (x / sum(x))*10^6) -> rna_seq_tpm
add_id <- function(df){
  cbind(ORF_id_trans=p_site.counts$ORF_id_trans,data.frame(df)) -> df_1
  return(df_1)
}
add_id(rna_seq_tpm) -> rna_seq_tpm_1
median_tpm <- apply(rna_seq_tpm, 1,median)
add_id(median_tpm)  %>% rename(median_tpm=df) -> median_tpm_1
fwrite_c(rna_seq_tpm_1,o("rna_seq_tpm.txt"))
fwrite_c(median_tpm_1,o("median_tpm.txt"))
## 计算RPKM
fread_c(rna_seq_lib_path) -> libsize
get_rpkm(rna_seq_counts,length=p_site.counts$Length,libsize=libsize) -> rna_seq_rpkm
add_id(rna_seq_rpkm) -> rna_seq_rpkm_1
median_rpkm <- apply(rna_seq_rpkm, 1,median)
add_id(median_rpkm) %>% rename(median_rpkm=df) -> median_rpkm_1
fwrite_c(rna_seq_rpkm_1,o("rna_seq_rpkm.txt"))
fwrite_c(median_rpkm_1,o("median_rpkm.txt"))

#3、计算翻译效率
readxl::read_excel("/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Mouse_Public/Ribo_RNA_corre_20250220.xlsx",
                   col_names = c("Sample_ribo","Sample_rna","Study"),sheet = 1) -> sample_corre
match(colnames(ribo_rpkm),sample_corre$Sample_ribo) -> idx
sample_corre$Sample_rna[idx] -> sample_order_rna
rna_seq_rpkm %>% dplyr::select(sample_order_rna) -> rna_rpkm_same_order
TE <- ribo_rpkm/(rna_rpkm_same_order+0.001)
median_TE <- apply(TE, 1,median)
add_id(TE) -> TE_1
add_id(median_TE) %>% rename(median_TE=df) -> median_TE_1

fwrite_c(TE_1,o("TE.txt"))
fwrite_c(median_TE_1,o("median_TE.txt"))

