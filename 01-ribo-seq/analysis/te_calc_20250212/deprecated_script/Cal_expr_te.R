# setwd("./01-ribo-seq/analysis/te_calc_20250212/")
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
# 指定路径，读入文件
p_site_count_path <- "./output/anno_orfs//ribo_counts/p_site.counts.txt"
rna_seq_counts_path <- "./output/anno_orfs//rna_counts/"
rna_seq_lib_path <- "./libsize/rna/libsize.txt"
output_path <- "./output/anno_orfs//te"
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

# 首先计算特定样本的翻译效率（由于对长度和文库的归一化，也就是分母都是一样的，因此计算翻译效率只需要将counts相除即可）
readxl::read_excel("../../rawdata/Mouse_Public/Ribo_RNA_corre_20250220.xlsx",
                   col_names = c("Sample_ribo","Sample_rna","Study"),sheet = 1) -> sample_corre
match(colnames(p_site.counts_1),sample_corre$Sample_ribo) -> idx
sample_corre$Sample_rna[idx] -> sample_order_rna
rna_seq_counts %>% dplyr::select(sample_order_rna) -> rna_seq_counts_same_order

TE <- p_site.counts_1/(rna_seq_counts_same_order+0.001)
median_TE <- apply(TE, 1,median)
summary(median_TE)

cbind(ORF_id_trans=p_site.counts$ORF_id_trans,TE) -> TE_1
cbind(ORF_id_trans=p_site.counts$ORF_id_trans,data.frame(median_TE=median_TE)) -> median_TE_1

fwrite_c(TE_1,o("TE.txt"))
fwrite_c(median_TE_1,o("median_TE.txt"))

# 其次计算ORF的表达量
## 计算TPM
rna_seq_counts_same_order/p_site.counts$Length -> rna_seq_norm_by_length
rna_seq_norm_by_length %>% apply(., 2, function(x) (x / sum(x))*10^6) -> rna_seq_tpm
cbind(ORF_id_trans=p_site.counts$ORF_id_trans,data.frame(rna_seq_tpm)) -> rna_seq_tpm_1
fwrite_c(rna_seq_tpm_1,o("rna_seq_tpm.txt"))
## 计算RPKM
fread_c(rna_seq_lib_path) -> libsize
col_n <-  ncol(rna_seq_norm_by_length)
row_n <-  nrow(rna_seq_norm_by_length)
rna_seq_norm_by_length/matrix(rep(libsize$V2, each = row_n), 
                              nrow = row_n, ncol = col_n)*10^6 -> rna_seq_rpkm
cbind(ORF_id_trans=p_site.counts$ORF_id_trans,data.frame(rna_seq_rpkm)) -> rna_seq_rpkm_1
fwrite_c(rna_seq_rpkm_1,o("rna_seq_rpkm.txt"))