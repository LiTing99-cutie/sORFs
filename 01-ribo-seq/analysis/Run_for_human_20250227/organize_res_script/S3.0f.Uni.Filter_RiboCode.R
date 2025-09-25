#设置参数指定FDR阈值以及是否过滤codon为ATG
#20241218由于RiboCode的结果文件$sample.txt在不指定起始密码子时没有start_codon列，因此去掉基于start_codon列的过滤
#Usage S3.0f.Uni.Filter_RiboCode.R nonCano.sorf.txt 0.05 0

args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

orfs_path <- args[1]
fdr_cutoff <- as.numeric(args[2])
only_canonical_start_codon <- as.numeric(args[3])
only_canonical_type <- as.numeric(args[4])
output_file <- args[5]

fread_c(orfs_path) -> df
df$adjusted_p_value <- p.adjust(df$pval_combined,method = "BH")
df %>% subset(ORF_length<=450) -> df
if(only_canonical_type){
df %>% subset(ORF_type!="annotated") -> df
}
df %>% subset(adjusted_p_value<=fdr_cutoff) -> df
if(only_canonical_start_codon){
  # df %>% subset(start_codon=="ATG") -> df
  df$start_codon = "ATG"
}

fwrite_c(df,output_file)