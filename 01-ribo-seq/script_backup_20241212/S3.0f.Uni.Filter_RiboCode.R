#设置参数指定FDR阈值以及是否过滤codon为ATG
#Usage S3.0f.Uni.Filter_RiboCode.R nonCano.sorf.txt 0.05 0

args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

orfs_path <- args[1]
fdr_cutoff <- as.numeric(args[2])
only_canonical_start_codon <- as.numeric(args[3])
output_file <- args[4]

fread_c(orfs_path) -> df
df %>% subset(ORF_length<=450 & ORF_type!="annotated") -> df
df$adjusted_p_value <- p.adjust(df$pval_combined,method = "BH")
df %>% subset(adjusted_p_value<=fdr_cutoff) -> df
if(only_canonical_start_codon){
  df %>% subset(start_codon<="ATG") -> df
}

fwrite_c(df,args[4])