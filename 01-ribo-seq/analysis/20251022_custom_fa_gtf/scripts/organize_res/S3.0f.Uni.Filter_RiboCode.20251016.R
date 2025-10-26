# 设置参数指定FDR阈值以及是否过滤codon为ATG
# 20241218 由于RiboCode结果$sample.txt在不指定起始密码子时无start_codon列，因此去掉基于start_codon列的过滤
# Usage:
#   Rscript S3.0f.Uni.Filter_RiboCode.R nonCano.sorf.txt 0.05 0 0 out.tsv 1 6 450
#   (参数：orfs_path fdr_cutoff only_canonical_start_codon only_canonical_type output_file enable_len_filter min_aa max_aa)

args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

orfs_path <- args[1]
fdr_cutoff <- as.numeric(args[2])
only_canonical_start_codon <- as.numeric(args[3])
only_canonical_type <- as.numeric(args[4])
output_file <- args[5]

# 新增：可选长度过滤（按AA输入，比较时换算为nt=AA*3；未启用则不做长度限制）
enable_len_filter <- ifelse(length(args) >= 6, as.integer(args[6]), 0)
min_aa <- ifelse(length(args) >= 7 && nzchar(args[7]) && args[7] != "NA", as.numeric(args[7]), -Inf)
max_aa <- ifelse(length(args) >= 8 && nzchar(args[8]) && args[8] != "NA", as.numeric(args[8]),  Inf)

fread_c(orfs_path) -> df
df$adjusted_p_value <- p.adjust(df$pval_combined, method = "BH")

# 类型过滤
if (only_canonical_type) {
  df <- subset(df, ORF_type != "annotated")
}

# FDR过滤
df <- subset(df, adjusted_p_value <= fdr_cutoff)

# 起始密码子过滤（按你之前的兼容处理）
if (only_canonical_start_codon) {
  # 原逻辑：df %>% subset(start_codon=="ATG")
  df$start_codon <- "ATG"
}

# 长度过滤（启用时才执行；以AA输入，转为nt比较 ORF_length）
if (enable_len_filter == 1) {
  min_nt <- ifelse(is.finite(min_aa), min_aa * 3, -Inf)
  max_nt <- ifelse(is.finite(max_aa), max_aa * 3,  Inf)
  df <- subset(df, ORF_length >= min_nt & ORF_length <= max_nt)
}

fwrite_c(df, output_file)

# 简要日志
message(sprintf(
  "[INFO] n=%d | FDR<=%.3g | enable_len_filter=%d (min_aa=%s, max_aa=%s -> min_nt=%s, max_nt=%s)",
  nrow(df), fdr_cutoff, enable_len_filter,
  ifelse(is.finite(min_aa), as.character(min_aa), "NA"),
  ifelse(is.finite(max_aa), as.character(max_aa), "NA"),
  ifelse(enable_len_filter==1 && is.finite(min_aa), as.character(min_aa*3), "NA"),
  ifelse(enable_len_filter==1 && is.finite(max_aa), as.character(max_aa*3), "NA")
))
