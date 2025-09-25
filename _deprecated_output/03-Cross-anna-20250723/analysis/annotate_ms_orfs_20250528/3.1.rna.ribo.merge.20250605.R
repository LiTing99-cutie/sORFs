source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528")
fread_c("./output/S1/sep_info_20250603.txt") -> sep_info

# 函数
rename_cus <- function(str){
  # 获取最后一个'/'后的内容
  filename <- basename(str)
  
  # 移除指定后缀
  clean_name <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", filename)
  
  return(clean_name)
}
fc_output_add_gn <- function(df){
  # 需要有一列列名为GeneID
  ts_meta_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/transcript_meta_output/ts.meta.txt"
  fread_c(ts_meta_path) -> ts_meta
  ts_meta %>% select(2,3) %>% dplyr::rename(GeneID=V2,Gene_name=V3) %>% 
    distinct(GeneID,Gene_name) %>% distinct(Gene_name,.keep_all = T)-> ts_meta_1
  merge(df,ts_meta_1) -> df_add_gene_name
}
# 整理in house的基因表达数据
fread_c("/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/expr/rpkm_N_C_A.txt") -> gene_expr
colnames(gene_expr) <- c("N","C","GeneID","A")
fc_output_add_gn(gene_expr) -> gene_expr_add_gene_name
# 合并基因表达数据
merge(sep_info,gene_expr_add_gene_name,by="Gene_name") -> sep_info_1
# 合并Ribo-seq表达数据
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/fc_output_ribo_seq/ribo_seq_combined_RPKM.txt") -> ribo_rpkm
colnames(ribo_rpkm) <- rename_cus(colnames(ribo_rpkm))
ribo_rpkm[,2:5] %>% rowMeans() -> in_house_mean_rpkm
data.frame(GeneID=ribo_rpkm$GeneID,Ribo_rpkm=in_house_mean_rpkm) -> ribo_rpkm_1
fc_output_add_gn(ribo_rpkm_1) -> ribo_rpkm_2

# 再次看下RNA和Ribo的表达量的相关性【修正了之前的一个小错误】
all(gene_expr_add_gene_name$GeneID==ribo_rpkm_2$GeneID)
merge(gene_expr_add_gene_name,ribo_rpkm_2) -> tmp
cor(tmp$A,tmp$Ribo_rpkm,method = "pearson")
cor(tmp$A,tmp$Ribo_rpkm,method = "spearman")
tmp %>% filter(A>0 & Ribo_rpkm>0) -> tmp_1
cor(tmp_1$A,tmp_1$Ribo_rpkm,method = "spearman")

# 合并
merge(sep_info_1,ribo_rpkm_2,by="Gene_name") -> sep_info_2

create_path("./output/S3/")
fwrite_c(sep_info_2,"./output/S3/sep_info_20250605.txt")
