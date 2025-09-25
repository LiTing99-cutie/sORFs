source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
# setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401")
args <- commandArgs(TRUE)
fc_output_path <- args[1]
if(is.na(args[1])){
  fc_output_file <- "/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/featureCounts/rna-seq-counts.txt"
}
# 1. 计算表达量
add_rpkm <- function(counts,length){
  RPM <- (counts/sum(counts))*10^6
  RPKM <- (RPM/length)*10^3
  return(RPKM)
}
add_tpm <- function(counts, length) {
  # 1. 计算每百万转录本中每千碱基的reads数 (RPK)
  rpk <- counts / (length / 1000)  # 注意长度单位是kb
  # 2. 计算所有基因RPK的总和（按百万归一化）
  per_million_scale_factor <- sum(rpk) / 1e6
  # 3. 计算TPM
  tpm <- rpk / per_million_scale_factor
  return(tpm)
}

fread_c(fc_output_file)  -> rna_counts
colnames(rna_counts)
colSums(rna_counts[,7:ncol(rna_counts)])
apply(rna_counts[,7:ncol(rna_counts)],2,add_rpkm,length=rna_counts$Length) %>% as.data.frame()-> rpkm
apply(rna_counts[,7:ncol(rna_counts)],2,add_tpm,length=rna_counts$Length)%>% as.data.frame() -> tpm

# 2.提取样本名【需要在特定情况下修改】
colnames(rpkm) %>% sub("/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/mapping/","",.) %>% 
  sub(".R1_Aligned.sortedByCoord.out.bam","",.) -> colnames_organized
colnames_organized
colnames_organized <- c("C_1","C_2","N")
colnames(rpkm) <- colnames_organized
colnames(tpm) <- colnames_organized
head(rpkm)
filter(rpkm,C_1+C_2>0) -> rpkm_c_expressed
## 两个重复的相关性特别好，0.999997
cor.test(rpkm_c_expressed$C_1,rpkm_c_expressed$C_2)
cor(rpkm_c_expressed$C_1,rpkm_c_expressed$C_2)
## 核内核外的相关性比较差
cor.test(rpkm_c_expressed$C_1,rpkm_c_expressed$N)
## 合并
mutate(rpkm,C=(C_1+C_2)/2) %>% mutate(C_1=NULL,C_2=NULL) -> rpkm_N_C

# 比较lncRNA和protein-coding基因的NC ratio
rpkm_N_C$Geneid <- rna_counts$Geneid
## 过滤低表达基因
filter(rpkm_N_C,N+C>=0.2) -> rpkm_N_C_expressed
## 导入转录本基因信息
fread("/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt",header = F,data.table = F) -> info
head(info)
distinct(info,V1,V4) %>% merge(rpkm_N_C_expressed,by.x="V1",by.y="Geneid") -> rpkm_N_C_expressed_add_gt
table(rpkm_N_C_expressed_add_gt$V4) %>% sort()
## 比较lncRNA和protein-coding基因
filter(rpkm_N_C_expressed_add_gt,V4 %in% c("lncRNA","protein_coding")) -> rpkm_N_C_expressed_lnc_pc
## 计算log2((N+0.1)/(C+0.1))
rpkm_N_C_expressed_lnc_pc %>% mutate(NC_ratio=log2((N+0.1)/(C+0.1))) -> rpkm_N_C_expressed_lnc_pc
table(rpkm_N_C_expressed_lnc_pc$V4)
box_violin_plot(rpkm_N_C_expressed_lnc_pc,x = "V4",y = "NC_ratio",fill_col ="V4",log10 = F) -> p
p
## 调整violin plot的参数bw
box_violin_plot_v1 <- function(data,x,y,fill_col,log10=T,bw=0.1){
  # log10 是否启用对数坐标轴，注意如果数据中有零值，那么会作图的时候会自动去掉这些值
  ggplot(data,aes_string(x=x,y=y,fill=fill_col))+
    geom_violin(trim = F,bw=bw)+
    geom_boxplot(outlier.shape = NA,fill="white",width=0.2)+
    scale_fill_manual(values = brewer.pal(n=n_distinct(data[,x]),name="Set3"))+
    theme_3() -> p
  if(log10){
    p+scale_y_log10(labels = label_number())->p
  }
  return(p)
}
box_violin_plot_v1(rpkm_N_C_expressed_lnc_pc,x = "V4",y = "NC_ratio",fill_col ="V4",log10 = F,bw=0.5) ->p
p+labs(y="Normalized N/C ratio",fill="Type",x=NULL) -> p
compare(p,list(c("lncRNA","protein_coding")))

###### TPM ######
# 2.提取样本名【需要在特定情况下修改】
filter(tpm,C_1+C_2>0) -> tpm_c_expressed
## 两个重复的相关性特别好，0.999997
cor.test(tpm_c_expressed$C_1,tpm_c_expressed$C_2)
cor(tpm_c_expressed$C_1,tpm_c_expressed$C_2)
## 核内核外的相关性比较差
cor.test(tpm_c_expressed$C_1,tpm_c_expressed$N)
## 合并
mutate(tpm,C=(C_1+C_2)/2) %>% mutate(C_1=NULL,C_2=NULL) -> tpm_N_C

# 比较lncRNA和protein-coding基因的NC ratio
tpm_N_C$Geneid <- rna_counts$Geneid
## 过滤低表达基因
filter(tpm_N_C,N+C>=0.2) -> tpm_N_C_expressed
## 导入转录本基因信息
distinct(info,V1,V4) %>% merge(tpm_N_C_expressed,by.x="V1",by.y="Geneid") -> tpm_N_C_expressed_add_gt
table(tpm_N_C_expressed_add_gt$V4) %>% sort()
## 比较lncRNA和protein-coding基因
filter(tpm_N_C_expressed_add_gt,V4 %in% c("lncRNA","protein_coding")) -> tpm_N_C_expressed_lnc_pc
## 计算log2((N+0.1)/(C+0.1))
tpm_N_C_expressed_lnc_pc %>% mutate(NC_ratio=log2((N+0.1)/(C+0.1))) -> tpm_N_C_expressed_lnc_pc
table(tpm_N_C_expressed_lnc_pc$V4)
box_violin_plot(tpm_N_C_expressed_lnc_pc,x = "V4",y = "NC_ratio",fill_col ="V4",log10 = F) -> p
p


colSums(rna_counts[,7:ncol(rna_counts)])
colSums(rpkm)
rpkm_N_C_expressed_lnc_pc %>% group_by(V4) %>% summarise(median(NC_ratio))
tpm_N_C_expressed_lnc_pc %>% group_by(V4) %>% summarise(median(NC_ratio))
colSums(rpkm_N_C_expressed_lnc_pc[,c(3,4)])
colSums(tpm[,c(1,2,3)])
colSums(tpm_N_C[,c(1,2)])
colSums(rpkm_N_C[,c(1,2)])
