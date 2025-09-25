#!/usr/bin/env Rscript
# 功能：从 featureCounts 输出计算 RPKM，并比较核内（N）和核外（C）基因表达差异
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/06-RNA-seq")
###### libsize为UNIQUELY MAPPED READS而不是featureCounts assigned reads ######
fc_output_file <- "/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/featureCounts/rna-seq-counts.txt"
fread_c(fc_output_file)  -> rna_counts
libsize <- fread_c("./02-output/expr/sample_name.lib.txt")
colnames(libsize) <- c("Sample","libsize")
libsize

# 重命名列
colnames(rna_counts[,7:ncol(rna_counts)]) %>% sub("/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/mapping/","",.) %>% 
  sub(".R1_Aligned.sortedByCoord.out.bam","",.) -> colnames_organized
colnames_organized
colnames_organized_1 <- c("C_1","C_2","N")

# 计算RPKM
get_rpkm <- function(counts,libsize){
  rkm <- counts[,7:ncol(counts)]/counts$Length*1000
  rpkm <- data.frame(t(t(rkm) / libsize) * 10^6 )
  return(rpkm)
}
fread_c(fc_output_file)  -> rna_counts
# 确保libsize样本的顺序和counts中样本的顺序是一致的【已确认】
all(libsize$Sample==colnames_organized)
get_rpkm(rna_counts,libsize$libsize) -> rpkm
colnames(rpkm) <- colnames_organized_1

# 计算相关性
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
compare(p,list(c("lncRNA","protein_coding"))) -> p_1
create_path("./02-output/expr/plot/")
ggsave(p_1,filename = "./02-output/expr/plot/lnc_pc_ncRatio.pdf",width = 5,height = 5)

# 导出基因的表达信息
rpkm_N_C %>% mutate(A=(N+C)/2) -> rpkm_N_C_A
fwrite_c(rpkm_N_C_A, "./02-output/expr/rpkm_N_C_A.txt")
