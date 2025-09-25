source("~/bin/lit_utils.R")
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
lib_text()

# 1. 计算表达量
add_rpkm <- function(rna_counts,length){
  counts <- rna_counts
  RPM <- (counts/sum(counts))*10^6
  RPKM <- (RPM/length)*10^3
  return(RPKM)
}

fread_c("./rna_counts_1.txt")  -> rna_counts_1
fread_c("./rna_counts_2.txt")  -> rna_counts_2
fread_c("./rna_counts_3.txt")  -> rna_counts_3
fread_c("./rna_counts_4.txt")  -> rna_counts_4

merge(rna_counts_1,rna_counts_2,by="Geneid") %>% merge(rna_counts_3,by="Geneid") %>% 
  merge(rna_counts_4,by="Geneid") -> merged_rna_counts
merged_rna_counts %>% select(Geneid,merged.1.bam,merged.2.bam,merged.3.bam,merged.4.bam) -> merged_rna_counts
merged_rna_counts$Mean_counts <- apply(merged_rna_counts[,2:5],1,mean)
add_rpkm(merged_rna_counts$Mean_counts,rna_counts_1$Length) -> rpkm
merged_rna_counts$Mean_rpkm <- rpkm
fwrite_c(merged_rna_counts,"./output/merged_rna_counts.txt")

# 2.查看不同批次的数据的相关性
merged_rna_counts %>% filter(Mean_counts>0) -> expressed_rna_counts
library(corrplot)
cor_matrix <- cor(expressed_rna_counts[,2:5], method = "pearson")
pdf("output/sample.gene_expr.corr.pdf",height = 5,width = 5)
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)
dev.off()
