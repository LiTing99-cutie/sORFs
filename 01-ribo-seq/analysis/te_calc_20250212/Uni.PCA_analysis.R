source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()
args <- commandArgs(TRUE)
TE_path <- args[1]
output_path <- args[2]
base <- args[3]
if(is.na(args[1])){
  TE_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/output/anno_orfs/te/rna_seq_rpkm.txt"
  output_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/output/anno_orfs/rna_based_pca_corr_rpkm"
  base <- "rna"
}
if(!dir.exists(output_path)){
  dir.create(output_path,recursive = T)
}
# 行为ORF，列为样本
fread_c(TE_path) -> TE
# 首列用作行名
TE %>% column_to_rownames("ORF_id_trans") -> TE
# PCA分析时的样本分组【不同物种时需要替换】
readxl::read_excel("/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Mouse_Public/Ribo_RNA_corre_20250220.xlsx",
                   col_names = c("Sample_ribo","Sample_rna","Study"),sheet = 1) -> sample_corre
if(base=="rna"){
  sample_corre$Study[match(colnames(TE),sample_corre$Sample_rna)] -> group
}else {
  sample_corre$Study[match(colnames(TE),sample_corre$Sample_ribo)] -> group
}
# 加载必要包
library(FactoMineR)
library(factoextra)
library(corrplot)
# 1. 数据准备 ------------------------------------------------
# 假设数据格式：行名为基因ID，列为样本，数值为表达量
raw_data <- TE

# 2. 变异系数筛选 --------------------------------------------
# 定义变异系数计算函数（处理零均值情况）
calc_cv <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  if(mu <= 1e-6) return(NA)  # 过滤极低表达基因
  sd(x, na.rm = TRUE)/mu
}

# 计算所有基因的变异系数
cv_values <- apply(raw_data, 1, calc_cv) 

# 筛选top 10,000高变异基因
selected_genes <- cv_values %>%
  sort(decreasing = TRUE, na.last = TRUE) %>%
  head(10000) %>%
  names()

# 3. 数据预处理 ----------------------------------------------
# 提取目标基因表达矩阵
processed_data <- raw_data[selected_genes, ] %>% 
  as.matrix() %>%
  t()                     # 转置为样本×基因矩阵

# 4. PCA分析 -------------------------------------------------
pca_result <- PCA(processed_data,
                  scale.unit = TRUE,   # 启用自动标准化
                  ncp = 10,           # 保留前10个主成分
                  graph = FALSE)

# 样本分布图（按实验类型着色）
# 假设存在样本分组信息：metadata包含sample_type列
# 请根据实际metadata修改group_colors映射
fviz_pca_ind(pca_result,
             col.ind = group,
             palette = c("#F8766D", "#00BFC4","black"),
             addEllipses = FALSE, 
             ellipse.level = 0.95,
             legend.title = "Experiment Type",
             title = "Sample Distribution in PCA Space")+
  theme_3() -> p1
ggsave(p1,filename = o("pca.pdf"),height = 6,width = 6)
# 5. 检查样本相关性
apply(TE, 1, sum) -> TE_sum_per_orf
TE[TE_sum_per_orf >0,] -> TE_over_0
cor_matrix <- cor(TE, method = "spearman")
pdf(o("corr_1.pdf"),height = 10,width = 10)
pheatmap::pheatmap(cor_matrix,
                   clustering_method = "ward.D2",
                   main = "Sample Correlation Heatmap",
                   display_numbers = T,fontsize = 10, angle_col = 45)
dev.off()
pdf(o("corr_2.pdf"),height = 10,width = 10)
corrplot::corrplot(cor_matrix, type = "upper", tl.col = "black", order = "hclust", tl.srt = 45, addCoef.col = "white")
dev.off()

pdf(o("boxplot.pdf"),height = 10,width = 10)
boxplot(TE_over_0,las=2,cex.axis=0.6,ylim=c(0,5))
dev.off()

