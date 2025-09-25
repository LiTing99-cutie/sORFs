source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()

setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528")
rna_rpkm_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/fc_output_rna_seq/rna_seq_combined_RPKM.txt"
fread_c(rna_rpkm_file) -> rna_rpkm
metadata_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/sample_info_rna_seq.rds"
readRDS(metadata_file) -> metadata

rename_cus <- function(str){
  # 获取最后一个'/'后的内容
  filename <- basename(str)
  
  # 移除指定后缀
  clean_name <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", filename)
  
  return(clean_name)
}
colnames(rna_rpkm) <- rename_cus(colnames(rna_rpkm))
##### 预处理 #####
rna_rpkm %>% tibble::column_to_rownames("GeneID") %>% as.matrix() -> rpkm_matrix

rpkm_matrix -> expr_matrix
# 假设表达矩阵为expr_matrix，行为基因，列为样本
keep_genes <- rowSums(expr_matrix > 0) > 0  # 保留至少在1个样本中表达的基因
filtered_expr <- expr_matrix[keep_genes, ]
log_expr <- log2(filtered_expr + 1)  # 避免log(0)

##### 热图 #####

# 只展示前30个高变基因
high_var_genes <- names(sort(apply(log_expr, 1, var), decreasing = TRUE)[1:30])
pheatmap(log_expr[high_var_genes, ])

# 只展示前2000个高变基因
high_var_genes_2k <- names(sort(apply(log_expr, 1, var), decreasing = TRUE)[1:2000])
log_expr[high_var_genes_2k, ] -> log_expr_top_var_2k

##### 热图1 #####

library(pheatmap)

# 准备元数据注释（从metadata数据框中选择需要展示的列）
metadata %>% tibble::column_to_rownames("Sample") -> metadata_1
annotation_col <- metadata_1[, c("Source","LibraryType","SingleEnd", "Developmental_stage_1", "Ananomical_region","Developmental_stage")]
annotation_col <- annotation_col[colnames(log_expr_top_var_2k), , drop = FALSE]

annotation_colors <- list(
  Source = brewer.pal(3, "Set1")[1:2] %>% setNames(c("Public", "In_house")),
  LibraryType = brewer.pal(3, "Pastel1") %>% setNames(c("0", "1", "2")),
  SingleEnd = c("0" = "#66C2A5", "1" = "#FC8D62"),  # 绿色 vs 橙色
  Developmental_stage_1 = brewer.pal(3, "Set2") %>% setNames(unique(annotation_col$Developmental_stage_1)),
  Ananomical_region = brewer.pal(5, "Set3") %>%
    setNames(unique(annotation_col$Ananomical_region))
)
# 绘制热图
pheatmap(log_expr_top_var_2k,
         scale = "row",  # 按行标准化
         clustering_distance_rows = "pearson",
         clustering_distance_cols = "pearson",
         clustering_method = "average",
         annotation_col = annotation_col,
         show_rownames = FALSE,  # 基因太多时不显示行名
         main = "Normalized Expression Heatmap",
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100)) -> expr_heatmap
pdf(file = "./output/S2/plot/expr_heatmap.pdf",height = 10,width = 20)
expr_heatmap
dev.off()
###### 测试不同的参数组合###### 
# 定义绘图函数
try <- function(para_1, para_2) {
  pheatmap(log_expr_top_var_2k,
           scale = "row",
           clustering_distance_rows = para_1,
           clustering_distance_cols = para_1,
           clustering_method = para_2,
           annotation_col = annotation_col,
           show_rownames = FALSE,
           main = paste("Dist:", para_1, "\nMethod:", para_2),
           annotation_colors = annotation_colors,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           silent = TRUE)  # 静默模式，不打印信息
}

# 参数组合
para_1 <- c("pearson", "spearman", "euclidean")
para_2 <- c("ward.D", "complete", "average")

# 生成所有组合的热图并转换为grob列表
plot_list <- lapply(para_1, function(p1) {
  lapply(para_2, function(p2) {
    p <- try(p1, p2)
    grid::grid.grabExpr(print(p))  # 关键修正：显式print后捕获
  })
}) |> unlist(recursive = FALSE)

# 排列所有热图（3x3网格）
final_plot <- plot_grid(plotlist = plot_list, ncol = 3, labels = "AUTO")
print(final_plot)

pdf(file = "./output/S2/plot/final_plot.pdf",height = 30,width = 60)
final_plot
dev.off()



##### 热图2 #####

# library(ComplexHeatmap)
# library(circlize)
# # 定义颜色映射
# col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# # 创建列注释
# col_ha <- HeatmapAnnotation(
#   df = metadata_1[, c("Source", "Developmental_stage_1")],
#   col = list(
#     Source = c("Public" = "grey", "In_house" = "orange"),
#     Developmental_stage = c("Adult" = "purple", "Fetal" = "green")
#   )
# )
# # 按行Z-score标准化
# expr_zscore <- t(scale(t(log_expr_top_var_2k)))
# # 绘制热图
# Heatmap(expr_zscore, 
#         name = "log2(expr+1)",
#         col = col_fun,
#         top_annotation = col_ha,
#         show_row_names = FALSE,
#         column_names_gp = gpar(fontsize = 8))

##### 相关性 #####
cor_matrix <- cor(log_expr, method = "spearman")  # 使用spearman或pearson

library(ggcorrplot)

ggcorrplot(cor_matrix,
           hc.order = TRUE,
           type = "full",
           lab = T,
           ggtheme = theme_minimal(),
           title = "Sample-sample Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradientn(
    limits = c(0, 1),              # 关键点：限制颜色范围
    colors = c("white", "red"),
    breaks = seq(0, 1, 0.2),      # 图例刻度
    guide = guide_colorbar(        # 调整图例
      title.position = "top",
      barwidth = unit(3, "cm"))
  ) -> corr_plot
ggsave(corr_plot,filename = "output/S2/plot/corr_plot.pdf",height = 30,width = 30)

# library(heatma()ply)
# heatmaply_cor(
#   cor_matrix,
#   k_col = 2,  # 列聚类数
#   k_row = 2,  # 行聚类数
#   margins = c(100, 100, 40, 20),
#   label_names = c("Sample1", "Sample2", "Correlation"),
#   col_side_colors = metadata[, "Developmental_stage"]
# )
##### PCA #####
library(ggplot2)

# 计算PCA
pca_res <- prcomp(t(log_expr), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, metadata, by = "Sample")

# 绘制PCA
ggplot(pca_df, aes(PC1, PC2, color = Developmental_stage_1,shape = Ananomical_region)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  labs(title = "PCA Plot with Developmental Stage") -> pca_plot
ggsave(pca_plot,filename = "output/S2/plot/pca_plot.pdf",height = 8,width = 8)

##### 合并所有的图片 #####
# library(patchwork)
# (pca_plot + corr_plot) / expr_heatmap
