# 注意，运行这个脚本的时候，如果只需要rerun作图部分，可以直接导入cor_matrix开始

source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
rna_rpkm_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/fc_output_rna_seq/rna_seq_combined_RPKM.txt"
fread_c(rna_rpkm_file) -> rna_rpkm
metadata_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/sample_info_rna_ribo.txt"
fread_c(metadata_file) -> metadata
ribo_rpkm_file <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/fc_output_ribo_seq/ribo_seq_combined_RPKM.txt"
fread_c(ribo_rpkm_file) -> ribo_rpkm

# 修改metadata的内容
metadata[4:7,"RNA_1"] <- "In_house_p21"
metadata[4:6,"RNA"] <- "In_house_p21"

# 修改rna_rpkm的内容
colnames(rna_rpkm) <- rename_cus(colnames(rna_rpkm))
rna_rpkm$In_house_p21 <- ((rna_rpkm$`L1HJD0800366-p21_1_C.R1`+rna_rpkm$`L1HJD0800367-p21_2_C.R1`)/2+rna_rpkm$`L1HJD2200084-P21.R1`)/2
rna_rpkm %>% mutate(`L1HJD0800366-p21_1_C.R1`=NULL,`L1HJD0800367-p21_2_C.R1`=NULL,`L1HJD2200084-P21.R1`=NULL) -> rna_rpkm

paste0(metadata$Ananomical_region,":",metadata$Developmental_stage) %>% table()
metadata$Group_1 <- paste0(metadata$Ananomical_region,":",metadata$Developmental_stage)

rename_cus <- function(str){
  # 获取最后一个'/'后的内容
  filename <- basename(str)
  
  # 移除指定后缀
  clean_name <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", filename)
  
  return(clean_name)
}
org_1 <- function(rpkm,type){
  # 按照元数据分组取RPKM的平均
  if(type=="ribo"){
    # 处理样本名
    colnames(rpkm) <- sub("_GSM.*RNA-Seq", "", colnames(rpkm))
  }
  colnames(rpkm) <- rename_cus(colnames(rpkm))
  rpkm %>%  pivot_longer(
    cols = -GeneID,              # 转换所有非gene列
    names_to = "Sample",       # 样本名列名
    values_to = "RPKM"         # RPKM值列名
  ) -> rpkm_long
  if(type=="ribo"){
    rpkm_long %>% merge(metadata,by.x="Sample",by.y = "Ribo") -> tmp
  }else{
    rpkm_long %>% merge(metadata,by.x="Sample",by.y = "RNA") -> tmp
  }
  tmp %>% group_by(Group_1,GeneID) %>% summarise(RPKM_mean=mean(RPKM)) ->  group_rpkm
  group_rpkm %>% pivot_wider(names_from = "Group_1",values_from = "RPKM_mean") -> group_rpkm_w
  return(group_rpkm_w)
}

org_1(rna_rpkm,"rna") -> group_rpkm_1
org_1(ribo_rpkm,"ribo") -> group_rpkm_2
# 查看RNA-seq与Ribo-seq之间的相关性【不同研究的两个层次的相关性序列】
cor_cus <- function(GROUP,cutoff=0){
  data.frame(group_rpkm_1[,GROUP],group_rpkm_2[,GROUP]) -> tmp
  colnames(tmp) <- c("RNA","Ribo")
  tmp %>% filter(RNA>cutoff & Ribo>cutoff) -> test
  cor(test$RNA,test$Ribo,method = "spearman")
}
# 如果改成gene name为key会去掉一些
cor_cus("WholeBrain:21pcw",0)
colnames(group_rpkm_2) %>% .[-1] %>% sapply(., cor_cus) %>% sort() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column()

# 查看RNA-seq与Ribo-seq之间的相关性【不同研究的两个层次的相关性矩阵】
matrix1 <- group_rpkm_1[,2:39]
matrix2 <- group_rpkm_2[,2:39]
all(colnames(matrix2)==colnames(matrix1))
intersect(which(rowSums(matrix1)>0),which(rowSums(matrix2)>0)) -> kept_row_n
matrix1[kept_row_n,] -> matrix1_1
matrix2[kept_row_n,] -> matrix2_1
# 计算两个矩阵的列间相关性（38x38）
cor_matrix <- cor(matrix1_1, matrix2_1, method = "spearman")  # 也可以选 "spearman" 或 "kendall"
create_path("./output/S2/stat_output/")
fwrite_c(cor_matrix,"./output/S2/stat_output/cor_matrix_sample_tissue_stage.txt")
library(ggcorrplot)
ggcorrplot(cor_matrix,
           # hc.order = TRUE,
           type = "full",
           lab = T,
           ggtheme = theme_minimal(),
           title = "Sample-sample Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradientn(
    limits = c(min(cor_matrix),max(cor_matrix)),              # 关键点：限制颜色范围
    colors = c("white", "red"),
    breaks = seq(min(cor_matrix),max(cor_matrix), 0.01),      # 图例刻度
    guide = guide_colorbar(        # 调整图例
      title.position = "top",
      barwidth = unit(3, "cm"))
  ) -> corr_plot
ggsave(corr_plot,filename = "output/S2/plot/corr_plot_rna_ribo.pdf",height = 30,width = 30)

cor_matrix <- cor(matrix1_1, matrix2_1, method = "pearson")  # 也可以选 "spearman" 或 "kendall"
library(ggcorrplot)
ggcorrplot(cor_matrix,
           # hc.order = TRUE,
           type = "full",
           lab = T,
           ggtheme = theme_minimal(),
           title = "Sample-sample Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradientn(
    limits = c(min(cor_matrix),max(cor_matrix)),              # 关键点：限制颜色范围
    colors = c("white", "red"),
    breaks = seq(min(cor_matrix),max(cor_matrix), 0.01),      # 图例刻度
    guide = guide_colorbar(        # 调整图例
      title.position = "top",
      barwidth = unit(3, "cm"))
  ) -> corr_plot
ggsave(corr_plot,filename = "output/S2/plot/corr_plot_rna_ribo_pearson.pdf",height = 30,width = 30)

##### 改变分组 ######
org_1 <- function(rpkm, type, group_col = "Group_1") {
  # 按照元数据分组取RPKM的平均
  # 参数说明：
  #   rpkm: 表达矩阵（行为基因，列为样本）
  #   type: "rna" 或 "ribo"，指定数据类型
  #   group_col: 元数据中用于分组的列名（默认为"Group_1"）
  
  # 处理样本名
  if (type == "ribo") {
    colnames(rpkm) <- sub("_GSM.*RNA-Seq", "", colnames(rpkm))
  }
  colnames(rpkm) <- rename_cus(colnames(rpkm))
  
  # 转换为长格式
  rpkm_long <- rpkm %>% 
    pivot_longer(
      cols = -GeneID,
      names_to = "Sample",
      values_to = "RPKM"
    )
  
  # 合并元数据
  merge_by <- ifelse(type == "ribo", "Ribo", "RNA")
  tmp <- rpkm_long %>% 
    merge(metadata, by.x = "Sample", by.y = merge_by)
  
  # 按指定分组列计算均值
  group_rpkm <- tmp %>% 
    group_by(across(all_of(c(group_col, "GeneID")))) %>% 
    summarise(RPKM_mean = mean(RPKM), .groups = "drop")
  
  # 转换回宽格式
  group_rpkm_w <- group_rpkm %>% 
    pivot_wider(
      names_from = all_of(group_col),
      values_from = "RPKM_mean"
    )
  
  return(group_rpkm_w)
}
org_1(rna_rpkm,"rna","Developmental_stage_1") -> group_rpkm_1
org_1(ribo_rpkm,"ribo","Developmental_stage_1") -> group_rpkm_2
# 查看RNA-seq与Ribo-seq之间的相关性【不同研究的两个层次的相关性矩阵】
n <- ncol(group_rpkm_1)
matrix1 <- group_rpkm_1[,2:n]
matrix2 <- group_rpkm_2[,2:n]
intersect(which(rowSums(matrix1)>0),which(rowSums(matrix2)>0)) -> kept_row_n
matrix1[kept_row_n,] -> matrix1_1
matrix2[kept_row_n,] -> matrix2_1
# 计算两个矩阵的列间相关性（38x38）
cor_matrix <- cor(matrix1_1, matrix2_1, method = "spearman")  # 也可以选 "spearman" 或 "kendall"
fwrite_c(cor_matrix,"./output/S2/stat_output/cor_matrix_stage.txt")
ggcorrplot(cor_matrix,
           # hc.order = TRUE,
           type = "full",
           lab = T,
           ggtheme = theme_minimal(),
           title = "Sample-sample Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradientn(
    limits = c(min(cor_matrix),max(cor_matrix)),              # 关键点：限制颜色范围
    colors = c("white", "red"),
    breaks = seq(min(cor_matrix),max(cor_matrix), 0.01),      # 图例刻度
    guide = guide_colorbar(        # 调整图例
      title.position = "top",
      barwidth = unit(3, "cm"))
  ) -> corr_plot_1
ggsave(corr_plot_1,filename = "output/S2/plot/corr_plot_rna_ribo_adult_fetal.pdf",height = 8,width = 8)
