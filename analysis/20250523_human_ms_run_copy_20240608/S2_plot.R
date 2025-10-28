source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
source("~/bin/lit_utils.R")
lib_text()
lib_plot()

setwd("/rd1/user/lit/project/sORFs/analysis/20250523_human_ms_run/")
output_path <- "./stat_output/S2"
create_path(output_path)

fread_c("./stat_output/S1/sample_metadata_ordered.txt") -> sample_metadata
# 去掉PCP第4个重复
sample_metadata %>% filter(Replicate!=4) -> sample_metadata
fread_c(o("all_sample_sep_unique_psm.txt")) -> all_sample_sep_unique_psm
# 去掉PCP第4个重复
all_sample_sep_unique_psm %>% filter(!grepl("4_PCP",all_sample_sep_unique_psm$Sample)) -> all_sample_sep_unique_psm
##### 2. Plot #####
##### 2.0 不同酶切方法以及富集方法鉴定到的小肽数量 #####
library(ggplot2)
all_sample_sep_unique_psm %>% distinct(Sample,Protein,.keep_all = T) -> df_1
df_1 %>% group_by(Sample) %>% summarise(count = n(), .groups = 'drop') %>% merge(sample_metadata,by="Sample") -> df_2
ggplot(df_2, aes(x = Enrichment, y = count, fill = Eyzyme)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.7) +
  stat_summary(
    fun = mean, geom = "point", shape = 18, size = 3, color = "red",
    position = position_dodge(0.8), show.legend = FALSE
  ) +
  labs(x = "Enrichment",
       y = "Identified SEP number",
       fill = "Eyzyme") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") -> p
ggsave(p,filename=o("enrichmen_enzyme_identified_sep_number.pdf"),width = 8,height = 5)
##### 2.1 Saturation Plot #####
set.seed(123)
# sample_order <- sample(sample_metadata$Sample)
sample_order <- sample_metadata %>% arrange(Enrichment,Eyzyme,Replicate) %>% .$Sample
# 1. 蛋白累积分析（带柱状图）
res_protein <- cumu_combo_plot(
  df = all_sample_sep_unique_psm,
  sample_order = sample_order,
  target_col = "Protein",
  show_bar = TRUE
)
ggsave("protein_cumulative_combo.pdf", width = 10, height = 5)
ggsave(res_protein$version_sampleNum,filename=o("cumu_pro_sample_n.pdf"),height = 5,width = 18)
ggsave(res_protein$version_sampleName,filename=o("cumu_pro_sample_name.pdf"),height = 5,width = 18)

# 2. 肽段累积分析（带柱状图）
res_peptide <- cumu_combo_plot(
  df = all_sample_sep_unique_psm,
  sample_order = sample_order,
  target_col = "Peptide",
  show_bar = TRUE
)
ggsave(res_peptide$version_sampleNum,filename=o("cumu_pep_sample_n.pdf"),height = 5,width = 18)
ggsave(res_peptide$version_sampleName,filename=o("cumu_pep_sample_name.pdf"),height = 5,width = 18)

##### 2.2 Overlap between different methods #####
all_sample_sep_unique_psm %>% distinct(Protein,Sample) -> all_sample_sep
all_sample_sep %>% merge(sample_metadata,by = "Sample") -> all_sample_sep_m_meta
colnames(all_sample_sep_m_meta)
##### 2.2.1 Compare enzyme #####
compare_enzyme <- function(method){
  all_sample_sep_m_meta %>% filter(Enrichment==method) %>% filter(Eyzyme!="Null") -> df
  split(df[, "Protein"], df$Eyzyme) -> split_df
  venn_plot_n(split_df) -> p
  path <- paste0("eyzyme_compare_",method,".pdf")
  ggsave(p,filename=o(path),width = 8,height = 5)
  return(p)
}
lapply(c("MWCO","PAGE","C8","PCP"), compare_enzyme) -> l
##### 2.2.2 Compare enrichment method #####
compare_method <- function(enzyme){
  all_sample_sep_m_meta %>% filter(Eyzyme==enzyme) -> df
  split(df[, "Protein"], df$Enrichment) -> split_df
  venn_plot_n(split_df) -> p
  path <- paste0("method_compare_",enzyme,".pdf")
  ggsave(p,filename=o(path),width = 8,height = 5)
  return(p)
}
lapply(c("Trypsin","Trypsin_LysN","Trypsin_Chymotrypsin"), compare_method) -> l

##### 2.2.3 Compare different replicates #####
split(all_sample_sep_m_meta[, "Protein"], all_sample_sep_m_meta$Replicate) -> split_df
venn_plot_n(split_df)->p
ggsave(p,filename=o("replicate_compare.pdf"),width = 8,height = 5)

##### 2.2.4 绘制一个比较大的热图 #####

# Assuming your data is in a dataframe called all_sample_sep_m_meta

# First, create a binary matrix of protein presence (1) vs absence (0)
library(pheatmap)

all_sample_sep_m_meta %>% arrange(Enrichment,Eyzyme,Replicate) -> all_sample_sep_m_meta

# Create a presence/absence matrix
binary_matrix <- all_sample_sep_m_meta %>%
  mutate(present = 1) %>%
  select(Sample, Protein, present) %>%
  pivot_wider(names_from = Sample, values_from = present, values_fill = 0) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

# Create annotation data for columns (samples)
sample_annot <- all_sample_sep_m_meta %>%
  select(Sample,Replicate, Eyzyme,Enrichment) %>%
  distinct() %>%
  mutate(Replicate = as.factor(Replicate)) %>%
  column_to_rownames("Sample")

# Create the heatmap
pheatmap::pheatmap(binary_matrix,
         cluster_cols = FALSE,  # Cluster samples
         cluster_rows = TRUE,  # Cluster proteins
         color = c("white", "black"),  # 0 = white, 1 = black
         legend_breaks = c(0, 1),
         legend_labels = c("Absent", "Present"),
         annotation_col = sample_annot,
         show_rownames = FALSE,  # Typically too many proteins to show all names
         show_colnames = FALSE,  # Typically too many proteins to show all names
         main = "Protein Presence/Absence Across Samples") -> p

ggsave(p,filename=o("protein_sample_detection_heatmap.pdf"),width = 8,height = 8)

##### 2.2.4.1 绘制一个比较大的热图，但是使用complexheatmap #####
library(ComplexHeatmap)
library(circlize)

all_sample_sep_m_meta %>% arrange(Enrichment,Eyzyme,Replicate) -> all_sample_sep_m_meta

# Assuming your data is in a dataframe called all_sample_sep_m_meta
# Create a binary detection matrix
binary_matrix <- all_sample_sep_m_meta %>%
  mutate(present = 1) %>%
  select(Sample, Protein, present) %>%
  pivot_wider(names_from = Sample, values_from = present, values_fill = 0) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

# Calculate detection frequency for each protein
detection_freq <- rowSums(binary_matrix)
detection_number <- colSums(binary_matrix)

# Get unique sample annotations
sample_annot <- all_sample_sep_m_meta %>%
  select(Sample,Enrichment,Eyzyme,Replicate) %>%
  distinct() %>%
  mutate(Replicate = as.factor(Replicate)) %>%
  column_to_rownames("Sample")

# Define colors for annotations
enzyme_colors <- setNames(brewer.pal(length(unique(sample_annot$Eyzyme)), "Set1"), 
                          unique(sample_annot$Eyzyme))
enrichment_colors <- setNames(brewer.pal(length(unique(sample_annot$Enrichment)), "Set2"), 
                              unique(sample_annot$Enrichment))
replicate_colors <- setNames(brewer.pal(length(unique(sample_annot$Replicate)), "Set3"), 
                             unique(sample_annot$Replicate))

# Create column annotations
col_ha <- HeatmapAnnotation(
  df = sample_annot,
  col = list(
    Eyzyme = enzyme_colors,
    Enrichment = enrichment_colors,
    Replicate = replicate_colors
  ),
  annotation_name_side = "left"
)

# Create row annotation for detection frequency
row_ha <- rowAnnotation(
  "Detection Frequency" = anno_barplot(detection_freq, 
                                       bar_width = 0.8,
                                       gp = gpar(fill = "steelblue",col=NULL))
)
col_ha_1 <- columnAnnotation(
  "Detection Number" = anno_barplot(detection_number, 
                                       bar_width = 0.8,
                                       gp = gpar(fill = "steelblue",col=NULL))
)

# Create the heatmap
ht <- Heatmap(
  binary_matrix,
  name = "Detection",
  col = c("white", "black"),
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = col_ha,
  right_annotation = row_ha,
  bottom_annotation = col_ha_1,
  column_title = "Protein Detection Across Samples",
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("Not detected", "Detected")
  ),
  use_raster=FALSE
)

ht_1 <- Heatmap(
  binary_matrix,
  name = "Detection",
  col = c("white", "black"),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = col_ha,
  right_annotation = row_ha,
  bottom_annotation = col_ha_1,
  column_title = "Protein Detection Across Samples",
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("Not detected", "Detected")
  ),
  use_raster=FALSE
)

# Draw the heatmap
# 1. 打开 PDF 设备，设置文件名和尺寸（单位：英寸）
pdf(o("protein_detection_heatmap.pdf"), width = 12, height = 8)
# 2. 绘制热图（假设 `ht` 是你的 Heatmap 对象）
draw(ht)
# 3. 关闭设备，保存文件
dev.off()

# Draw the heatmap
# 1. 打开 PDF 设备，设置文件名和尺寸（单位：英寸）
pdf(o("protein_detection_heatmap_no_cluster_cols.pdf"), width = 12, height = 8)
# 2. 绘制热图（假设 `ht` 是你的 Heatmap 对象）
draw(ht_1)
# 3. 关闭设备，保存文件
dev.off()

##### 2.2.5 不同测定鉴定到的次数 #####
all_sample_sep_m_meta %>% count(Protein) %>% count(n) %>% bar_plot_stated("n","nn")
