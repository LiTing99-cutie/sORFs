#!/usr/bin/env Rscript
# Ribo-seq数据分析可视化脚本

args <- commandArgs(T)
# ==================== 初始化设置 ====================
# 加载必要的库和函数
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()

# 设置工作目录和输出路径
setwd(args[1])
output_path <- args[2]
create_path(output_path)

# ==================== 数据文件路径 ====================
input_files <- list(
  unique_mapped = "Uniquely_mapped_reads_rate_number.txt",
  tag_count = "r_distri/summary_tag_count.tsv",
  length_dist = "reads_length_distribution.txt",
  periodicity = "ribotish/Frame_distr_file.txt",
  mapping_stats = "mapping_statistics.txt",
  raw_reads = "raw_reads.txt",
    offset="./ribotish/offset.tab.all.txt"
)



# ==================== 自定义函数 ====================

# 保存图形函数
save_plot <- function(plot, filename, width = 8, height = 12) {
  ggsave(
    filename = file.path(output_path, filename),
    plot = plot,
    width = width,
    height = height
  )
}

# ==================== 数据分析与可视化 ====================

# 1. 不同区域的reads比例分析
analyze_region_distribution <- function(reverse_y = TRUE, out = "reads_distri.pdf") {
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)

  # 1) 读入
  df <- fread_c(input_files$tag_count)

  # 2) 仅保留 CDS / UTR
  df_filtered <- df %>%
    filter(Type %in% c("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons"))

  # 3) 生成样本的目标顺序（批次 + 数值尾号：1..10 再 11..18）
  ord_levels <- df_filtered %>%
    distinct(Sample) %>%
    mutate(
      .batch = str_extract(Sample, "^p21_\\d+"),
      .idx   = as.integer(str_extract(Sample, "(?<=_)\\d+$"))
    ) %>%
    arrange(.batch, .idx) %>%
    pull(Sample)

  if (reverse_y) ord_levels <- rev(ord_levels)  # 需要从上到下递增时，设 TRUE

  # 4) 计算分数（对每个 Sample 的 Type 构成做标准化）
  df_fraction <- df_filtered %>%
    mutate(Sample = factor(Sample, levels = ord_levels)) %>%
    group_by(Sample) %>%
    mutate(
      Total    = sum(Tag_count, na.rm = TRUE),
      Fraction = ifelse(Total > 0, Tag_count / Total, 0)
    ) %>%
    ungroup()
  df_fraction$Type = factor(df_fraction$Type,levels=c("5'UTR_Exons","CDS_Exons",  "3'UTR_Exons"))
  # 5) 绘图（按给定顺序，不聚类不重排）
  p <- ggplot(df_fraction, aes(y = Sample, x = Fraction, fill = Type)) +
    geom_col(width = 0.8) +
    scale_y_discrete(limits = ord_levels) +   # ★ 强制使用我们指定的顺序
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1),
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
    breaks  = c("5'UTR_Exons","CDS_Exons",  "3'UTR_Exons"),
    labels  = c( "5' UTR","CDS", "3' UTR"),
    name    = "Region",
    values = c("CDS_Exons" = "#ffffbc",
               "5'UTR_Exons" = "#9dd1c7",
               "3'UTR_Exons" = "#bebad8")
  )+
    labs(x = "Fraction", y = "Sample", fill = "Region type") +
    theme_3()

  save_plot(p, out)
}

# 2. 三碱基周期性分析
analyze_periodicity <- function() {
  df <- fread_c(input_files$periodicity)
  
  # 处理数据
  df <- df[seq(2, nrow(df), 2), ] %>%
    dplyr::rename(Sample = V1) %>%
    mutate(
      Periodicity = V2 / (V2 + V3 + V4),
      CDS_P_site_number = V2 + V3 + V4
    )
  
  # 保存到全局环境供后续使用
  assign("sample_periodicity", df, envir = .GlobalEnv)
  
  # 绘制周期性柱状图
  p <- ggplot(df, aes(x = Periodicity, y = Sample)) +
    geom_col() +
    labs(x = "Periodicity", y = "Sample") +
    theme_3() +
    scale_x_continuous(expand = c(0, 0)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0.65, linetype = "dashed", color = "red")
  
  save_plot(p, "periodicity.pdf")
}

# 3. Reads长度分布分析
analyze_length_distribution <- function() {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(pheatmap)

  df <- fread_c(input_files$length_dist)
  colnames(df) <- c("Sample", "Length", "Count")

  # 目标顺序（批次 + 尾号的数值）
  ord_levels <- df %>%
    distinct(Sample) %>%
    mutate(
      .batch = str_extract(Sample, "^p21_\\d+"),
      .idx   = as.integer(str_extract(Sample, "(?<=_)\\d+$"))
    ) %>%
    arrange(.batch, .idx) %>%
    pull(Sample)

  # 过滤长度 + 归一化，并按 Sample 目标顺序排列
  df_norm <- df %>%
    filter(Length %in% 26:34) %>%
    group_by(Sample) %>%
    mutate(proportion = Count / sum(Count)) %>%
    ungroup() %>%
    mutate(Sample = factor(Sample, levels = ord_levels)) %>%
    arrange(Sample, Length)

  # 宽表 → 矩阵（行顺序即为 ord_levels）
  df_matrix <- df_norm %>%
    select(Sample, Length, proportion) %>%
    pivot_wider(
      names_from  = Length,
      values_from = proportion,
      values_fill = 0
    )

  rownames_df <- as.character(df_matrix$Sample)
  heat_data <- as.matrix(df_matrix[, -1])
  rownames(heat_data) <- rownames_df

  # 画图（不聚类 -> 保持给定顺序）
  p <- pheatmap(
    heat_data,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "deepskyblue4"))(100),
    fontsize_row = 10,
    fontsize_col = 10
  )

  save_plot(p, "reads_l_distri.pdf", width = 5, height = 10)
}

# 4. 合并所有统计数据
merge_all_stats <- function() {
# 读取各数据文件
  df_1 <- fread_c(input_files$unique_mapped)
  df_5 <- fread_c(input_files$mapping_stats)
  df_5$Sample <- word(df_5$Sample, start = -3, sep = "/")
  # 合并唯一比对和周期统计数据
  tmp <- merge(df_1, df_5)
  tmp_1 <- merge(tmp, sample_periodicity[, c(1, 5, 6)])
  
  # 处理原始reads数
  total_reads <- fread_c(input_files$raw_reads) %>%
    select(1, 4) %>%
    setNames(c("file", "Raw_reads_number")) %>%
    mutate(
      Sample = str_extract(file, "(?<=/)[^/]+(?=\\.fq\\.gz$)"),
      Raw_reads_number = as.numeric(gsub(",", "", Raw_reads_number))
    )
  
  # 合并数据
  tmp_2 <- merge(total_reads[, c(2, 3)], tmp_1)
  tmp_2$CDS_P_site_rate <- tmp_2$CDS_P_site_number / tmp_2$Raw_reads_number
  
  # 合并p位点数据
    fread_c(input_files$offset) -> offset
    fread_c(input_files$length_dist) -> length_distri
    merge(offset,length_distri,by=c("V1","V2")) %>% group_by(V1) %>% summarise(P_site_number=sum(V3.y))  %>% 
    rename("Sample"="V1") -> p_sites_number
  tmp_3 <- merge(p_sites_number, tmp_2)
tmp_3
  tmp_3$P_site_rate <- tmp_3$P_site_number / tmp_3$Raw_reads_number
  
  # 重排列顺序
  col_order <- c(
    "Sample", "Raw_reads_number","rRNA","tRNA","snoRNA", "Uniquely_mapped_reads_number",
    "Uniquely_mapped_reads_rate","Periodicity", "CDS_P_site_number", "CDS_P_site_rate",
     "P_site_number", "P_site_rate"
  )
  
  # 确保所有列都存在
  existing_cols <- intersect(col_order, names(tmp_3))
  final_stats <- tmp_3[, existing_cols]
  
  # 保存结果
  fwrite_c(final_stats, o("stat.all.txt"))
}

# ==================== 执行分析 ====================
analyze_region_distribution()
analyze_periodicity()
analyze_length_distribution()
merge_all_stats()

message("分析完成！结果保存在: ", output_path)