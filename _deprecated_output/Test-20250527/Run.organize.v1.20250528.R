#!/usr/bin/env Rscript
# Ribo-seq数据分析可视化脚本

# ==================== 初始化设置 ====================
# 加载必要的库和函数
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()

# 设置工作目录和输出路径
setwd("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250527")
output_path <- "01-output/plot/"
create_path(output_path)

# ==================== 数据文件路径 ====================
input_files <- list(
  unique_mapped = "01-output/stat/Uniquely_mapped_reads_rate_number.txt",
  tag_count = "01-output/stat/summary_tag_count.tsv",
  length_dist = "01-output/stat/reads_length_distribution.txt",
  periodicity = "01-output/stat/ribotish/Frame_distr_file.txt",
  mapping_stats = "01-output/stat/mapping_statistics.txt",
  raw_reads = "01-output/stat/raw_reads.txt",
  p_sites = "01-output/stat/p_sites_number.txt",
  orf_number = "01-output/stat/orf_number.txt"
)

# ==================== 自定义函数 ====================
# 自定义绘图主题
my_theme <- theme_3()
    

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
analyze_region_distribution <- function() {
  df <- fread_c(input_files$tag_count)
  
  # 过滤只保留CDS和UTR区域
  df_filtered <- df %>%
    filter(Type %in% c("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons"))
  
  # 计算每个样本的reads比例
  df_fraction <- df_filtered %>%
    group_by(Sample) %>%
    mutate(
      Total = sum(Tag_count),
      Fraction = Tag_count / Total
    )
  
  # 绘制堆叠柱状图
  p <- ggplot(df_fraction, aes(y = Sample, x = Fraction, fill = Type)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 0.90, 1),
      labels = c("0%", "25%", "50%", "75%", "90%", "100%"),
      expand = c(0, 0)
    ) +
    scale_fill_brewer(palette = "Set3") +
    labs(x = "Fraction", y = "Sample", fill = "Region Type") +
    my_theme()
  
  save_plot(p, "reads_distri.pdf")
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
    my_theme() +
    scale_x_continuous(expand = c(0, 0)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0.65, linetype = "dashed", color = "red")
  
  save_plot(p, "periodicity.pdf")
}

# 3. Reads长度分布分析
analyze_length_distribution <- function() {
  df <- fread_c(input_files$length_dist)
  colnames(df) <- c("Sample", "Length", "Count")
  
  # 过滤特定长度范围
  df <- filter(df, Length %in% seq(26, 34))
  
  # 计算标准化比例
  df_norm <- df %>%
    group_by(Sample) %>%
    mutate(proportion = Count / sum(Count)) %>%
    ungroup()
  
  # 转换为宽格式矩阵
  df_matrix <- df_norm %>%
    select(Sample, Length, proportion) %>%
    pivot_wider(
      names_from = Length,
      values_from = proportion,
      values_fill = 0
    )
  
  # 设置行名并转换为矩阵
  rownames_df <- df_matrix$Sample
  df_matrix <- df_matrix[, -1]
  rownames(df_matrix) <- rownames_df
  heat_data <- as.matrix(df_matrix)
  
  # 绘制热图
  p <- pheatmap(
    heat_data,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "deepskyblue4"))(100),
    main = "Read Length Distribution (Normalized)",
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
  df_5$Sample <- str_extract(df_5$Sample, "(?<=call-orfs\\/)[^\\/]+(?=\\/output)")
  
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
  p_sites_number <- fread_c(input_files$p_sites) %>%
    setNames(c("file", "P_site_number")) %>%
    mutate(Sample = str_extract(file, "(?<=/)[^/]+(?=_Aligned\\.sortedByCoord\\.out\\.withPeri\\.bam$)"))
  
  tmp_3 <- merge(p_sites_number[, c(2, 3)], tmp_2)
  tmp_3$P_site_rate <- tmp_3$P_site_number / tmp_3$Raw_reads_number
  
  fread_c(input_files$orf_number) -> orf_number
  colnames(orf_number) <- c("Sample","Identified_ORF_number")
  tmp_3 %>% merge(orf_number) -> tmp_4
  tmp_4$`Identified_ORF_number_per_M_raw_reads` <- tmp_4$Identified_ORF_number/tmp_4$Raw_reads_number*10^6
  
  # 重排列顺序
  col_order <- c(
    "Sample", "Raw_reads_number","rRNA","tRNA","snoRNA", "Uniquely_mapped_reads_number",
    "Uniquely_mapped_reads_rate","Periodicity", "CDS_P_site_number", "CDS_P_site_rate",
     "P_site_number", "P_site_rate" ,"Identified_ORF_number","Identified_ORF_number_per_M_raw_reads"
  )
  
  # 确保所有列都存在
  existing_cols <- intersect(col_order, names(tmp_4))
  final_stats <- tmp_4[, existing_cols]
  
  # 保存结果
  fwrite_c(final_stats, "01-output/stat/stat.all.txt")
}

# ==================== 执行分析 ====================
analyze_region_distribution()
analyze_periodicity()
analyze_length_distribution()
merge_all_stats()

message("分析完成！结果保存在: ", output_path)