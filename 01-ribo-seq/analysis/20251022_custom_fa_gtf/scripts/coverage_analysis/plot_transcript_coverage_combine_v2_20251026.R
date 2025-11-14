#!/usr/bin/env Rscript

# 加载包
library(ggplot2)
library(readr)
library(dplyr)
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()

args <- commandArgs(T)

# 检查参数
if (length(args) < 3) {
  stop("Usage: Rscript script.R <transcript_results.csv> <gene_results.csv> <output.pdf>")
}

# 读取转录本和基因数据
transcript_results <- read_csv(args[1]) %>%
  mutate(type = "Transcript")

gene_results <- read_csv(args[2]) %>%
  mutate(type = "Gene")

# 合并数据
combined_results <- bind_rows(transcript_results, gene_results)

# 绘制覆盖曲线
p <- ggplot(combined_results, aes(x = n_samples, y = mean_detected, color = type, fill = type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = mean_detected - std_detected, 
                  ymax = mean_detected + std_detected), 
              alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Transcript" = "#2E86AB", "Gene" = "#A23B72")) +
  scale_fill_manual(values = c("Transcript" = "#2E86AB", "Gene" = "#A23B72")) +
  labs(title = "Coverage vs Number of Samples",
       x = "Number of Samples",
       y = "Number of Detected Features",
       color = "Feature Type",
       fill = "Feature Type") +
  theme_3() +
  theme(legend.position = "top")

# 保存为PDF
ggsave(args[3], plot = p, width = 10, height = 6)

cat("Plot saved to:", args[3], "\n")