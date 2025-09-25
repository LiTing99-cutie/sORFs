#!/usr/bin/env Rscript

# 加载包
library(ggplot2)
library(readr)
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
args <- commandArgs(T)
# 读取数据
results <- read_csv(args[1])

# 绘制简单的覆盖曲线
p <- ggplot(results, aes(x = n_samples, y = mean_detected)) +
  geom_line(size = 1.2, color = "blue") +
  geom_point(size = 3, color = "blue") +
  geom_ribbon(aes(ymin = mean_detected - std_detected, 
                  ymax = mean_detected + std_detected), 
              alpha = 0.3, fill = "lightblue") +
  labs(title = "Transcript Coverage vs Number of Samples",
       x = "Number of Samples",
       y = "Number of Detected Transcripts") +
  theme_3()
# 保存为PDF
ggsave(args[2], plot = p, width = 10, height = 6) 