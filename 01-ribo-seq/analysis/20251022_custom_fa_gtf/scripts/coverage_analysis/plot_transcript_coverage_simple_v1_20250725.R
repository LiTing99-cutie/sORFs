#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# 可选：如果你依赖自定义主题/字体
safely_source <- function(p){
  if (file.exists(p)) try(source(p), silent = TRUE)
}
safely_source("/home/user/data2/lit/bin/lit_utils.R")
safely_source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
if (exists("lib_text")) lib_text()
if (exists("lib_plot")) lib_plot()

args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_saturation_tx_gene.R <tx_csv> <gene_csv> <out.pdf>")
}
tx_file   <- args[1]
gene_file <- args[2]
out_pdf   <- args[3]

# 读取两个结果（要求包含列：n_samples, mean_detected, std_detected）
tx   <- read_csv(tx_file, show_col_types = FALSE)  %>% mutate(level = "Transcript")
gene <- read_csv(gene_file, show_col_types = FALSE) %>% mutate(level = "Gene")

# 检查必需列
need_cols <- c("n_samples","mean_detected","std_detected")
for (nm in need_cols) {
  if (!(nm %in% names(tx)))   stop(paste("Missing column in tx:", nm))
  if (!(nm %in% names(gene))) stop(paste("Missing column in gene:", nm))
}

dat <- bind_rows(tx, gene) %>%
  mutate(
    ymin = mean_detected - std_detected,
    ymax = mean_detected + std_detected
  )

# 作图
p <- ggplot(dat, aes(x = n_samples, y = mean_detected, color = level)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = level),
              alpha = 0.25, color = NA) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  labs(
    title = "饱和曲线（转录本 vs 基因） / Saturation Curves",
    x = "样本数 / Number of Samples",
    y = "检测数量（均值±SD） / Detected (mean ± SD)",
    color = "层级 / Level",
    fill  = "层级 / Level"
  ) +
  (if (exists("theme_3")) theme_3() else theme_classic()) +
  theme(legend.position = "right")

ggsave(out_pdf, plot = p, width = 10, height = 6)
