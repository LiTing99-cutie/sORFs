# 加载必要的库
source("~/bin/lit_utils.R")
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
library(FNN)
library(dplyr)

# 定义核心函数
process_peaks <- function(input_path1, input_path2, output_path) {
  # 读取输入文件
  fread_c(input_path1) -> df1
  fread_c(input_path2) -> df2

  # 构造矩阵
  mat1 <- as.matrix(df1[, c("mz", "rt")])
  mat2 <- as.matrix(df2[, c("mz", "rt")])

  # 查找最近邻
  nn <- get.knnx(data = mat1, query = mat2, k = 1)

  # 拼接结果
  res <- cbind(
    df2,
    index = nn$nn.index,
    df1[nn$nn.index, ],
    dist = nn$nn.dist
  )
  colnames(res) <- c("mz2", "rt2", "index_2", "mz1", "rt1", "dist")
  res %>% mutate(index_1 = row_number()) -> res

  # 写入输出文件
  fwrite_c(res, output_path)
}

# 从终端获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) != 3) {
  stop("Usage: Rscript process_peaks.R <input_path1> <input_path2> <output_path>")
}

# 提取参数
input_path1 <- args[1]
input_path2 <- args[2]
output_path <- args[3]

# 调用函数
process_peaks(input_path1, input_path2, output_path)