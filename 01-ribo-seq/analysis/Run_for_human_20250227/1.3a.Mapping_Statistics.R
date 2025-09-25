source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

path="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/"
output_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/stat"

# 获取所有 fq.stat.txt 文件的路径
files <- list.files(path=path,pattern = "fq.stat.txt", full.names = TRUE, recursive = TRUE)

# 初始化一个空的列表来存储每个文件处理后的 df
dfs <- list()

# 遍历所有文件
for (file in files) {
  stat <- read.table(file, header = TRUE)
  
  # 处理 stat$num_seqs
  stat$num_seqs <- as.numeric(gsub(",", "", stat$num_seqs))
  
  # 计算 last_r_c, each_step_rm_n, each_step_rm_per
  stat$last_r_c <- c(0, stat$num_seqs[1:(nrow(stat)-1)])
  stat$each_step_rm_n <- stat$last_r_c - stat$num_seqs
  stat$each_step_rm_per <- paste0(round(stat$each_step_rm_n / stat$num_seqs[1] * 100, 2), "%")
  
  # 提取 sample 名称
  sample <- gsub(path, "", file)
  sample <- gsub("/output/fq.stat.txt", "", sample)
  # sample <- file
  
  # 创建当前文件的 df
  df <- matrix(c(sample, stat$each_step_rm_per[2:nrow(stat)]), nrow = 1) %>% as.data.frame()
  # colnames(df) <- c("Sample", "QC", "rRNA", "tRNA", "snoRNA")
  colnames(df) <- c("Sample", "rRNA", "tRNA", "snoRNA")
  
  # 将当前 df 添加到列表中
  dfs[[length(dfs) + 1]] <- df
}

# 使用 do.call 和 rbind 将所有的 df 连接起来
df_all <- do.call(rbind, dfs)

# 打印结果
# print(df_all)

# 写出到文件中
write.table(df_all,o("mapping_statistics.txt"),sep = '\t',row.names = F,quote = F)