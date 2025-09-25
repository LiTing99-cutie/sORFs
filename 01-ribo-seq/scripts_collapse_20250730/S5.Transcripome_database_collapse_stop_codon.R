source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

custom <- function(group){
  if(any(group$start_codon=="ATG")){
    group %>% filter(start_codon=="ATG") %>% filter(start_pos==min(start_pos))
  }else{
    group %>% filter(start_pos==min(start_pos))
  }
}

readRDS("./annot/RiboORF/mm39/candidateORF.genepred.rds") -> candidateORF.genepred
colnames(candidateORF.genepred) <- c("ID", "chrom", "strand", "start", "end", "cds_start", "cds_end", "exon_count", "exon_starts", "exon_ends")
str_split_fixed(candidateORF.genepred$ID,pattern = ":|\\|",n = 9) -> id_split
candidateORF.genepred$transcript_id <- id_split[,1]
candidateORF.genepred$start_pos <- id_split[,6]
candidateORF.genepred$stop_pos <- id_split[,7]
candidateORF.genepred$start_codon <- id_split[,9]

Sys.time()
library(furrr)
options(future.globals.maxSize = 1e9)  # 调整为合适的大小
options(future.psock.connectTimeout = 300)
# 设置并行线程数
plan(multisession, workers = 10)  # 根据你的 CPU 核心数调整 workers 参数
# 分组处理并行化
candidateORF.genepred.callapse.stop_codon <- candidateORF.genepred %>%
  group_by(transcript_id, stop_pos) %>%
  group_split() %>%  # 分组为列表
  future_map_dfr(~ custom(.x))  # 并行处理每个组并合并为数据框
# 关闭并行计划
plan(sequential)
Sys.time()

fwrite(candidateORF.genepred.callapse.stop_codon,
       file = "./annot/RiboORF/mm39/candidateORF.genepred.callapse.stop_codon.txt",col.names = F,sep = '\t')