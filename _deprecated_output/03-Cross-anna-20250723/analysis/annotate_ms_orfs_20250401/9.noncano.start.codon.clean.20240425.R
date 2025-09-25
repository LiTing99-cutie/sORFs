source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401")

##### 起始密码子附近的kozak序列 ##### 
get_window_15_start_codon <- function(df){
  df %>% filter(Strand=="+") ->  all_sep_forward
  data.frame(all_sep_forward$Chr,all_sep_forward$Start-15,all_sep_forward$Start+18,all_sep_forward$ORF_id_trans,".",
             all_sep_forward$Strand) -> window_15_start_codon_f
  df %>% filter(Strand=="-") ->  all_sep_reverse
  data.frame(all_sep_reverse$Chr,all_sep_reverse$End-18,all_sep_reverse$End+15,all_sep_reverse$ORF_id_trans,".",
             all_sep_reverse$Strand) -> window_15_start_codon_r
  colnames(window_15_start_codon_r) <- colnames(window_15_start_codon_f)
  rbind(window_15_start_codon_f,window_15_start_codon_r) -> window_15_start_codon
  return(window_15_start_codon)
}
## 方法一
# 输入fasta所在的路径
seqlogo_cus <- function(fasta_path,output_path,keep_start_codon_or_not){
  library(seqLogo)
  library(Biostrings)
  # 读取 FASTA 文件
  seqs <- readDNAStringSet(fasta_path)
  # 直接计算碱基频率矩阵（PFM）
  pfm <- consensusMatrix(seqs)[1:4, ]  # 只保留A/C/G/T
  rownames(pfm) <- c("A", "C", "G", "T")
  # 检查PFM（确保无NA/Inf值）
  print(pfm[, 1:5])  # 查看前5列
  if(!keep_start_codon_or_not){
    pfm[,c(1:15,19:33)] -> pfm
  }
  # 转换为PWM（需确保所有列和为1）
  pwm <- makePWM((pfm) / colSums(pfm))  # 加伪计数避免除零错误
  # 绘制Sequence Logo
  pdf(output_path, width = 7, height = 4)  # 可自定义尺寸
  seqLogo(pwm)
  dev.off()
}
## 方法二
# 读取 FASTA 文件
ggseqlogo_cus <- function(fasta_path,output_path,keep_start_codon_or_not){
  library(ggseqlogo)
  library(Biostrings)
  seqs <- readDNAStringSet(fasta_path)
  ggseqlogo(as.character(seqs)) -> p
  if(!keep_start_codon_or_not){
    sub <- sapply(as.character(seqs), function(s) paste0(substr(s, 1, 15), substr(s, 19, 33)))
    ggseqlogo(sub) -> p
  }
  ggsave(p,filename = output_path, width = 7, height = 4)
}
inte_func <- function(tab_info_path,output_path,name,only_cano){
  # 参数一：整理好的蛋白质元信息路径，需要包括起始密码子的位置以及起始密码子的类型（如果需要用到only_cano）
  # 参数二：输出路径
  # 参数三：文件名前缀
  # 参数四：是否只统计起始密码子，需要包括这一列（输入的值为0，1，2）
  fread_c(tab_info_path)->sep
  if(only_cano==1){
    get_window_15_start_codon(sep %>% filter(Scodon=="ATG")) -> window_15_start_codon_sep
  }else if(only_cano==0){
    get_window_15_start_codon(sep %>% filter(Scodon!="ATG")) -> window_15_start_codon_sep
  }else{
    get_window_15_start_codon(sep) -> window_15_start_codon_sep
  }
  create_path(output_path)
  output_bed <- paste0(output_path,name,".bed")
  output_fa <- paste0(output_path,name,".fa")
  output_pdf <- paste0(output_path,name,".pdf")
  fwrite(window_15_start_codon_sep,output_bed,col.names = F,sep = '\t')
  command <- paste0("bedtools getfasta -s -fi /home/user/data/lit/database/public/genome/hg38/hg38.fa -bed ",
                    output_bed," > ",output_fa)
  system(command)
  ggseqlogo_cus(output_fa,output_pdf,TRUE)
  ggseqlogo_cus(output_fa,paste0(output_path,name,".nosc.pdf"),FALSE)
  seqlogo_cus(output_fa, paste0(output_path,name,".seqlogo.pdf"),TRUE)
}
inte_func("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/uniprot.human.sep.tab_info.txt",
          "./output/S9/anno_sep/",
          "prefix",
          2
)
inte_func("./output/S7/all_putative_sep_undetected_sampled.tab_info.txt",
          "./output/S9/sampled_puatative_sep/",
          "prefix",
          2
)
for (i in 0:2){
  inte_func("./output/sep_add_basic_ms_ribo_info_group_retained.txt",
            paste0("./output/S9/ms_detected_novel_sep",i,"/"),
            "prefix",
            i
  )
}
