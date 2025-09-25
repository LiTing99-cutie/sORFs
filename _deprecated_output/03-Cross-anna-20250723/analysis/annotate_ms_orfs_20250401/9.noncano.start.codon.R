source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
# 想看下非经典的起始密码子被支持的肽段在SEP上的分布
fread_c("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/psm_retained_group_based_on_transcript.txt") -> psm
fread_c("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/sep_unique_psm.txt") -> psm_1
cols <- c("Protein","Peptide")
rbind(psm[,cols],psm_1[,cols]) %>% distinct(Protein,Peptide) -> pro_pep
n_distinct(pro_pep$Protein)
fread("./output/sep_add_basic_ms_ribo_info_group_retained.txt") -> all_sep
head(all_sep)
pro_pep %>% rename(ORF_id_trans=Protein) -> tmp
update_orf_id_trans(tmp) -> pro_pep_1
merge(pro_pep_1,all_sep,by="ORF_id_trans") -> all_sep_ms_pep

# 示例数据
protein_sequence <- "ABCDCC"
peptide_sequence <- "ABCD"
# 函数：计算每个位置被支持的次数
calculate_coverage_relative_to_N <- function(protein, peptide) {
  protein_len <- nchar(protein)
  peptide_len <- nchar(peptide)
  # 初始化结果向量
  coverage <- rep(0, protein_len)
  # 在蛋白质序列中查找肽段的所有可能匹配
  for (i in 1:(protein_len - peptide_len + 1)) {
    substring <- substr(protein, i, i + peptide_len - 1)
    if (substring == peptide) {
      coverage[i:(i + peptide_len - 1)] <- coverage[i:(i + peptide_len - 1)] + 1
    }
  }
  # 创建输出格式
  result <- c()
  for (pos in 1:protein_len) {
    result <- c(result, pos, coverage[pos])
  }
  matrix(result,nrow = 2) %>% t() %>% as.data.frame() -> df
  colnames(df) <- c("Position","Coverage")
  df$Position <- df$Position-1
  return(df)
}
calculate_coverage_relative_to_C<- function(protein, peptide){
  df <- calculate_coverage_relative_to_N(protein,peptide) -> df
  df$Position <- df$Position - nchar(protein)+1
  return(df)
}
# 使用函数
calculate_coverage_relative_to_C(protein_sequence,peptide_sequence)
calculate_coverage_relative_to_N(protein_sequence,peptide_sequence)

cumu_pos_cov <- function(protein_sequence,peptide_sequence,CorN){
  # 输入数据库框的两列
  if(CorN=="C"){
    purrr::pmap(list(protein_sequence,peptide_sequence),calculate_coverage_relative_to_C) %>% Reduce(rbind,.) -> pos_cov
  }else if(CorN=="N"){
    purrr::pmap(list(protein_sequence,peptide_sequence),calculate_coverage_relative_to_N) %>% Reduce(rbind,.) -> pos_cov
  }
  pos_cov %>% group_by(Position) %>% summarise(Cov_sum=sum(Coverage)) -> pos_cov_sum
  return(pos_cov_sum)
}

# 所有SEP N端开始的位置index以及覆盖度
cumu_pos_cov(all_sep_ms_pep$Seq,all_sep_ms_pep$Peptide,"N") -> pos_cov_sum_N
cumu_pos_cov(all_sep_ms_pep$Seq,all_sep_ms_pep$Peptide,"C") -> pos_cov_sum_C
pos_cov_sum_N %>% head(26) -> df_1
pos_cov_sum_C %>% tail(26) -> df_2
max <- max(df_1$Cov_sum,df_2$Cov_sum)
bar_plot_basic(df_1,"Position","Cov_sum")+ylim(0,1500)+labs(y="Peptide coverage") -> p1
bar_plot_basic(df_2,"Position","Cov_sum")+ylim(0,1500)+labs(y="Peptide coverage") -> p2
wrap_plots(list(p1,p2))

# 非经典起始密码子的SEP N端开始的位置index以及覆盖度
plot_cus <- function(df){
  cumu_pos_cov(df$Seq,df$Peptide,"N") -> pos_cov_sum_N
  cumu_pos_cov(df$Seq,df$Peptide,"C") -> pos_cov_sum_C
  pos_cov_sum_N %>% head(26) -> df_1
  pos_cov_sum_C %>% tail(26) -> df_2
  max <- max(df_1$Cov_sum,df_2$Cov_sum)
  bar_plot_basic(df_1,"Position","Cov_sum")+ylim(0,1500)+labs(y="Peptide coverage") -> p1
  bar_plot_basic(df_2,"Position","Cov_sum")+ylim(0,1500)+labs(y="Peptide coverage") -> p2
  wrap_plots(list(p1,p2)) -> p
  return(p)
}
plot_cus(all_sep_ms_pep %>% filter(Scodon!="ATG")) %>% ggsave(.,filename = "./output/S9/plot/noncano.peptide.cov.pdf",width = 6,height = 4)
plot_cus(all_sep_ms_pep %>% filter(Scodon=="ATG")) %>% ggsave(.,filename = "./output/S9/plot/cano.peptide.cov.pdf",width = 6,height = 4)
plot_cus(all_sep_ms_pep) %>% ggsave(.,filename = "./output/S9/plot/all.peptide.cov.pdf",width = 6,height = 4)

# # 查看有多少肽段从N端的第一个位置开始覆盖（M没有被CLIP）或者从N端的第二个位置开始覆盖（M被CLIP）
# cus_2 <- function(df){
#   purrr::pmap(list(df$Seq,df$Peptide),calculate_coverage_relative_to_N) %>% Reduce(rbind,.) -> pos_cov_1
#   pos_cov_1 %>% group_by(Position) %>% summarise(Cov_sum=sum(Coverage)) -> pos_cov_sum_1
#   
#   purrr::pmap(list(df$Seq,df$Peptide),calculate_coverage_relative_to_N) -> l
#   # 是否肽段从N段的第一个位置开始覆盖
#   cus_1 <- function(result){
#     return(result %>% head(1) %>% .$Coverage ==1)
#   }
#   # 是否肽段从N段的第二个位置开始覆盖
#   cus_1_1 <- function(result){
#     return(result[2,] %>% .$Coverage ==1 & result %>% head(1) %>% .$Coverage !=1)
#   }
#   sapply(l, cus_1) %>% sum()
#   df$Whether_pep_start_from_N_terminal <- sapply(l, cus_1)
#   df$Whether_pep_start_from_N_terminal_2nd_pos <- sapply(l, cus_1_1)
#   df %>% filter(Whether_pep_start_from_N_terminal!=1) %>% .$ORF_id_trans %>% unique() %>% length()
#   n_distinct(df$ORF_id_trans)
#   return(list(pos_cov_sum_1,df))
# }
# all_sep_ms_pep %>% filter(Scodon=="ATG") -> df
# cus_2(df) -> res_l
# res_l[[2]] %>% .$ORF_id_trans %>% unique_n
# res_l[[2]] %>% filter(Whether_pep_start_from_N_terminal==TRUE) %>% .$ORF_id_trans %>% unique_n
# res_l[[2]] %>% filter(Whether_pep_start_from_N_terminal_2nd_pos==TRUE) %>% .$ORF_id_trans %>% unique_n
# res_l[[2]] %>% filter(Whether_pep_start_from_N_terminal!=TRUE & Whether_pep_start_from_N_terminal_2nd_pos!=TRUE) %>% .$ORF_id_trans %>% unique_n
# 
# # 查看有多少肽段从N端的第一个位置开始覆盖
# all_sep_ms_pep %>% filter(Scodon!="ATG") -> df
# cus_2(df) -> res_l
# res_l[[2]] %>% .$ORF_id_trans %>% unique_n
# res_l[[2]] %>% filter(Whether_pep_start_from_N_terminal==TRUE) %>% .$ORF_id_trans %>% unique_n

# 导出相对位置
all_sep %>% filter(Strand=="+") ->  all_sep_forward
data.frame(all_sep_forward$Chr,all_sep_forward$Start-15,all_sep_forward$Start+18,all_sep_forward$ORF_id_trans,".",
           all_sep_forward$Strand) -> window_15_start_codon_f
all_sep %>% filter(Strand=="-") ->  all_sep_reverse
data.frame(all_sep_reverse$Chr,all_sep_reverse$End-18,all_sep_reverse$End+15,all_sep_reverse$ORF_id_trans,".",
           all_sep_reverse$Strand) -> window_15_start_codon_r
colnames(window_15_start_codon_r) <- colnames(window_15_start_codon_f)
rbind(window_15_start_codon_f,window_15_start_codon_r) -> window_15_start_codon
fwrite(window_15_start_codon,"./output/S9/window_15_start_codon.bed",col.names = F,sep = '\t')


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
get_window_15_start_codon(all_sep %>% filter(Scodon=="ATG")) -> window_15_start_codon_cano
get_window_15_start_codon(all_sep %>% filter(Scodon!="ATG")) -> window_15_start_codon_uncano
fwrite(window_15_start_codon_cano,"./output/S9/window_15_start_codon_cano.bed",col.names = F,sep = '\t')
fwrite(window_15_start_codon_uncano,"./output/S9/window_15_start_codon_uncano.bed",col.names = F,sep = '\t')

# 得到fasta文件
command <- "bedtools getfasta -s -fi /home/user/data/lit/database/public/genome/hg38/hg38.fa -bed ./output/S9/window_15_start_codon_cano.bed > ./output/S9/window_15_start_codon_cano.fa"
system(command)

# 两种方法可以来绘制
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
seqlogo_cus("output/S9/window_15_start_codon.fa","output/S9/window_15_start_codon.pdf",TRUE)

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
ggseqlogo_cus("output/S9/window_15_start_codon.fa","output/S9/window_15_start_codon.pdf",FALSE)
ggseqlogo_cus("output/S9/window_15_start_codon_cano.fa","output/S9/window_15_start_codon_cano.pdf",FALSE)
ggseqlogo_cus("output/S9/window_15_start_codon_uncano.fa","output/S9/window_15_start_codon_uncano.pdf",FALSE)

##### 添加注释蛋白以及过滤掉的蛋白做正对照和负对照 #####
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S7/all_putative_sep_undetected_sampled.tab_info.txt")->sep_novel_undetected
name <- "window_15_start_codon_sep_novel_undetected"
get_window_15_start_codon(sep_novel_undetected) -> window_15_start_codon_sep_novel_undetected
output_bed <- paste0("./output/S9/",name,".bed")
output_fa <- paste0("./output/S9/",name,".fa")
output_pdf <- paste0("./output/S9/",name,".pdf")
fwrite(window_15_start_codon_sep_novel_undetected,output_bed,col.names = F,sep = '\t')
command <- paste0("bedtools getfasta -s -fi /home/user/data/lit/database/public/genome/hg38/hg38.fa -bed ",
                  output_bed," > ",output_fa)
system(command)
ggseqlogo_cus(output_fa,output_pdf,TRUE)
ggseqlogo_cus(output_fa,paste0("./output/S9/",name,".nosc.pdf"),FALSE)
seqlogo_cus(output_fa, paste0("./output/S9/",name,".seqlogo.pdf"),TRUE)

inte_func <- function(tab_info_path,output_path,name){
  fread_c(tab_info_path)->sep
  get_window_15_start_codon(sep) -> window_15_start_codon_sep
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
          "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S9/",
          "anno_sep"
          )