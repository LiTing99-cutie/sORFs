source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/")
##### 肽段在SEP上相对位置的分布【可视化】 ##### 
fread_c("./S4/novel_sep_ms_pep.txt")-> all_sep_ms_pep
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
# N端开始的位置index以及覆盖度
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
# 创建输出路径
create_path("./S4/plot")
plot_cus(all_sep_ms_pep %>% filter(Scodon!="ATG")) %>% ggsave(.,filename = "./S4/plot/noncano.peptide.cov.pdf",width = 6,height = 4)
plot_cus(all_sep_ms_pep %>% filter(Scodon=="ATG")) %>% ggsave(.,filename = "./S4/plot/cano.peptide.cov.pdf",width = 6,height = 4)
plot_cus(all_sep_ms_pep) %>% ggsave(.,filename = "./S4/plot/all.peptide.cov.pdf",width = 6,height = 4)
