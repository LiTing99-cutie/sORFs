args <- commandArgs(TRUE)

suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

read.table(args[1]) -> nonCano.sorf.tab

# 去除同一位置编码出的同一seq
rm_dup <- function(nonCano.sorf.tab){
  nonCano.sorf.tab -> df
  colnames(df) <- c("ORF_id","Seq")
  # 同一个位置编码出的同一个seq只保留一个（包括overlapped基因或者同一基因的不同转录本）
  df$Location <- str_extract(df$ORF_id,"[+-]chr\\w+:\\d+[+-]\\d+")
  # ID代表独特的小肽
  df$ID <- paste0(df$Location,":",df$Seq)
  df %>% distinct(ID,.keep_all = T) -> df_rm_dup
  return(df_rm_dup)
}

fwrite(rm_dup(nonCano.sorf.tab)[,c("ORF_id","Seq","ID")],args[2],col.names =F,sep = '\t')
