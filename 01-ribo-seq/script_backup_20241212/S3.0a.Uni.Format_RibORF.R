#Usage riborf.R repre.valid.pred.pvalue.parameters.txt repre.valid.ORF.genepred.txt nonCano.sorf.genepred.formatted.txt
args <- commandArgs(TRUE)
suppressMessages(library(data.table))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
fread(args[1],data.table = F) -> orfs
fread(args[2],data.table = F) -> genePred
# 筛选起始密码子为ATG；预测分数大于0.7；类型为非经典；长度为6-150
orfs <- orfs %>%
  separate(col = "orfID", into = c("ID", "Chr", "Strand", "Rank", "Span","RelaStart", "RelaEnd", "Type", "Scodon"), sep = "[|:]", remove = FALSE)
subset(orfs,Scodon=="ATG")-> orfs
subset(orfs,Type!="canonical" & pred.pvalue>=0.7) -> nonCanonical_orfs
nonCanonical_orfs %>% subset(length<=453 & length>=21) -> nonCanonical_sorfs
# 替换ID名称
new_id <- paste0(nonCanonical_sorfs$ID,nonCanonical_sorfs$Strand,nonCanonical_sorfs$Chr,
":",nonCanonical_sorfs$codon5,"-",nonCanonical_sorfs$codon3)
nonCanonical_sorfs$new_id <- new_id
merge(nonCanonical_sorfs[,c("orfID","new_id")],genePred,by.x="orfID",by.y="V1") -> m
m <- m[,-1]
# 导出文件
fwrite(m,file = args[3],sep = '\t',col.names = F)


