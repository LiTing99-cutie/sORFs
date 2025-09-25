#20241212修改脚本让一些参数可变
#Usage riborf.R repre.valid.pred.pvalue.parameters.txt repre.valid.ORF.genepred.txt stat.cutoff.txt custom 0 nonCano.sorf.genepred.formatted.txt
args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

orfs_path <- args[1]
gpe_path <- args[2]
stat_cutoff_path <- args[3]
# can be custom or youden; default custom
score_cutoff <- args[4]
only_canonical_start_codon <- as.numeric(args[5])
output_file <- args[6]

fread_c(orfs_path) -> orfs
fread_c(gpe_path) -> genePred
# 类型为非经典；长度为6-150
orfs <- orfs %>%
  separate(col = "orfID", into = c("ID", "Chr", "Strand", "Rank", "Span","RelaStart", "RelaEnd", "Type", "Scodon"), sep = "[|:]", remove = FALSE)
subset(orfs,Type!="canonical" & length<=453 & length>=21) -> nonCanonical_sorfs
# 是否只选取经典起始密码子
if(only_canonical_start_codon){
  subset(nonCanonical_sorfs,Scodon=="ATG")-> nonCanonical_sorfs
}
# 选择预测分数的cutoff
fread_c(stat_cutoff_path) -> A
if(score_cutoff=="custom"){
  filter(A,False.pos.rate<=0.1 & True.pos.rate>=0.9) %>% filter(True.pos.rate==max(True.pos.rate)) %>% filter(False.pos.rate==min(False.pos.rate)) -> B
  cutoff <- B$cutoff
}else if(score_cutoff=="youden"){
  fpr <- A[,6]
  tpr <- A[,7]
  which.max(tpr-fpr) -> idx
  cutoff <- A$cutoff[idx]
}
cat("cutoff is:",cutoff,"\n")
subset(nonCanonical_sorfs,pred.pvalue>=cutoff) -> nonCanonical_sorfs

# 替换ID名称
new_id <- paste0(nonCanonical_sorfs$ID,nonCanonical_sorfs$Strand,nonCanonical_sorfs$Chr,
                 ":",nonCanonical_sorfs$codon5,"-",nonCanonical_sorfs$codon3)
nonCanonical_sorfs$new_id <- new_id
merge(nonCanonical_sorfs[,c("orfID","new_id")],genePred,by.x="orfID",by.y="V1") -> m
m <- m[,-1]
# 导出文件
fwrite(m,file = output_file,sep = '\t',col.names = F)


