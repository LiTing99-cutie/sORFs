library(dplyr)
library(magrittr)
library(stringr)
library(data.table)

setwd("/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163093/output/")

read.table("./Ribo-ORFs/merge/RiboCode/nonCano.sorf.rmDup.tab") -> RiboCode.nonCano.sorf.rmDup.tab
read.table("./Ribo-ORFs/merge/PRICE/nonCano.sorf.rmDup.tab") -> PRICE.nonCano.sorf.rmDup.tab
read.table("./Ribo-ORFs/merge/RibORF/nonCano.sorf.rmDup.tab") -> RibORF.nonCano.sorf.rmDup.tab

merge(data.frame("ID"=RiboCode.nonCano.sorf.rmDup.tab$V3,"RiboCode"=1),
      data.frame("ID"=PRICE.nonCano.sorf.rmDup.tab$V3,"PRICE"=1),by="ID",all=T) %>% 
  merge(data.frame("ID"=RibORF.nonCano.sorf.rmDup.tab$V3,"RibORF"=1),by="ID",all=T) -> merge_res
merge_res[is.na(merge_res)] <- 0
merge_res <- merge_res[,-1]
pdf(file = "./Ribo-ORFs/merge/upset.plot.pdf")
UpSetR::upset(merge_res,sets = colnames(merge_res))
dev.off()

library(VennDiagram)
venn.diagram(x=list(RiboCode=RiboCode.nonCano.sorf.rmDup.tab$V3,
                    PRICE=PRICE.nonCano.sorf.rmDup.tab$V3,
                    RibORF=RibORF.nonCano.sorf.rmDup.tab$V3),
             filename = NULL,
             ) -> venn.plot
pdf(file = "./Ribo-ORFs/merge/venn.pdf")
grid.draw(venn.plot)
dev.off()
