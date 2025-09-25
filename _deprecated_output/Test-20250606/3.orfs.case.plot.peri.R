source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()

bedgraph_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/test_output_20250606/target.bedtools.cov.add.pos.frame.bedgraph"
fread_c(bedgraph_path) -> bedgraph
frame_1_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/test_output_20250606/target.frame1.txt"
fread_c(frame_1_path) -> frame_1

colnames(bedgraph) <- c("Chr","Start","End","ORF_id_trans","PH_1","Strand","PH_2","Count","Pos","Frame")

bedgraph$Frame <- factor(bedgraph$Frame,levels = 1:3)

bedgraph %>% group_by(ORF_id_trans,Strand) %>% summarise(count_sum=sum(Count)) %>% arrange(desc(count_sum)) %>% head(5)

head(frame_1)

filter(frame_1,Total_Count>0) %>% mutate(Frame1_Count/Total_Count) -> peri
mean(peri$`Frame1_Count/Total_Count`)
hist(peri$`Frame1_Count/Total_Count`)

nrow(frame_1)
nrow(peri)
###### case plot #####
case_1_forward <- "ENST00000505232.5+chr4:41257653-41263240"
case_1_reverse <- "ENST00000455785.7-chr1:25901015-25901084"
case_plot <- function(case){
  filter(bedgraph,ORF_id_trans==case) -> df
  ggbarplot(df,x = "Pos",y="Count",color="Frame",fill = "Frame")+
    # 参考2022-NN配色，深棕色，深黄色，浅黄色
    scale_fill_manual(values = c("#772128","#eb944c","#faedaf"))+
    scale_color_manual(values = c("#772128","#eb944c","#faedaf"))+
    labs(x="Relative Position",y="Frequency",title = case) -> p
  return(p)
}
case_plot(case_1_forward)
case_plot(case_1_reverse)
###### peri #####
filter(peri,Group==case_1_forward)[,4]
filter(peri,Group==case_1_reverse)[,4]
