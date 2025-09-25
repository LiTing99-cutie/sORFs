# 注释每一个ORF的类型

args <- commandArgs(TRUE)
# 1.加载R包
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
# nonCano.sorf.filtered.add_meta.txt
sorf_path <- args[1]
output_file <- args[2]
# 20240212修改为可变参数
gpe_path <- args[3]

# 2.添加CDS_start以及CDS_end信息
sORF <- fread(sorf_path)
gpe <- fread(gpe_path,header = F,data.table = F)
gpe[,c("V1","V6","V7")] -> gpe_selected
colnames(gpe_selected) <- c("ENS_id","CDS_start","CDS_end")
merge(sORF,gpe_selected) -> sORF

# 3.对ORF相对于经典CDS的位置进行注释
# 定义函数
## Other_CDS_variant是起始密码子与经典CDS的起始密码子重合，但是终止密码子在经典CDS的终止密码子上游
# 20240212 增加了canonical类型
get_orf_type <- function(Trans_type,Strand,Start,End,CDS_start,CDS_end){
  if(Trans_type!="protein_coding"){
    return("ncORF")}
  else {
    if (Strand == "+"){
      if (Start == CDS_start && End == CDS_end){
        return("Canonical")
      }else if (End <= CDS_start){
        return("uORF")
      }else if (Start < CDS_start && End > CDS_start && End < CDS_end){
        return("uoORF")
      }else if (Start < CDS_start && End == CDS_end){
        return("Extension")
      }else if (Start == CDS_start && End > CDS_end){
        return("Readthrough")
      }else if (Start < CDS_start && End > CDS_end){
        return("Containing")
      }else if (Start == CDS_start){
        return("Other_CDS_variant")
      }else if (Start >= CDS_end){
        return("dORF")
      }else if (End > CDS_end && Start > CDS_start && Start < CDS_end){
        return("doORF")
      }else if (Start > CDS_start && End < CDS_end){
        return("Internal")
      }else if (Start > CDS_start && End == CDS_end){
        return("Truncation")
      }
    }
    else if (Strand == "-"){
      if (Start == CDS_start && End == CDS_end){
        return("Canonical")
      }else if (Start >= CDS_end){
        return("uORF")
      }else if (End > CDS_end && Start > CDS_start && Start < CDS_end){
        return("uoORF")
      }else if (End > CDS_end && Start == CDS_start){
        return("Extension")
      }else if (End == CDS_end && Start < CDS_start){
        return("Readthrough")
      }else if (End > CDS_end && Start < CDS_start){
        return("Containing")
      }else if (End == CDS_end){
        return("Other_CDS_variant")
      }else if (End <= CDS_start){
        return("dORF")
      }else if (Start < CDS_start && End > CDS_start && End < CDS_end){
        return("doORF")
      }else if (Start > CDS_start && End < CDS_end){
        return("Internal")
      }else if (Start == CDS_start && End < CDS_end){
        return("Truncation")
      }
    }
  }
  return("other")
}
# 生成ORF_type
mapply(get_orf_type,sORF$Transcript_type,sORF$Strand,sORF$Start,sORF$End,sORF$CDS_start,sORF$CDS_end) -> ORF_type_reassign
# 添加回去
sORF$ORF_type <- ORF_type_reassign
# 进一步合并，并根据转录本的类型生成子类
## 子类1
sORF %>% mutate(ORF_type_1=case_when(
ORF_type=="Truncation" | ORF_type=="Extension" | ORF_type=="Readthrough" | ORF_type=="Other_CDS_variant" ~ "CDS_variant",
TRUE ~ ORF_type)) -> sORF
## 子类2
sORF %>% mutate(ORF_type_2=case_when(
grepl("pseudogene",sORF$Transcript_type) ~ "ncORF_pseudogene",
Transcript_type == "lncRNA" | Transcript_type == "retained_intron" | Transcript_type == "processed_transcript" ~ "ncORF_lncRNA",
TRUE ~ ORF_type)) -> sORF
sORF %>% mutate(ORF_type_2=case_when(
ORF_type_2=="ncORF" ~ "ncORF_other",
TRUE ~ ORF_type_2)) -> sORF
## 子类3
sORF %>% mutate(ORF_type_3=case_when(
ORF_type_2=="ncORF_lncRNA" ~ Transcript_type,
TRUE ~ ORF_type_2)) -> sORF
sORF %>% mutate(ORF_type=NULL) -> sORF

# 4.导出
fwrite(sORF,file = output_file,sep = '\t')
