setwd("/home/user/data3/lit/project/sORFs/01-ribo-seq/output/assembled_trans")

source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()

fread("./stringtie_output/mouse_brain_new_trans.gtf",
        skip = 2,data.table = F) -> new_gtf

# 提取 gene_id
new_gtf$gene_id <- sub('.*gene_id "([^"]+)";.*', '\\1', new_gtf$V9)

# 提取 transcript_id
new_gtf$transcript_id <- sub('.*transcript_id "([^"]+)";.*', '\\1', new_gtf$V9)

# 构建gene feature
filter(new_gtf,V3=="transcript") %>% group_by(gene_id) %>% 
  mutate(V4=min(V4),V5=max(V5),V3="gene") %>% distinct(gene_id,.keep_all = T) -> new_gtf_gene
n_distinct(new_gtf$gene_id)
new_gtf_gene$V9 <- paste0("gene_id"," ","\"",new_gtf_gene$gene_id,"\"","; ",
                          "gene_type"," ","\"","newly_assembled","\"","; ",
       "gene_name"," ","\"",new_gtf_gene$gene_id,"\"",";")

# 构建transcript feature
filter(new_gtf,V3=="transcript") -> new_gtf_transcript
new_gtf_transcript$V9 <- paste0("gene_id"," ","\"",new_gtf_transcript$gene_id,"\"","; ",
                     "transcript_id"," ","\"",new_gtf_transcript$transcript_id,"\"","; ",
                     "gene_type"," ","\"","newly_assembled","\"","; ",
                          "gene_name"," ","\"",new_gtf_transcript$gene_id,"\"","; ",
                     "transcript_type"," ","\"","newly_assembled","\"","; ",
                     "transcript_name"," ","\"",new_gtf_transcript$transcript_id,"\"",";")
# 构建exon feature
filter(new_gtf,V3=="exon") -> new_gtf_exon
new_gtf_exon$exon_number <- sub('.*exon_number "([^"]+)";.*', '\\1', new_gtf_exon$V9)
new_gtf_exon$V9 <- paste0("gene_id"," ","\"",new_gtf_exon$gene_id,"\"","; ",
                                "transcript_id"," ","\"",new_gtf_exon$transcript_id,"\"","; ",
                                "gene_type"," ","\"","newly_assembled","\"","; ",
                                "gene_name"," ","\"",new_gtf_exon$gene_id,"\"","; ",
                                "transcript_type"," ","\"","newly_assembled","\"","; ",
                                "transcript_name"," ","\"",new_gtf_exon$transcript_id,"\"","; ",
                                "exon_number"," ","\"",new_gtf_exon$exon_number,"\"",";")
new_gtf_exon$exon_number <- NULL
# 合并
rbind(new_gtf_gene,new_gtf_transcript,new_gtf_exon) -> new_gtf_reformat
new_gtf_reformat$transcript_order <- str_split_fixed(new_gtf_reformat$transcript_id,"\\.",3)[,3] %>% 
  as.numeric()
new_gtf_reformat$V3 <- factor(new_gtf_reformat$V3,levels = c("gene","transcript","exon"))
new_gtf_reformat %>% filter(V7=="+") %>% 
  arrange(gene_id,transcript_order,V3,V4) -> tmp_1
new_gtf_reformat %>% filter(V7=="-") %>% 
  arrange(gene_id,transcript_order,V3,desc(V4)) -> tmp_2

rbind(tmp_1,tmp_2) -> new_gtf_reformat_ordered

write.table(new_gtf_reformat_ordered[,1:9],file = "./stringtie_output/mouse_brain_new_trans_reformat_ordered.gtf",
            sep='\t',col.names =F,quote = F,row.names = F)
