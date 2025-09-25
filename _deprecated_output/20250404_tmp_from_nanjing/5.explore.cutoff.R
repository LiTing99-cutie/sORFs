fread_c("./sep_add_overlap_pro.txt") -> sep_add_overlap_pro
sep_add_overlap_pro[is.na(sep_add_overlap_pro)] <- 0
fread_c("./sep_add_overlap_pro_filter.txt") -> sep_add_overlap_pro_filter
fread_c("./ts.meta.txt") -> ts_meta

sep_add_overlap_pro %>% filter(ORF_type_1=="ncORF") %>% count(Gene_type)
sep_add_overlap_pro %>% filter(ORF_type_1=="ncORF") %>% filter(Gene_type=="protein_coding") -> tmp
# 342
nrow(tmp)
# 35
tmp$ORF_id_trans %in% sep_add_overlap_pro_filter$ORF_id_trans %>% sum()

tmp[tmp$ORF_id_trans %in% sep_add_overlap_pro_filter$ORF_id_trans,]

# 南京服务器没有修好，粗略重新assign一下类别
## ncORF，但是基因protein-coding，过滤掉那些没有特异性肽段的
## ncORF，但是基因是ncORF，不进一步过滤
sep_add_overlap_pro %>% filter(ORF_type_1=="ncORF") %>% filter(Gene_type!="protein_coding")
## CDS_variant，过滤掉那些没有特异性肽段的
## internal，过滤掉那些没有特异性肽段的
## dORF,doORF,uORF,uoORF，不进一步过滤

table(sep_add_overlap_pro$ORF_type_1) 
bar_plot(sep_add_overlap_pro,"ORF_type_1")
sep_add_overlap_pro %>% filter(Overlapped_pro<=0.9) %>% nrow()
sep_add_overlap_pro %>% filter(Overlapped_pro==0) %>% nrow()
table(sep_add_overlap_pro_filter$ORF_type_1)
table(sep_add_overlap_pro %>% filter(Overlapped_pro!=1) %>% .$ORF_type_1)
sep_add_overlap_pro %>% filter(ORF_type_1=="ncORF") %>% filter(Gene_type=="protein_coding") %>% nrow()

sum(sep_add_overlap_pro$Ribo_evidence)
nrow(sep_add_overlap_pro)
