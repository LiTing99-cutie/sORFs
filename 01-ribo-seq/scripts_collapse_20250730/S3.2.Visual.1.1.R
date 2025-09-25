args <- commandArgs(TRUE)
# 1.加载R包
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()
# nonCano.sorf.filtered.add_meta.orf_type_annotated.txt
sORF_path <- args[1]
# nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.filtered.txt
all_sORF_path <- args[2]
output_path <- args[3]

# 1. 每个sORF被不同方法和样本独立检测到的次数
fread_c(sORF_path) -> sORF
ggplot(data=sORF,aes(x=Method_n)) +geom_bar(width = 0.8)+theme_3() -> p
ggsave(p,filename = o("sORF_method_n.bar.pdf"),height = 5,width = 5)
ggplot(data=sORF,aes(x=Sample_n))+geom_bar(width = 0.8)+theme_3() -> p
ggsave(p,filename = o("sORF_sample_n.bar.pdf"),height = 5,width = 5) 
# 2. 每个样本和方法检测到sORF的个数
fread_c(all_sORF_path) %>% count(Sample,Method) -> sample_method_orfs
bar_plot_basic(sample_method_orfs,"Sample","n",fil_col = "Method") -> p
ggsave(p,filename = o("per_sample_per_method_number.bar.pdf"),height = 5,width = 5)

# 3. 鉴定到的sORF的类型
pie_plot(sORF,"ORF_type_1") -> p
ggsave(p,filename = o("sORF_type.pie.pdf"),height = 5,width = 5)
