# 在每一个Ribo-ORFs文件夹下运行
mkdir -p merge/RibORF merge/PRICE merge/RiboCode
cd merge
sample=$(echo $PWD | awk -F'/' '{print $(NF-3)}')
echo $sample

genePred=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.genePred.txt
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa

# 使用base中的python
source activate base
translate_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/Uni.translate_gtf.sh
rmdup_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/Uni.rmdup.R

# PRICE
cd PRICE
price_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/Generate_genepred_PRICE.py
## 构建一个gpe
python $price_script ../../PRICE/*.orfs.tsv $genePred nonCano.formatted.gpe
## 根据gpe得到fasta文件
genePredToGtf file nonCano.formatted.gpe nonCano.formatted.gtf
bash $translate_script nonCano.formatted.gtf $fa
mv prot.fa nonCano.fa
## 得到长度后根据长度过滤
seqkit fx2tab -l nonCano.fa | awk -v OFS='\t' '{print $1,$3}' > nonCano.pro.l.txt
awk -v OFS='\t' '$2>=6 && $2<=150' nonCano.pro.l.txt > nonCano.sorf.pro.l.txt
seqkit grep -n -f <(cut -f 1 nonCano.sorf.pro.l.txt) nonCano.fa > nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
## 去掉重复位置编码出的fasta
Rscript $rmdup_script nonCano.sorf.tab nonCano.sorf.rmDup.tab
awk -v OFS='\t' '{print $1,$2,$3,"PRICE"}' nonCano.sorf.rmDup.tab > nonCano.sorf.meta.txt

# RiboCode
cd ../RiboCode
compare_script_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/ref/Ribo-seq-Tool-Comparison-Scripts-v2.0/Scripts_for_RiboCode_Analysis
output_1=../../RiboCode/$sample.txt
output_2=../../RiboCode/$sample.gtf
## 过滤gtf
cat <(head -n1 $output_1) <(awk '$2 != "annotated" && $10 <= 450 ' $output_1) > nonCano.sorf.txt
tail -n +2 nonCano.sorf.txt | cut -f 1 > nonCano.sorf.id.txt
grep -F -f nonCano.sorf.id.txt $output_2 > nonCano.sorf.gtf
## 统一ORF_id
python $compare_script_path/Formatting_RiboCode_gtf.py nonCano.sorf.gtf nonCano.sorf.formatted.gtf
## 得到fasta文件
python $compare_script_path/Generate_Fasta_RiboCode.py nonCano.sorf.txt nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
## 去掉重复位置编码出的fasta
Rscript $rmdup_script nonCano.sorf.tab nonCano.sorf.rmDup.tab
awk -v OFS='\t' '{print $1,$2,$3,"RiboCode"}' nonCano.sorf.rmDup.tab > nonCano.sorf.meta.txt

# RibORF
cd ../RibORF
riborf_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/Format_RibORF.R
output_1=../../RibORF/repre.valid.pred.pvalue.parameters.txt
output_2=../../RibORF/repre.valid.ORF.genepred.txt
Rscript $riborf_script "$output_1" "$output_2" "nonCano.sorf.formatted.gpe"
genePredToGtf file nonCano.sorf.formatted.gpe nonCano.sorf.formatted.gtf
bash $translate_script nonCano.sorf.formatted.gtf $fa
mv prot.fa nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
## 去掉重复位置编码出的fasta
Rscript $rmdup_script nonCano.sorf.tab nonCano.sorf.rmDup.tab
awk -v OFS='\t' '{print $1,$2,$3,"RibORF"}' nonCano.sorf.rmDup.tab > nonCano.sorf.meta.txt

# merge
cd ..
find ./ -name nonCano.sorf.meta.txt | xargs cat > tmp.txt
awk -v OFS='\t' '{print $1,$2,$3,$4,"'$sample'"}' tmp.txt  > nonCano.sorf.meta.merge.txt
rm -rf tmp.txt
