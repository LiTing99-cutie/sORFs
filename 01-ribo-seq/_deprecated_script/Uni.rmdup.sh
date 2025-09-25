# 在每一个Ribo-ORFs文件夹下运行
cd merge
sample=$(echo $PWD | awk -F'/' '{print $(NF-3)}')
echo $sample

rmdup_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/Uni.rmdup.R

# PRICE
cd PRICE
Rscript $rmdup_script nonCano.sorf.tab nonCano.sorf.rmDup.tab
awk -v OFS='\t' '{print $1,$2,$3,"PRICE"}' nonCano.sorf.rmDup.tab > nonCano.sorf.meta.txt

# RiboCode
cd ../RiboCode
Rscript $rmdup_script nonCano.sorf.tab nonCano.sorf.rmDup.tab
awk -v OFS='\t' '{print $1,$2,$3,"RiboCode"}' nonCano.sorf.rmDup.tab > nonCano.sorf.meta.txt

# RibORF
cd ../RibORF
Rscript $rmdup_script nonCano.sorf.tab nonCano.sorf.rmDup.tab
awk -v OFS='\t' '{print $1,$2,$3,"RibORF"}' nonCano.sorf.rmDup.tab > nonCano.sorf.meta.txt

# merge
cd ..
find ./ -name nonCano.sorf.meta.txt | xargs cat > tmp.txt
awk -v OFS='\t' '{print $1,$2,$3,$4,"'$sample'"}' tmp.txt  > nonCano.sorf.meta.merge.txt
rm -rf tmp.txt
