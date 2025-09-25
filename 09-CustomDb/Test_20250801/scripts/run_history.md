# 需要根据iso-seq输出生成gtf以及翻译出来的蛋白质序列
mkdir -p ../processed/mkAnno_for_moPepGen
outDir=../processed/mkAnno_for_moPepGen
isoseq_out=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/processed/classify/collapsed_classification.filtered_lite_classification.txt
isoseq_out_gff=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/processed/collapsed.sorted.filtered_lite.gff
less $isoseq_out|grep full-splice_match|cut -f8 > $outDir/fsm.transcript.id.txt
gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
<!-- grep -f $outDir/fsm.transcript.id.txt $gtf > $outDir/fsm.transcript.id.gtf -->
## 加速
LC_ALL=C parallel -j50 --pipepart --block 20M -a "$gtf" -k grep -F -f "$outDir/fsm.transcript.id.txt" > "$outDir/fsm.transcript.id.gtf"
less $isoseq_out|egrep -v "full-splice_match|incomplete-splice_match"|cut -f1|tail -n +2> $outDir/novel.pb.id.txt
LC_ALL=C parallel -j50 --pipepart --block 20M -a "$isoseq_out_gff" -k grep -F -f "$outDir/novel.pb.id.txt" > "$outDir/novel.pb.id.gff"
gffread $outDir/novel.pb.id.gff -T -o $outDir/novel.pb.id.gtf
cat "$outDir/fsm.transcript.id.gtf" $outDir/novel.pb.id.gtf > $outDir/custom.gtf
gtfToGenePred -geneNameAsName2 -genePredExt $outDir/custom.gtf $outDir/custom.15.gpe
gtfToGenePred $outDir/custom.gtf $outDir/custom.10.gpe

pc_translations_fa=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.pc_translations.fa.gz
seqkit seq -w 0 $pc_translations_fa > /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.pc_translations.w0.fa
pc_translations_fa_1=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.pc_translations.w0.fa
LC_ALL=C parallel -j50 --pipepart --block 20M -a "$pc_translations_fa_1" -k grep -A 1 -F -f "$outDir/fsm.transcript.id.txt" > "$outDir/fsm.pc_translations.fa"
grep -v ^- "$outDir/fsm.pc_translations.fa" > "$outDir/fsm.pc_translations.1.fa" 
[ -d $outDir/index ] && rm -rf $outDir/index
moPepGen generateIndex \
    --genome-fasta /home/user/data/lit/database/public/genome/hg38/hg38.fa \
    --annotation-gtf $outDir/custom.gtf \
    --proteome-fasta $outDir/fsm.pc_translations.1.fa \
    --output-dir $outDir/index
vep_txt=/home/user/data3/lit/project/sORFs/07-Genome/results/vep/human_brain_21pcw_vep_annotated.txt
<!-- moPepGen parseVEP \
    -i $vep_txt \
    --index-dir $outDir/index \
    --output-path $outDir/parsed_vep -->
mkdir $outDir/parsed_vep/
moPepGen parseVEP \
  -i $vep_txt \
  --index-dir $outDir/index \
  --output-path $outDir/parsed_vep/human_brain_21pcw.gvf \
  --source gSNP,gIndel
