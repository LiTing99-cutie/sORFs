# 1.1 注释vcf
nohup bash pnps/1.build.annotation.sh > ../log/1.build.annotation.log 2>&1 &
# nohup bash evolution_para/1.annotate_vcf.test.sh > ../log/annotate_vcf.test.log 2>&1 &
# tail -f ../log/annotate_vcf.log
# 1.2 提取目标区域的vcf
# vcf=/home/user/data3/lit/project/sORFs/11-denovo-list/public_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
# chr1    99335318        99335387        99309590-99464377_203   .       +
nohup bash pnps/2.extract.vcf.sh &> ../log/extract.vcf.log &
nohup bash pnps/3.annotate.vcf.sh &> ../log/annotate.vcf.log &

# 1.3 计算目标区域的同义位点和非同义位点数量
# vcf=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/pnps/orf_regions/hg38_custom.eff.vcf
# orf_bed=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/intergenic.lnc.noncano.cano.orfs.bed
bash pnps/4.intersect.orf.bed.vcf.sh

### 20251129 计算de novo基因的pn/ps ### 
nohup bash pnps/1-4.combine.sh &> ../log/1-4.combine.log &

### 20251129 计算背景基因的p-sites和n-sites ### 
# test
cds_fa=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/all_orfs.cds.fa
faTrans -stop $cds_fa ../processed/seve_list_compare/input_for_pnps/all_orfs.cds.noM.fa
prot_fa=../processed/seve_list_compare/input_for_pnps/all_orfs.cds.noM.fa
out_dir=../processed/pnps/orfs_fa
mkdir -p $out_dir
orf_id="PB.18548.1:chr17:+|48|1683:845:959|noncoding|TTG"
seqkit grep -p "$orf_id" $prot_fa > $out_dir/"$orf_id".fa
cat $out_dir/"$orf_id".fa $out_dir/"$orf_id".fa > $out_dir/"$orf_id".ortho.prot.fa
seqkit replace -p '.+' -r 'human_{nr}' $out_dir/"$orf_id".ortho.prot.fa > $out_dir/"$orf_id".ortho.prot.rename.fa
seqkit grep -p "$orf_id" $cds_fa > $out_dir/"$orf_id".cds.fa
cat $out_dir/"$orf_id".cds.fa $out_dir/"$orf_id".cds.fa > $out_dir/"$orf_id".ortho.cds.fa
seqkit replace -p '.+' -r 'human_{nr}' $out_dir/"$orf_id".ortho.cds.fa > $out_dir/"$orf_id".ortho.cds.rename.fa
script=/home/user/data3/lit/project/sORFs/10-feature-egi/scripts/kaks/kaks.bin.sh
mkdir -p $out_dir/"$orf_id"
bash $script -p $out_dir/"$orf_id".ortho.prot.rename.fa -c $out_dir/"$orf_id".ortho.cds.rename.fa -o $out_dir/"$orf_id" -m all -k NG
awk 'NR>1{print "'$orf_id'""\t"$0}' $out_dir/"$orf_id"/kaks.res

# run
# nohup bash pnps/get.background.sh &> ../log/get.background.log &
nohup bash pnps/get.background.parallel.sh &> ../log/get.background.parallel.log &
nohup bash pnps/get.background.parallel.GMYN.sh &> ../log/get.background.parallel.GMYN.log &