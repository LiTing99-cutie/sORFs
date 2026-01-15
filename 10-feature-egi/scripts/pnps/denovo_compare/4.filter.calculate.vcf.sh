source activate biotools
out_dir=$(realpath ../processed/pnps/20251211_fine_tune)
mkdir -p $out_dir
eff_vcf=$out_dir/orf_regions/hg38_denovo_20251211.eff.vcf 
bcftools view $eff_vcf -Oz -o ${eff_vcf}.gz
bcftools index -c ${eff_vcf}.gz
# 188697
bcftools index -n ${eff_vcf}.gz
# way 1
mkdir $out_dir/filtered_1
less ${eff_vcf} |grep -v "^#" | grep ";DB;" | grep POSITIVE_TRAIN_SITE > $out_dir/filtered_1/eff.vcf
## 55137
wc -l $out_dir/filtered_1/eff.vcf
# way 2
mkdir $out_dir/filtered_2
bcftools view -H -f PASS ${eff_vcf}.gz > $out_dir/filtered_2/eff.vcf
## 139669
wc -l $out_dir/filtered_2/eff.vcf
# way 3
mkdir $out_dir/filtered_3
# bcftools view -i 'AF >= 0.05' ${eff_vcf}.gz -o $out_dir/filtered_3/eff.vcf.gz
bcftools view -i 'AF >= 0.05' ${eff_vcf}.gz -o $out_dir/filtered_3/eff.vcf
## 27021
less $out_dir/filtered_3/eff.vcf |grep -v "^#"|wc -l

# vcf=$1
# out_dir=$2
# CDS_bed=$3
awk 'BEGIN{OFS="\t"} {sub(/_uq[0-9]+$/, "", $4); print}' $out_dir/gene.bed > $out_dir/gene.rmSuffix.bed
bash pnps/denovo_compare/Uni.calculate.pn.ps.sh $out_dir/filtered_1/eff.vcf  $out_dir/filtered_1 $out_dir/gene.rmSuffix.bed
bash pnps/denovo_compare/Uni.calculate.pn.ps.sh $out_dir/filtered_2/eff.vcf  $out_dir/filtered_2 $out_dir/gene.rmSuffix.bed
bash pnps/denovo_compare/Uni.calculate.pn.ps.sh $out_dir/filtered_3/eff.vcf  $out_dir/filtered_3 $out_dir/gene.rmSuffix.bed
mkdir $out_dir/filtered_0
bash pnps/denovo_compare/Uni.calculate.pn.ps.sh $eff_vcf  $out_dir/filtered_0 $out_dir/gene.rmSuffix.bed
