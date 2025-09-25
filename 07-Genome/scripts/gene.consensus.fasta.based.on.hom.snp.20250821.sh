ref=/home/user/data/lit/database/public/genome/hg38/hg38.fa
vcf=/home/user/data3/lit/project/sORFs/07-Genome/processed/vcf_filter/human_brain_21pcw_snps_filtered_pass.vcf.gz
mkdir -p ../results/custom_fa
out_dir=../results/custom_fa
bcftools view -i 'GT="1/1"' $vcf -Oz -o $out_dir/hom.vcf.gz
bcftools index $out_dir/hom.vcf.gz
cat $ref | bcftools consensus $out_dir/hom.vcf.gz > $out_dir/custom_ref.fa