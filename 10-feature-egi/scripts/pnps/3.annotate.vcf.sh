source activate biotools
snpEff_dir=/home/user/data3/lit/anaconda3/envs/biotools/share/snpeff-5.2-1
myGenome=hg38_custom
out_dir=$(realpath ../processed/pnps/orf_regions)
vcf=$out_dir/orf_regions.vcf.gz
# 4) 注释
echo "开始注释VCF (这可能需要很长时间)... at $(date '+%Y-%m-%d %H:%M:%S')"
time java -Xmx10G -jar $snpEff_dir/snpEff.jar eff \
    -c $snpEff_dir/snpEff.config.$myGenome $myGenome ${vcf} > \
    ${out_dir}/${myGenome}.eff.vcf \
    -csvStats ${out_dir}/${myGenome}.csv \
    -stats ${out_dir}/${myGenome}.html && echo "snpEff done"