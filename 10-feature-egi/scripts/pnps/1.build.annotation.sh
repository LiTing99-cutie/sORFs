# 注释vcf
source activate biotools
snpEff_dir=/home/user/data3/lit/anaconda3/envs/biotools/share/snpeff-5.2-1
myGenome=hg38_custom
fa=/home/user/data/lit/database/public/genome/hg38/hg38.fa
# gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
custom_gtf=../processed/seve_list_compare/input_for_pnps/all_orfs.gtf
out_dir=../processed/pnps/
vcf=/home/user/data3/lit/project/sORFs/11-denovo-list/public_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
mkdir -p $out_dir
echo "合并GTF文件... at $(date '+%Y-%m-%d %H:%M:%S')"
cat $custom_gtf > $out_dir/merged.gtf
# 1) 目录与文件
echo "准备snpEff数据库... at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $snpEff_dir/data/$myGenome
mkdir -p $snpEff_dir/data/genomes
cp $fa $snpEff_dir/data/genomes/$myGenome.fa
cp $fa $snpEff_dir/data/$myGenome/sequences.fa
cp $out_dir/merged.gtf $snpEff_dir/data/$myGenome/genes.gtf

# 2) 配置
cp $snpEff_dir/snpEff.config $snpEff_dir/snpEff.config.$myGenome
printf "data.dir = $snpEff_dir/data\n$myGenome.genome : $myGenome\n" >> $snpEff_dir/snpEff.config.$myGenome

# 3) 构建
echo "构建snpEff数据库... at $(date '+%Y-%m-%d %H:%M:%S')"
time java -Xmx8g -jar $snpEff_dir/snpEff.jar build -c $snpEff_dir/snpEff.config.$myGenome -noCheckProtein -noCheckCds -gtf22 -v $myGenome && echo "done"

# 4) 注释
# echo "开始注释VCF (这可能需要很长时间)... at $(date '+%Y-%m-%d %H:%M:%S')"
# time java -Xmx10G -jar $snpEff_dir/snpEff.jar eff \
#     -c $snpEff_dir/snpEff.config.$myGenome $myGenome ${vcf} > \
#     ${out_dir}/${myGenome}.eff.vcf \
#     -csvStats ${out_dir}/${myGenome}.csv \
#     -stats ${out_dir}/${myGenome}.html && echo "snpEff done"



