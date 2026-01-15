
source activate biotools
##### 准备输入的gtf和bed文件 ######
out_dir=$(realpath ../processed/pnps/20251211_fine_tune)
mkdir -p $out_dir
all_orfs_gtf=/home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/all_ORFs.gt30.gtf
all_orfs_bed=/home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/all_ORFs.gt30.ORF_no_intron.bed
denovo_gtf=/home/user/data3/lit/project/sORFs/11-denovo-list/processed/new_list_20251128/denovo_list.gtf
denovo_CDS_bed=/home/user/data3/lit/project/sORFs/11-denovo-list/processed/new_list_20251128/denovo_list.bed
cat $all_orfs_gtf $denovo_gtf > $out_dir/gene.gtf
cat $all_orfs_bed $denovo_CDS_bed > $out_dir/gene.bed
custom_gtf=$out_dir/gene.gtf
CDS_bed=$out_dir/gene.bed
myGenome=hg38_denovo_20251211

snpEff_dir=/home/user/data3/lit/anaconda3/envs/biotools/share/snpeff-5.2-1
fa=/home/user/data/lit/database/public/genome/hg38/hg38.fa
vcf=/home/user/data3/lit/project/sORFs/11-denovo-list/public_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
##### 第一步、构建注释数据库 ######
# 1) 目录与文件
echo "准备snpEff数据库... at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $snpEff_dir/data/$myGenome
mkdir -p $snpEff_dir/data/genomes
cp $fa $snpEff_dir/data/genomes/$myGenome.fa
cp $fa $snpEff_dir/data/$myGenome/sequences.fa
cp $custom_gtf $snpEff_dir/data/$myGenome/genes.gtf

# 2) 配置
cp $snpEff_dir/snpEff.config $snpEff_dir/snpEff.config.$myGenome
printf "data.dir = $snpEff_dir/data\n$myGenome.genome : $myGenome\n" >> $snpEff_dir/snpEff.config.$myGenome

# 3) 构建
echo "构建snpEff数据库... at $(date '+%Y-%m-%d %H:%M:%S')"
time java -Xmx8g -jar $snpEff_dir/snpEff.jar build -c $snpEff_dir/snpEff.config.$myGenome -noCheckProtein -noCheckCds -gtf22 -v $myGenome && echo "done"

##### 第二步、提取对应区域的vcf ######
# 输入文件
out_dir_1=$(realpath $out_dir/orf_regions)
mkdir -p $out_dir_1

# 检查BED文件格式和染色体命名
echo "Checking BED file..."
head -3 $CDS_bed
echo ""

# 检查VCF染色体命名
echo "Checking VCF chromosomes..."
bcftools view -H $vcf | head -1000 | cut -f1 | sort -u | head -5
echo ""

# 一步到位
# -R 只保留BED文件中指定区域的variants
# -z 输出格式: z=bgzip压缩VCF, v=未压缩, b=BCF
bcftools view -R $CDS_bed $vcf -O z -o $out_dir_1/orf_regions.vcf.gz && \
bcftools index -c $out_dir_1/orf_regions.vcf.gz && \
echo "Extracted $(bcftools view -H $out_dir_1/orf_regions.vcf.gz | wc -l) variants"

##### 第三步、注释对应区域的vcf ######
vcf=$out_dir_1/orf_regions.vcf.gz
echo "开始注释VCF... at $(date '+%Y-%m-%d %H:%M:%S')"
time java -Xmx10G -jar $snpEff_dir/snpEff.jar eff \
    -c $snpEff_dir/snpEff.config.$myGenome $myGenome ${vcf} > \
    ${out_dir_1}/${myGenome}.eff.vcf \
    -csvStats ${out_dir_1}/${myGenome}.csv \
    -stats ${out_dir_1}/${myGenome}.html && echo "snpEff done"