
source activate biotools
custom_gtf=/home/user/data3/lit/project/sORFs/11-denovo-list/processed/new_list_20251128/denovo_list.gtf
CDS_bed=/home/user/data3/lit/project/sORFs/11-denovo-list/processed/new_list_20251128/denovo_list.bed
myGenome=hg38_denovo
out_dir=$(realpath ../processed/pnps_denovo)
mkdir -p $out_dir

snpEff_dir=/home/user/data3/lit/anaconda3/envs/biotools/share/snpeff-5.2-1
fa=/home/user/data/lit/database/public/genome/hg38/hg38.fa
vcf=/home/user/data3/lit/project/sORFs/11-denovo-list/public_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
##### 第一步、构建注释数据库 ######
# 1) 目录与文件
# echo "准备snpEff数据库... at $(date '+%Y-%m-%d %H:%M:%S')"
# mkdir -p $snpEff_dir/data/$myGenome
# mkdir -p $snpEff_dir/data/genomes
# cp $fa $snpEff_dir/data/genomes/$myGenome.fa
# cp $fa $snpEff_dir/data/$myGenome/sequences.fa
# cp $custom_gtf $snpEff_dir/data/$myGenome/genes.gtf

# # 2) 配置
# cp $snpEff_dir/snpEff.config $snpEff_dir/snpEff.config.$myGenome
# printf "data.dir = $snpEff_dir/data\n$myGenome.genome : $myGenome\n" >> $snpEff_dir/snpEff.config.$myGenome

# # 3) 构建
# echo "构建snpEff数据库... at $(date '+%Y-%m-%d %H:%M:%S')"
# time java -Xmx8g -jar $snpEff_dir/snpEff.jar build -c $snpEff_dir/snpEff.config.$myGenome -noCheckProtein -noCheckCds -gtf22 -v $myGenome && echo "done"

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
echo "开始注释VCF (这可能需要很长时间)... at $(date '+%Y-%m-%d %H:%M:%S')"
time java -Xmx10G -jar $snpEff_dir/snpEff.jar eff \
    -c $snpEff_dir/snpEff.config.$myGenome $myGenome ${vcf} > \
    ${out_dir_1}/${myGenome}.eff.vcf \
    -csvStats ${out_dir_1}/${myGenome}.csv \
    -stats ${out_dir_1}/${myGenome}.html && echo "snpEff done"

##### 第四步、计算同义突变和非同义突变的数量 ######
vcf=${out_dir_1}/${myGenome}.eff.vcf

# ====== 中间和输出路径 ======
tmp_dir=${out_dir}/tmp_strict
mkdir -p "$tmp_dir"

variants_bed=${tmp_dir}/variants.info.bed
intersect_tsv=${tmp_dir}/variants_in_orfs.strict.tsv
out_tsv=${out_dir}/orf_syn_nonsyn.strict.tsv

echo "[1/2] 从 VCF 提取变异位置 + REF/ALT/INFO ..."
# A: chr  start(0-based)  end(1-based)  REF  ALT  INFO
bcftools view -H "$vcf" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2,$4,$5,$8}' \
  > "$variants_bed"

echo "[2/2] 与 ORF CDS bed 取交集 ..."
# B: chr start end ORF_ID ...
# 输出：1-6 是变异信息，7- 为 ORF bed 的列（第10列是 ORF_ID）
bedtools intersect -a "$variants_bed" -b "$CDS_bed" -wa -wb \
  > "$intersect_tsv"

echo "交集文件：$intersect_tsv"
echo "接下来运行 Python 严格统计脚本："
echo "  python count_orf_syn_nonsyn_strict.py $intersect_tsv $CDS_bed $out_tsv"
python pnps/count_orf_syn_nonsyn_strict.py $intersect_tsv $CDS_bed $out_tsv