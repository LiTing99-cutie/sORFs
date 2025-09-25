###### 24-36长度中三碱基周期性大于等于50%且显著性小于0.01的reads数量 ######
# 对于同一个样本，计算不同层次的Ribo-seq Reads结果
conda activate ribocode
# 得到24-36长度中三碱基周期性大于等于50%且显著性小于0.01的reads长度对应的bam
sample=human_brain_ribo_1
pushd $PWD/human_brain_output_20250227/2020-Nature/human_brain_ribo_1/output/Ribo_ORFs/RiboCode
bam_1=../../alignment/${sample}_Aligned.toTranscriptome.out.bam
bam_2=../../alignment/${sample}_Aligned.sortedByCoord.out.bam
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RiboCode_annot
metaplots -a $RiboCode_annot -r $bam_1 -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o config_0.5_0.01
grep "^#" config_0.5_0.01_pre_config.txt|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.tab.txt
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250509/extract_specific_length_reads.sh
bash $script $bam_2 offset.tab.txt ../../alignment/${sample}_Aligned.sortedByCoord.out.withPeri.bam
popd

# find ./ -name "*withPeri.bam" | xargs -I {} sh -c 'cnt=$(samtools view -c "{}"); echo -e "{}\t$cnt"' > $output_dir/p_sites_number.txt
find ./ -name "*withPeri.bam" | xargs -I {} sh -c 'cnt=$(samtools view -c "{}"); echo -e "{}\t$cnt"'

###### 5. 总体的三碱基周期性 ######
output_dir=stat
define_annotation_gencode_v41_human
conda activate ribocode
for bam in $(find $PWD -name "human_brain_ribo_1*Aligned.sortedByCoord.out.bam");do
[ -f $bam.bai ] || samtools index $bam
ribotish quality -b $bam -g $gtf -p 30
done
# 编码基因所有codon的不同frame的比例
script=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/ribosome_profiling/parse_ribotish_qual.py
mkdir -p $output_dir/ribotish/
for sample in human_brain_ribo_1;do
echo "$sample"
python $script --sample_name $sample \
     --txt_path $PWD/human_brain_output_20250227/2020-Nature/${sample}/output/alignment/${sample}_Aligned.sortedByCoord.out_qual.txt \
	 --offset_path $PWD/human_brain_output_20250227/2020-Nature/${sample}/output/alignment/${sample}_Aligned.sortedByCoord.out.bam.para.py \
	 --RPF_start_distr_file $output_dir/ribotish/$sample.RPF_start_distr_file.txt \
	 --RPF_stop_distr_file $output_dir/ribotish/$sample.RPF_stop_distr_file.txt \
	 --frame_distr_file $output_dir/ribotish/$sample.frame_distr_file.txt
done
# find ./ -name "*frame_distr_file.txt" |xargs cat > $output_dir/ribotish/Frame_distr_file.txt
find ./ -name "*frame_distr_file.txt" |xargs cat