##### 零、尝试运行 #####

nohup bash Ribo.Run.Uni.test.sh &> log/Ribo.Run.Uni.test.log &

# modify adapter; start from trim galore
nohup bash Ribo.Run.Uni.test.sh &> log/Ribo.Run.Uni.test.log &

# modify bowtie index path
nohup bash Ribo.Run.Uni.test.sh &> log/Ribo.Run.Uni.test.log &

# source activate ribocode
nohup bash Ribo.Run.Uni.test.sh &> log/Ribo.Run.Uni.test.log &

# 2022_NSMB_MatthewK_development_neocortex这篇文章的reads有点奇怪
adapter=" TGGAATTCTCGGGTGCCAAGG -a  GTTCAGAGTTCTACAGTCCGACGATC"
raw_fastq=/home/user/data3/licq/peptidomics/PublicData/mouse_brain/2022_NSMB_MatthewK_development_neocortex/ribo/SRR14048711.fastq.gz
trim_galore --adapter $adapter -j 8 -q 20 --length 20 $raw_fastq --gzip -o output/trimmed_fastq --basename modi_adapter &> log/trim_galore.modi_adapter.log
trim_galore --adapter "file:./multiple_adapters.fa" -j 8 -q 20 --length 20 $raw_fastq --gzip -o output/trimmed_fastq --basename modi_adapter_1 &> log/trim_galore.modi_adapter_1.log
fastqc -o output/fastqc -t 10 output/trimmed_fastq/modi_adapter_1_trimmed.fq.gz &> log/trimmed_fastqc.log

##### 一、批量运行 #####
nohup bash mouse_brain_run.20241004.sh &> log/mouse_brain_run.20241004.log &

# 检查批量运行的结果
## 发现tRNA等用的是人中的序列，纠正

# 重新organize一下ncRNA的文件夹
ln -s $PWD/Pre-Run/ncRNA $PWD/annot/ncRNA/human

# 建立好鼠的污染RNA的index之后重新运行
mkdir mouse_brain_output_20241011
nohup bash mouse_brain_run.20241011.sh &> log/mouse_brain_run.20241011.log &

# 检查是否为链特异性文库
find ./ -name "*_Aligned.sortedByCoord.out.bam" |tail -n +3 | xargs -i /home/user/BGM/lit/anaconda3/envs/py2/bin/infer_experiment.py -i {} \
-q 255 -r /home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.bed > mouse_brain_output_20241011/whetherStranded.txt

##### 二、运行其他两个软件的结果 #####
# screen 测试 导出到2015_Science_ChoJ_hippocampus/SRR2163083/output/Ribo-ORFs-clean
bash Multiple-tools-test-clean.sh 
# screen 运行 
## 把之前测试样本的结果copy一份到old，删除原来的结果，运行后再把RiboCode的结果加回
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR
2163083/output/Ribo-ORFs /home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/
output/Ribo-ORFs-old
find $PWD/mouse_brain_output_20241011/ -name output | parallel -j 10 --joblog log/PRICE.RibORF.prl.log 'log_path={//}/log;bash PRICE.RibORF.uni.sh {} &> $log_path/PRICE.RibORF.log'
## 运行完毕，加回RiboCode的结果
cp ./mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/output/Ribo-ORFs-old/RiboCode ./mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/output/Ribo-ORFs

##### 三、整合 #####
find ./mouse_brain_output_20241011/ -name "Ribo-ORFs" | parallel -j 10 --joblog log/Merge_res.prl.log 'cd {};bash /home/user/data3/lit/project/sORFs/01-ribo-seq/Merge_res.sh'
mkdir -p ./mouse_brain_output_20241011/Ribo_ORFs_merge
find ./mouse_brain_output_20241011/ -name "nonCano.sorf.meta.merge.txt"|xargs cat > ./mouse_brain_output_20241011/Ribo_ORFs_merge/nonCano.sorf.meta.merge.sample.txt

find ./mouse_brain_output_20241011/ -name "Ribo-ORFs" | parallel -j 10 --joblog log/Rmdup.prl.log 'cd {};bash /home/user/data3/lit/project/sORFs/01-ribo-seq/Uni.rmdup.sh'
find ./mouse_brain_output_20241011/ -name "nonCano.sorf.meta.merge.txt"|xargs cat > ./mouse_brain_output_20241011/Ribo_ORFs_merge/nonCano.sorf.meta.merge.sample.txt

find ./mouse_brain_output_20241011/ -name "Ribo-ORFs" | parallel -j 10 'cd {};bash /home/user/data3/lit/project/sORFs/01-ribo-seq/Uni.merge_res.raw.sh'
find ./mouse_brain_output_20241011/ -name "nonCano.sorf.meta.merge.raw.3_ways.txt"|xargs cat > ./mouse_brain_output_20241011/Ribo_ORFs_merge/nonCano.sorf.meta.merge.raw.3_ways.all_samples.txt

# 在R脚本中进行去重和meta table的生成
seqkit tab2fx mouse_brain_output_20241011/Ribo_ORFs_merge/nonCano.sorf.tab > mouse_brain_output_20241011/Ribo_ORFs_merge/nonCano.sorf.fa

##### 四、过滤 #####
bash Test.filter.sh

##### 2024-12-06 放宽松参数 重新运行上面的二到四 #####
# call ORFs
find $PWD/mouse_brain_output_20241011/ -name output | \
parallel -j 10 --joblog log/Uni.loose.parameter.20241210.log 'log_path={//}/log;bash Uni.loose.parameter.sh {} Ribo_ORFs_loose_para_20241206 &> $log_path/Uni.loose.parameter.20241210.log'
# organize results
find $PWD/mouse_brain_output_20241011/ -name Ribo_ORFs_loose_para_20241206 | \
parallel -j 10 --joblog log/Organize_res.20241212.log 'bash S3.0.Uni.Organize_res_v1.sh {}'
# merge,filter and add meta columns
bash S3.1.Uni.Merge_Filter_Annotate.v1.sh Ribo_ORFs_loose_para_20241206 ./mouse_brain_output_20241011/Ribo_ORFs_loose_para_20241206
# visualization
output_path=./mouse_brain_output_20241011/Ribo_ORFs_loose_para_20241206
mkdir -p  $output_path/visual/
Rscript S3.2.Visual.1.1.R $output_path/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt \
	$output_path/filter/nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.filtered.txt  $output_path/visual/

# 生成fa文件去搜库
mkdir -p database/custom/Ribo_ORFs_loose_para_20241206 && pushd database/custom/Ribo_ORFs_loose_para_20241206
output_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/Ribo_ORFs_loose_para_20241206
uniprot_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform_2024_10_25.fasta
less $output_path/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt | tail -n +2 | awk '{print $2,$4}' |seqkit tab2fx > nonCano.sorf.20241213.loose_para.fa
cat $uniprot_fasta nonCano.sorf.20241213.loose_para.fa > uniprot.nonCano.sorf.20241213.loose_para.fa
popd

##### 2024-12-18 对卡松参数前的结果也重新运行上面的三到四 #####
output_path=./mouse_brain_output_20241011/Ribo-ORFs-20241218
# organize results
find $PWD/mouse_brain_output_20241011/ -name Ribo-ORFs | \
parallel -j 10 --joblog log/Organize_res.20241218.log 'bash S3.0.Uni.Organize_res_v2.sh {} 0.05 1 fixed'
# merge,filter and add meta columns
bash S3.1.Uni.Merge_Filter_Annotate.v1.sh Ribo-ORFs $output_path
# visualization
mkdir -p  $output_path/visual/
Rscript S3.2.Visual.1.1.R $output_path/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt \
	$output_path/filter/nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.filtered.txt  $output_path/visual/

##### 基于转录组的库 #####
# 生成fa文件去搜库
cp annot/sORF_database/Transcriptome_based/ database/custom/
pushd database/custom/Transcriptome_based
uniprot_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform_2024_10_25.fasta
less trans_based_specific_sorfs.txt | tail -n +2 | awk '{print $1,$3}' |seqkit tab2fx > nonCano.sorf.20241218.trans_based_specific.fa
cat $uniprot_fasta nonCano.sorf.20241218.trans_based_specific.fa > uniprot.nonCano.sorf.20241218.trans_based_specific.fa

less trans_based_sorfs.txt | tail -n +2 | awk '{print $1,$3}' |seqkit tab2fx > nonCano.sorf.20241218.trans_based.fa
cat $uniprot_fasta nonCano.sorf.20241218.trans_based.fa > uniprot.nonCano.sorf.20241218.trans_based.fa
popd
