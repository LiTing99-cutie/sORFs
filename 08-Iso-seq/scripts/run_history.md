### Iso-seq分析：前处理
nohup bash isoseq_bulk_preprocess_v1_20250710.sh &> ../log/isoseq_bulk_preprocess_v1_20250710.20250710.log &
- Conda环境：biotools
- 输入文件：
  - ../rawdata/r84130_250703_001_1_A01/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam
  - 参考基因组：/home/user/data/lit/database/public/genome/hg38/hg38.fa
- 线程数：30
- 日志文件：../log/isoseq_bulk_preprocess_v1_20250710.20250710.log
- 1和2成功运行，3失败了
### Iso-seq分析：前处理
nohup bash isoseq_bulk_preprocess_v2_20250711.sh &> ../log/isoseq_bulk_preprocess_v2_20250711.20250711.log &
- 备注：系统自带的/home/user/data3/rbase/opt/pacbio/smrtlink/smrtcmds/bin/pbmm2无法识别bam文件，使用conda install -c bioconda pbmm2下载了最新的pbmm2 1.17.0，注释前两行，重新运行3
### Iso-seq分析：转录本注释
nohup bash isoseq_bulk_postprocess_v1_20250710.sh &> ../log/isoseq_bulk_postprocess_v1_20250710.20250711.log &

### 20250716 饱和度分析
python3 saturation_curve_data_generator_v1_20250716.py /home/user/data3/lit/project/sORFs/08-Iso-seq/processed/classify/collapsed_classification.filtered_lite_classification.txt --max-samples 100 ../results/_ &> ../log/saturation_curve_data_generator_v1_20250716.20250716.log

### 20250716 运行sqanti
mkdir -p ../results/sqanti
cd ../results/sqanti
# 手动修改
sqanti3 init -c sqanti3_config.yaml -a cpus=50 dir=.
cd ../../scripts
nohup bash sqanti.20250716.sh &>../log/sqanti.20250716.log &

### 20250716 合并不同日期下机的数据，并且重新生成pbi索引
nohup bash merge.bam.20250716.sh &> ../log/merge.bam.20250716.log &

### 20250717 
##### sqanti中的filter以及rescue步骤运行失败，重新运行
nohup bash sqanti.20250717.sh &>../log/sqanti.20250717.log &
- 备注：修改了config.yaml中的路径
##### sqanti中的rescue步骤运行失败，重新运行
nohup bash sqanti.20250717.1.sh &>../log/sqanti.20250717.1.log &
- 备注：注释了config.yaml中的一行参数

### 20250717 合并不同日期下机的数据，并且重新生成pbi索引
nohup bash merge.bam.20250717.sh &> ../log/merge.bam.20250717.log &
- 备注：之前合并bam时没有合并header，导致后面结果有些异常
bash rename.header.20250717.sh
- 备注：上面的合并方式不对，无法生成index，直接使用bam_1的index替换合并后bam的index
bash rename.header.20250718.sh
- 备注：上面的合并方式出错了，不小心还是使用了原来的bam_1，使用bam_1的header替换合并后bam的header会有一些问题，无法生成index,因此就把合并后的bam header中的sample替换一下就可以了