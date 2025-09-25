# 相关的之前的位置的文件可以删除
file1=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250527/rawdata/p21_40_0425.R1.raw.fastq.gz
file2=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250509/raw_data/p21_40_0425.raw.fastq.gz
## 不同的长度，应该也可以一起分析
cat $file1 $file2 > p21_40_0425.fq.gz
file1=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250509/raw_data/p21_40_0422.raw.fastq.gz
file2=/home/user/data3/licq/BSEP/riboseq/rawdata/20250528/p21_40_0422.raw.fastq.gz
cat $file1 $file2 > p21_40_0422.fq.gz
for file in p21_0523_5 p21_0523_6 p21_0523_7;do
file1=/home/user/data3/licq/BSEP/riboseq/rawdata/20250605/$file.R1.raw.fastq.gz
file2=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250625/$file.R1.raw.fastq.gz
cat $file1 $file2 > $file.fq.gz
done
# 直接拷贝源文件，而不是软链接（如果要拷贝软链接，就需要加上-d参数）
ls /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/organize_all_test_data_20250515/*_{1,2,3,4}.fq.gz|xargs -I {} cp {} ./ 

# 20250723 整合数据
nohup bash Run.20250723.sh &