data3_dir=/home/user/data3/lit/project/sORFs
backup_dir=/home/user/data/lit/project/sORFs
### Iso-seq原始数据备份与清理
mkdir -p $backup_dir/08-Iso-seq/rawdata
cp -r $data3_dir/08-Iso-seq/rawdata/r84130_250703_001_1_A01 \
    $data3_dir/08-Iso-seq/rawdata/sequencing_add_20250716 \
    $backup_dir/08-Iso-seq/rawdata
rm -rf $data3_dir/08-Iso-seq/rawdata/
# 只保留results中的sqanti结果
rm -rf $data3_dir/08-Iso-seq/processed/
rm -rf $data3_dir/08-Iso-seq/ref_set/
nohup cp -r $data3_dir/08-Iso-seq-20250717/rawdata/merge $backup_dir/08-Iso-seq/rawdata &
rm -rf $data3_dir/08-Iso-seq-20250717/rawdata

### 基因组数据原始数据备份与清理
mkdir -p $backup_dir/07-Genome/rawdata 
nohup cp -r $data3_dir/07-Genome/rawdata  $backup_dir/07-Genome/ &
rm -rf $data3_dir/07-Genome/rawdata

### 核内核外数据原始数据备份与清理
mkdir -p $backup_dir/06-RNA-seq/rawdata/
nohup cp --dereference $data3_dir/06-RNA-seq/01-rawdata/organized_20250624/* $backup_dir/06-RNA-seq/rawdata/ &

### Ribo-seq数据原始数据备份与清理
mkdir -p $backup_dir/01-ribo-seq/rawdata
cd $backup_dir/01-ribo-seq/rawdata 
mkdir -p ref_other_lab in_house/test/human/ in_house/test/mouse/ in_house/formal public/Disome public/human_brain public/mouse_brain 
nohup cp -r $data3_dir/01-ribo-seq/rawdata/Test-20241221/* $backup_dir/01-ribo-seq/rawdata/ref_other_lab &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/public-disome-20250614/* $backup_dir/01-ribo-seq/rawdata/public/Disome &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/Mouse_E16_test/cleandata/* $backup_dir/01-ribo-seq/rawdata/in_house/test/mouse/ &
nohup cp -r --dereference $data3_dir/01-ribo-seq/rawdata/organize_all_test_data_20250515/E16* $backup_dir/01-ribo-seq/rawdata/in_house/test/mouse/ &


cat $data3_dir/01-ribo-seq/analysis/Test-20250306/md5.txt \
    $data3_dir/01-ribo-seq/analysis/Test-20250408/md5.txt \
    $data3_dir/01-ribo-seq/rawdata/in_house/test/mouse/clean_md5.txt |grep E16 > $backup_dir/01-ribo-seq/rawdata/in_house/test/mouse/raw_md5.txt

nohup cp -r $data3_dir/01-ribo-seq/rawdata/in_house_human_organized_20250625/p21_40* $backup_dir/01-ribo-seq/rawdata/in_house/test/human/ &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/in_house_human_organized_20250625/p21_0523_{1..4}* $backup_dir/01-ribo-seq/rawdata/in_house/test/human/ &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/in_house_human_organized_20250625/p21_0523_{5..7}* $backup_dir/01-ribo-seq/rawdata/in_house/formal/demo &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/data-20250509/raw_data/* $backup_dir/01-ribo-seq/rawdata/in_house/test/human/demo/ &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/data-20250527/rawdata/* $backup_dir/01-ribo-seq/rawdata/in_house/test/human/add/ &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/data-20250625/* $backup_dir/01-ribo-seq/rawdata/in_house/formal/add/ &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/data-20250714/rawdata/* $backup_dir/01-ribo-seq/rawdata/in_house/formal/demo/ &
nohup cp -r $data3_dir/data-download-20250721/MJ20250701254-MJ-D-20250629001-李春琼-纯文库-6个样本/rawdata/* $backup_dir/01-ribo-seq/rawdata/in_house/formal/add/ &
nohup cp -r $data3_dir/01-ribo-seq/rawdata/in_house_human_organized_20250625/p21_0626_*gz $backup_dir/01-ribo-seq/rawdata/in_house/formal/merge/ &

cat /home/user/data3/lit/project/sORFs/data-download-20250722/MJ20250701254-MJ-D-20250629001-李春琼-纯文库-9个样本/rawdata/raw_md5.txt \
    /home/user/data3/lit/project/sORFs/data-download-20250721/MJ20250701254-MJ-D-20250629001-李春琼-纯文库-6个样本/rawdata/raw_md5.txt  \
    /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250625/raw_md5.txt |sort -k2,2 > $backup_dir/01-ribo-seq/rawdata/in_house/formal/add/raw_md5.txt

cat /home/user/data3/lit/project/sORFs/06-RNA-seq/01-rawdata/20250613/md5.txt \
    /home/user/data3/lit/project/sORFs/06-RNA-seq/01-rawdata/MJ20250407316-MJ-R-20250418027/cleandata/clean_md5.txt > \
    /home/user/data/lit/project/sORFs/06-RNA-seq/rawdata/raw_md5.txt
