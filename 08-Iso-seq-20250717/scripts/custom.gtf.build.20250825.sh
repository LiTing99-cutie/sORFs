ref_gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
pb_gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf
classify_txt=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/processed/classify/collapsed_classification.filtered_lite_classification.txt
out_dir=../processed/custom_gtf_20250825
grep \"PB.2340.1\" $pb_gtf > $out_dir/test.pb.gtf
cat <(head -n1 $classify_txt) <(grep -e "PB.2340.1\b" $classify_txt) > $out_dir/test.txt

# 初始版本
python3 custom.db.20250825.py \
  --ref_gtf /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf \
  --pb_gtf  $out_dir/test.pb.gtf \
  --classify_txt $out_dir/test.txt

# 输出的转录本feature按照特定顺序排序
# 结果中gene id是associated gene id而不是原来的pb id
# 小数据测试
python3 custom.db.20250825.v2.py \
  --ref_gtf /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf \
  --pb_gtf  $out_dir/test.pb.gtf \
  --classify_txt $out_dir/test.txt

# 正式数据
python3 custom.db.20250825.v2.py \
  --ref_gtf /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf \
  --pb_gtf  $pb_gtf \
  --classify_txt $classify_txt

# 负链的输出中，外显子的坐标需要从大到小
# 输出的gtf中有gene feature行
python3 custom.db.20250825.v3.py \
  --ref_gtf /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf \
  --pb_gtf  $pb_gtf \
  --classify_txt $classify_txt

python3 custom.db.20250825.v4.py \
  --ref_gtf /home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf \
  --pb_gtf  $pb_gtf \
  --classify_txt $classify_txt
gtfToGenePred -geneNameAsName2 -genePredExt ../results/custom.gtf.with_orf.gtf ../results/custom.gtf.with_orf.15.gpe
