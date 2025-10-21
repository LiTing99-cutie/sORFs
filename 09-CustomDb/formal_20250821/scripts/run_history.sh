# 0.为不同的ORF鉴定软件构建注释
nohup bash ribo.seq.tools.build.anno.20250821.sh &> ../log/ribo.seq.tools.build.anno.20250821.log &

# 1.创建自定义的数据库
## 1.1拷贝chunfu的脚本，进行自定义修改
cp /home/user/data3/rbase/translation_pred/models/src/process/run.sh build.custom.db.20250822.sh
cp /home/user/data3/rbase/translation_pred/models/src/process/ORF_combiner_generator.py ORF_combiner_generator.py
## 1.2使用iso-seq的gtf以及整合了纯合SNP的基因组
nohup bash build.custom.db.20250822.sh &> ../log/build.custom.db.20250822.log &

## 1.3整合污染
contam_fasta=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/contaminant_fasta/2022_JPR_contam.fasta
mkdir -p ../results/custom_db_20250826/
cat /home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.fa $contam_fasta > ../results/custom_db_20250826/human_brain_custom_db.fasta

# 2.计算protein的长度等信息
conda activate base
python fasta_to_tsv.py ../processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.fa -o ../processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.info.txt

# 3.创建整合ribo-seq信息的数据库【弃用】
nohup bash build.custom.db.use.ribo.20250827.sh &> ../log/build.custom.db.use.ribo.20250827.log &
## 3.1.整合污染
contam_fasta=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/contaminant_fasta/2022_JPR_contam.fasta
out_dir=../results/custom_db_20250904/
mkdir -p $out_dir
cat ../processed/annotation/RibORF_annot/candidate_ORFs_ribo_filtering/candidateORF.6aa.long.M.rmdup.pep.fa $contam_fasta > $out_dir/human_brain_custom_db_ribo_fil.fasta
## 重命名ID（公司需要，但是比较麻烦，最后直接让公司按照顺序生成编号去搜库）
<!-- awk '/^>/{gsub(/\+/, "plus"); gsub(/-/, "minus"); gsub(/[=|:() ]/,"."); gsub(/[.]+/,"."); print; next}1' \
../results/custom_db_20250904/human_brain_custom_db_ribo_fil.fasta > ../results/custom_db_20250904/human_brain_custom_db_ribo_fil.RenameID.fasta -->

# 4. 从冗余ID中选择代表性ID

python pick_representative.20250910.py \
  -i ../processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.dup.txt \
  -e /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.txt \
  -o ../processed/annotation/RibORF_annot/candidate_ORFs/representative.tsv \
  --seed 42

## 4.1 得到picked id以及representative id的对应
awk -v OFS='\t' '{print $3,$2}' ../processed/annotation/RibORF_annot/candidate_ORFs/representative.tsv > ../processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt

# 5. 基于出核信息和翻译信息进行过滤
nohup bash build.custom.db.filtering.20250917.sh &> ../log/build.custom.db.filtering.20250917.log &

# 6. 整理数据库中所有的ORF
id_convert_script=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/scripts/id.convert.20251017.sh
proj_path=/home/user/data3/lit/project/sORFs
file_to_be_convert=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.info.txt
id_map=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt
file_converted=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.info.rePicked.txt
bash $id_convert_script -c 1 --header yes $file_to_be_convert "$id_map" $file_converted

proj_path=/home/user/data3/lit/project/sORFs
file_converted=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.info.rePicked.txt
augment_orf_table_script=$proj_path/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/scripts/augment_orf_table.v2.20251017.sh
make_orf_id_custom_script=$proj_path/09-CustomDb/formal_20250821/scripts/make_orf_id_custom.py
output_path=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id_convert
mkdir -p $output_path
bash $augment_orf_table_script -n ID $file_converted $output_path
conda activate base
python $make_orf_id_custom_script --in $output_path/augmented.tsv \
  --out $output_path/augmented.add.custom.id.tsv \
  --old-col ID \
  --codon5 ORF_start \
	--codon3 ORF_end \
  --aa-col Seq