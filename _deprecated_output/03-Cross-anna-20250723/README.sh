# 20250211
project_path=/home/user/data3/lit/project/sORFs/01-ribo-seq
add_assemble_gpe=$project_path/output/assembled_trans/stringtie_output/gencode.vM29.add_assemble.genePred.txt
# run Generate_sorf_genepred.ipynb
cd gen_genepred
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
genePredToGtf file essen_info_to_gen_genepred.sorfs.gpe essen_info_to_gen_genepred.sorfs.gtf
translate_script=$project_path/S3.0c.Uni.translate_gtf.sh
bash $translate_script essen_info_to_gen_genepred.sorfs.gtf $fa
seqkit fx2tab cds.fa > cds.tab
seqkit fx2tab prot.fa > prot.tab

cd ..
Rscript $project_path/analysis/mouse_run_20250125/S3.1d.Run.annotate_ORF_type.v1.1.R \
	$PWD/output/merge_2_reformat_1.sorf.txt \
	$PWD/output/merge_2_reformat_1.sorf.add_type.txt \
	$add_assemble_gpe

# 对所有的Ribo-seq鉴定得到的sorf都生成对应的gpe（为了得到起始密码子）
jupyter nbconvert --to script Generate_sorf_genepred.ipynb
rm -rf Generate_sorf_genepred.py
conda activate base
python Uni.generate_sorf_genepred.py \
	./gen_genepred/essen_info.sorfs.all.txt \
	$add_assemble_gpe \
	./gen_genepred/sorfs.all.gpe

cd gen_genepred
prefix=sorfs.all
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
genePredToGtf file $prefix.gpe $prefix.gtf
translate_script=$project_path/S3.0c.Uni.translate_gtf.sh
bash $translate_script $prefix.gtf $fa
seqkit fx2tab cds.fa > $prefix.cds.tab
seqkit fx2tab prot.fa > $prefix.prot.tab

##### 20250303 修改之前的计算方法 #####
conda activate base
bash Uni.tmp.sh $PWD/output/MS_minus_sampled.txt $PWD/output/te_epxr/MS_minus_sampled yes
bash Uni.tmp.sh $PWD/output/MS_sep/merged_final.txt $PWD/output/te_epxr/MS_plus yes
gtf=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
bash Uni.i.gtf.o.te_expr.sh $gtf $PWD/output/te_epxr/Anno_pro yes

##### 20250402 重命名 #####
mv Uni.gen.genepred.f_sorf_id.py Uni.gen.genepred.i_sorf_id.py
mv Uni.generate_sorf_genepred.py Uni.gen.genepred.i_sorf_info.py