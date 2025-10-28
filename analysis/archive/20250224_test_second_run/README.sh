# 解释脚本目的

### R脚本导出了output/trypsin_all_method_sep_unique.txt
database=/rd1/user/lit/project/sORFs/custom_database/Ribo_ORFs_add_assemble_20250125/uniprot.nonCano.sorf.20250125.fa
database_1=/rd1/user/lit/project/sORFs/custom_database/Transcriptome_based_20250209/uniprot.nonCano.sorf.20250209.trans_based.fa
uniprot_fa=/rd1/user/lit/project/sORFs/custom_database/uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta
seqkit grep -f output/trypsin_all_method_sep_unique.txt $database $database_1 |seqkit rmdup -n > output/trypsin_all_method_sep_unique.fa
cat $uniprot_fa output/trypsin_all_method_sep_unique.fa > output/uniprot.trypsin_all_method_sep_unique.fa

# 运行搜库
mkdir -p log output/search
nohup bash Run.search.sh &> log/search.log &

# 整理结果，需要比较前后差别的小肽
mkdir -p output/stat output/visual
path=/rd1/user/lit/project/sORFs/analysis/20250224_test_second_run
pushd /rd1/user/lit/project/sORFs/analysis/20250210_stat_mouse_formal_data
for run in default open;do
Rscript S1a1.Uni.Organize_all_samples.R \
    $path/output/search/$run \
    $path/output/stat/$run \
    $PWD/output/sample_metadata.rds
Rscript S1b1.Uni.Visualize_all_samples.R \
    $path/output/stat/$run \
    $path/output/visual/$run \
    $PWD/output/sample_order.rds \
    $PWD/output/sample_metadata.rds
done
popd

# 使用*pdf合并会有点异常，使用tiff和png合并更好
output_path=/rd1/user/lit/project/sORFs/analysis/20250224_test_second_run/output
cd ${output_path}/sampled.psm.filtered/merge
montage ../*pdf -tile 5x10 -geometry +0+0 output.pdf
montage ../*png -tile 10x5 -geometry +0+0 output.png

cd ${output_path}/sampled.psm.kept/merge
montage ../*.tiff -tile 10x5 -geometry +0+0 output.tiff

cd ${output_path}/sampled.psm.kept.over_2
montage *png -tile 10x5 -geometry +0+0 1-output.png

output_path=/rd1/user/lit/project/sORFs/analysis/20250224_test_second_run/output
for output in sampled.psm.filtered sampled.psm.kept sampled.psm.kept.over_2 sampled.psm.kept.over_3;do
cd ${output_path}/${output}
[ -f 1-output.png ] && rm -rf 1-output.png
montage *png -tile 10x5 -geometry +0+0 1-output.png
done

for output in first_run second_run;do
cd ${output_path}/${output}
[ -f 1-output.png ] && rm -rf 1-output.png
montage *png -tile 10x5 -geometry +0+0 1-output.png
done