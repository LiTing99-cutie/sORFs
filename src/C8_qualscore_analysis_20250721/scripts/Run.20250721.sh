docker run --rm -it -v /rd1/user/lit/project/sORFs/src/C8_qualscore_analysis_20250721:/data proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:skyline_daily_25.1.1.174-b787f12 /bin/bash
sample=CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020
wine msconvert input/${sample}_calibrated.mzML --mzXML --zlib=false --32 -o output/${sample}.mzXML
# 首先根据脚本得到过滤之前的pep.xml文件：/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min_rename/output/db_search_20250523/Trypsin/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_d_copy_20250703/interact-CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.pep.xml
# 路径为CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_d_copy_20250703而不是默认的搜库输出
cd /rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min_rename/output/db_search_20250523/Trypsin/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_d_copy_20250703
/rd1/user/lit/software/fragpipe/tools/percolator_3_6_5/linux/percolator --only-psms --no-terminate --post-processing-tdc --num-threads 40 --results-psms CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_percolator_target_psms.tsv --decoy-results-psms CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_percolator_decoy_psms.tsv --protein-decoy-pattern rev_ CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_edited.pin
# 将默认的0.5参数改成0
# percolator.min-prob=0.5
java -cp "/rd1/user/lit/software/fragpipe/lib/*" com.dmtavt.fragpipe.tools.percolator.PercolatorOutputToPepXML CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.pin CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020 CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_percolator_target_psms.tsv CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_percolator_decoy_psms.tsv interact-CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020 DDA 0 /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min_rename/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_calibrated.mzML
# 1.*ipynb
cd ../output
cp ./CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.mzXML/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.mzXML ../
sudo rm -rf ./CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.mzXML/
jar=/rd1/user/lit/project/sORFs/src/Qualscore/qualscore_v1.0_2.jar
mv ../CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.mzXML .
java -jar $jar -o ./qualscore.res.20250721.txt -a -lp modified.pep.xml
# 2.*ipynb
# 3.*ipynb