
pep_xml_1=output/MS/Fragpipe_output/2024_12_03_batch_run/Trypsin_LysC/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/interact-CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.pep.xml
pep_xml_2=output/MS/Fragpipe_output/2024_12_06_batch_run/Trypsin_LysC/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/interact.pep.xml

cd output/MS/Fragpipe_output/2025_01_07_default_uniprot_database_F2/

# 测试运行two-step (基于open search)
for enzyme in Trypsin_ArgC Trypsin_AspN Trypsin_GluC Trypsin Trypsin_LysN Trypsin_LysC Trypsin_Chymotrypsin Nonspecific Nocleavage;do
bash Uni.gene.enzyme.v1.sh $enzyme config/Open/${enzyme}_6_50_prot_1.workflow config/Open/6_50_prot_1.workflow
done

cd analysis/20250114_multiple_steps_run

proj_path=/rd1/user/lit/project/sORFs
workflow=$proj_path/config/Open/Trypsin_LysC_6_50_prot_1.workflow
database_path=$proj_path/custom_database/2024-11-21-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta.fas
sed 's|database.db-path=.*|database.db-path='"$database_path"'|' $workflow > ./fragpipe.workflow.tmp
sed 's|tab-run.write_sub_mzml=false|tab-run.write_sub_mzml=true|' ./fragpipe.workflow.tmp > ./fragpipe.workflow

find /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/ -name "*uncalibrated.mzML" > all.uncali.mzml.txt
less all.uncali.mzml.txt|grep F2 | xargs -I {} sh -c 'echo -e "{}\t$(basename {})\t\tDDA"' > fragpipe-files.fp-manifest

conda activate fragpipe

/rd1/user/lit/software/fragpipe/bin/fragpipe --headless \
--workflow ./fragpipe.workflow \
--manifest ./fragpipe-files.fp-manifest \
--workdir ./F2_step_1 \
--ram 180 \
--threads 40 \
--config-tools-folder /rd1/user/lit/project/sORFs/fragpipe_tools \
--config-diann /rd1/user/lit/project/sORFs/fragpipe_tools/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
--config-python ~/rd1/anaconda3/envs/py3/bin/python


# PTMShepherd [Work dir: /rd1/user/lit/project/sORFs/analysis/20250114_multiple_steps_run/./F2_step_1]
# java -Xmx180G -Dlibs.thermo.dir=/rd1/user/lit/project/sORFs/fragpipe_tools/tools/MSFragger-4.1/ext/thermo -cp /rd1/user/lit/software/fragpipe/tools/ptmshepherd-3.0.0.jar:/rd1/user/lit/software/fragpipe/tools/batmass-io-1.33.4.jar:/rd1/user/lit/software/fragpipe/tools/commons-math3-3.6.1.jar:/rd1/user/lit/software/fragpipe/tools/hipparchus-1.8/hipparchus-core-1.8.jar:/rd1/user/lit/software/fragpipe/tools/hipparchus-1.8/hipparchus-stat-1.8.jar:/rd1/user/lit/software/fragpipe/lib/commons-lang3-3.14.0.jar:/rd1/user/lit/software/fragpipe/lib/fragpipe-22.0.jar edu.umich.andykong.ptmshepherd.PTMShepherd "/rd1/user/lit/project/sORFs/analysis/20250114_multiple_steps_run/./F2_step_1/shepherd.config"

# PTM-Shepherd version 3.0.0
# (c) 2022 University of Michigan

# Using Java 17.0.13-internal on 184320MB memory

# Finding spectral data
#         Indexing data from CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_uncalibrated_mzML
# Fatal error: In dataset "CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_uncalibrated_mzML" could not find mzData for run CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_uncalibrated
# Process 'PTMShepherd' finished, exit code: 1
# Process returned non-zero exit code, stopping

# 测试运行two-step (基于closed search)
proj_path=/rd1/user/lit/project/sORFs
workflow=$proj_path/config/Default/Trypsin_LysC_6_50_prot_1.workflow
database_path=$proj_path/custom_database/2024-11-21-decoys-contam-uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta.fas
sed 's|database.db-path=.*|database.db-path='"$database_path"'|' $workflow > ./fragpipe.workflow.tmp
sed 's|tab-run.write_sub_mzml=false|tab-run.write_sub_mzml=true|' ./fragpipe.workflow.tmp > ./fragpipe.workflow

find /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/ -name "*uncalibrated.mzML" > all.uncali.mzml.txt
less all.uncali.mzml.txt|grep F2 | xargs -I {} sh -c 'echo -e "{}\t$(basename {})\t\tDDA"' > fragpipe-files.fp-manifest

conda activate fragpipe

/rd1/user/lit/software/fragpipe/bin/fragpipe --headless \
--workflow ./fragpipe.workflow \
--manifest ./fragpipe-files.fp-manifest \
--workdir ./F2_step_1 \
--ram 180 \
--threads 40 \
--config-tools-folder /rd1/user/lit/project/sORFs/fragpipe_tools \
--config-diann /rd1/user/lit/project/sORFs/fragpipe_tools/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
--config-python ~/rd1/anaconda3/envs/py3/bin/python

# 合并搜库结果
