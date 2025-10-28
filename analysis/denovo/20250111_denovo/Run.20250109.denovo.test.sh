#!/usr/bin/sh

################################################
#File Name: Run.20250109.denovo.test.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2025年01月09日 星期四 11时48分46秒
################################################

set -eo pipefail

# mgf 所在路径
## 5.7G
mgf_1=raw_data/MS/Guomics_SP_MSdata/converted_files/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.mgf
## 2.5G
mgf_2=raw_data/MS/Guomics_SP_MSdata/converted_files/another_way/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.mgf
# 1.7G
mgf_3=raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_uncalibrated.mgf
# 2.4G
mgf_4=raw_data/MS/Guomics_SP_MSdata/converted_files/try_1/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.mgf
# 5.1G
mgf_5=raw_data/MS/Guomics_SP_MSdata/converted_files/try_1/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.1.mgf

# 报错
python src/denovopipeline/src/main.py reformatMGF --input $mgf_4 --output output/denovo/F2_reformatted.mgf

# 修复
cp src/denovopipeline/src/ src/denovopipeline/src_1/
pushd src/denovopipeline/src_1/
find ./ -type f -name "*.py" -exec sed -i '/future_fstrings/d' {} +
find ./ -type f -name "*.py" |xargs grep -H future_fstrings
popd

# 修复后重新运行
python src/denovopipeline/src_1/main.py reformatMGF --input $mgf_4


cd /rd1/user/lit/project/sORFs/src/denovopipeline/pretrained_model
cp Casanovo/* ../resources/Casanovo/

conda activate casanovo
cd src/denovopipeline/
mgf=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/converted_files/try_1/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_reformatted.mgf
output_path=/rd1/user/lit/project/sORFs/output/denovo/casanovo/casanovo_F2

# 报错
python src_1/main.py denovo --input $mgf --output $output_path --casanovo 1 --casanovo_model resources/Casanovo/casanovo_10epochs_HumanHCDMassive.ckpt

# Also 报错
casanovo --mode=denovo --model=resources/Casanovo/casanovo_10epochs_HumanHCDMassive.ckpt --peak_path=$mgf --config=resources/Casanovo/casanovo_config.yaml --output=$output_path

# 从官网下载
cd /rd1/user/lit/project/sORFs/src/
git clone https://github.com/Noble-Lab/casanovo.git
cd casanovo/
conda create --name casanovo_env python=3.10
conda activate casanovo_env
pip install casanovo
casanovo --help
casanovo configure
# 下载预训练数据
mkdir ckpt && cd ckpt
wget https://github.com/Noble-Lab/casanovo/releases/download/v4.2.0/casanovo_v4_2_0.ckpt
wget https://github.com/Noble-Lab/casanovo/releases/download/v4.0.0/casanovo_massivekb.ckpt
wget https://github.com/Noble-Lab/casanovo/releases/download/v4.0.0/casanovo_nontryptic.ckpt
cpkt=/rd1/user/lit/project/sORFs/src/casanovo/ckpt/casanovo_v4_2_0.ckpt
config=/rd1/user/lit/project/sORFs/src/casanovo/casanovo.20250109.yaml
# de novo搜库
mgf_3=raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_uncalibrated.mgf
casanovo sequence -c $config --model $cpkt -o /rd1/user/lit/project/sORFs/output/denovo/casanovo/results.mztab $mgf_3

casanovo sequence sample_data/sample_preprocessed_spectra.mgf

###### 测试novor ######
conda activate denovopipeline
mkdir -p /rd1/user/lit/project/sORFs/output/denovo/novor/F2
output_path=/rd1/user/lit/project/sORFs/output/denovo/novor/F2
# 报错
python src_1/main.py denovo --input $mgf --output $output_path --denovogui 1

par=resources/DeNovoGUI-1.16.6/parameter.20250109.par
java -cp resources/DeNovoGUI-1.16.6/DeNovoGUI-1.16.6.jar com.compomics.denovogui.cmd.DeNovoCLI \
-spectrum_files "$mgf" \
-output_folder "$output_path" \
-id_params "$par" \
-directag 0 \
-pepnovo 0 \
-novor 1 

# 成功
java -cp resources/DeNovoGUI-1.16.6/DeNovoGUI-1.16.6.jar com.compomics.denovogui.cmd.IdentificationParametersCLI -out resources/DeNovoGUI-1.16.6/parameter.20250109.par

python -m ipykernel install --user --name=casanovo_env --display-name "casanovo_env"

# 
conda activate fragpipe
msconvert --mgf --filter "msLevel 2" --outfile raw_data/MS/Guomics_SP_MSdata/converted_files/F2_bruker.mgf raw_data/MS/Guomics_SP_MSdata/converted_files/F2_bruker.mzML

# 只有35364个记录？
mgf_5=raw_data/MS/Guomics_SP_MSdata/converted_files/F2_mzmine.mgf
mgf_6=raw_data/MS/Guomics_SP_MSdata/converted_files/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.2.mgf
mzml=raw_data/MS/Guomics_SP_MSdata/converted_files/F2_bruker.mzML
mgf_7=raw_data/MS/Guomics_SP_MSdata/converted_files/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562.3.mgf

grep -c "BEGIN IONS" $mgf_5

psm=output/MS/Fragpipe_output/2024_12_03_batch_run/Trypsin_LysC/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/psm.tsv

less $psm |cut -f 1 |cut -d "." -f 2 |tail -n +2 |sort -k 1,1n |tail -n 1

###### 最终用novor ######
conda activate base

par=resources/DeNovoGUI-1.16.6/parameter.20250109.par
output_path=/rd1/user/lit/project/sORFs/output/denovo/novor/F2_uncalibrated_mgf
mkdir -p $output_path
mgf_3=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_uncalibrated.mgf
java -cp resources/DeNovoGUI-1.16.6/DeNovoGUI-1.16.6.jar com.compomics.denovogui.cmd.DeNovoCLI \
-spectrum_files "$mgf_3" \
-output_folder "$output_path" \
-id_params "$par" \
-directag 0 \
-pepnovo 0 \
-novor 1

grep -v "1/K0" "$mgf_3" > /rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_uncalibrated.reformatted.mgf
mgf=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20241022licq_BSEP_PreExp_DDA_60min/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2-52_1_3562_uncalibrated.reformatted.mgf
java -cp resources/DeNovoGUI-1.16.6/DeNovoGUI-1.16.6.jar com.compomics.denovogui.cmd.DeNovoCLI \
-spectrum_files "$mgf" \
-output_folder "$output_path" \
-id_params "$par" \
-directag 0 \
-pepnovo 0 \
-novor 1

less $mgf_3 |grep TITLE |sed "s/TITLE=//" > output/denovo/F2_spectrum_title.txt
less $mgf_3 |grep SCANS |sed "s/SCANS=//" > output/denovo/scans.txt
