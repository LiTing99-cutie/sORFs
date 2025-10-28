project_path=/rd1/user/lit/project/sORFs/analysis/20250707_c8_organize
cd $project_path
mgf_processor.py -i /rd1/user/lit/project/sORFs/src/SPEQ_analysis_20250707/input/C8 -o $PWD/input --remove-prefixes RTINSECONDS= 1/K0=
# cp /rd1/user/lit/project/sORFs/src/SPEQ_analysis_20250707/output/C8_processed_20250707_cleaner/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf ./input
mgf=$PWD/input/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf
mkdir -p output/denovo
### novor ###
source activate base
novor_script_path=/rd1/user/lit/project/sORFs/src/denovopipeline/resources/DeNovoGUI-1.16.6
par=$novor_script_path/parameter.20250707.par
output_path=output/denovo/novor/
mkdir -p $output_path
java -cp $novor_script_path/DeNovoGUI-1.16.6.jar com.compomics.denovogui.cmd.DeNovoCLI \
-spectrum_files "$mgf" \
-output_folder "$output_path" \
-id_params "$par" \
-directag 0 \
-pepnovo 0 \
-novor 1
### primenovo ###
# From NanJing
### casanovo ###
# From NanJing