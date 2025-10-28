mkdir -p input/formatted_mgf
mkdir output

novor_script_path=/rd1/user/lit/project/sORFs/src/denovopipeline/resources/DeNovoGUI-1.16.6/
par=$novor_script_path/parameter.20250109.par
output_path=./output/novor/20250613
mkdir -p $output_path
mgf_fragpipe_output=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min_rename/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf
filename=$(basename $mgf_fragpipe_output)
grep -v "1/K0" "$mgf_fragpipe_output" > input/formatted_mgf/$filename
mgf_for_novor_input=input/formatted_mgf/$filename
java -cp $novor_script_path/DeNovoGUI-1.16.6.jar com.compomics.denovogui.cmd.DeNovoCLI \
-spectrum_files "$mgf_for_novor_input" \
-output_folder "$output_path" \
-id_params "$par" \
-directag 0 \
-pepnovo 0 \
-novor 1

mkdir -p mgf_inspection && cd mgf_inspection
mgf_fragpipe_output=/rd1/user/lit/project/sORFs/raw_data/MS/Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min_rename/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf
# 49040780
wc -l $mgf_fragpipe_output
less $mgf_fragpipe_output |grep PEPMASS=|cut -f 1 -d ' '|sed 's/PEPMASS=//' > mass.txt
less $mgf_fragpipe_output |grep CHARGE=|sed 's/CHARGE=//;s/+//' > z.txt
less $mgf_fragpipe_output |grep "1/K0="|sed 's/1\/K0=//' > K0.txt
paste mass.txt z.txt K0.txt > m.z.K0.txt

gen_m_z_k0(){
    less $1 |grep PEPMASS=|cut -f 1 -d ' '|sed 's/PEPMASS=//' > mass.txt
    less $1 |grep CHARGE=|sed 's/CHARGE=//;s/+//' > z.txt
    less $1 |grep "ION_MOBILITY="|sed 's/ION_MOBILITY=//' > K0.txt
    paste mass.txt z.txt K0.txt > m.z.K0.txt
}
mgf=/rd1/user/lit/project/sORFs/analysis/20250613_denovo/output/peaks_20250618/mgf_mzml/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.mgf
# 32975967
wc -l $mgf
mkdir peaks_output_mgf && cd peaks_output_mgf
gen_m_z_k0 $mgf