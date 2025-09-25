conda activate base
mgf_processor.py -i /home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/20250708_denovo_c8/input/ -o /home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/20250708_denovo_c8/output/ --remove-prefixes RTINSECONDS= SCANS= 1/K0=
mgf_processor.py -i /home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/20250708_denovo_c8/input/ -o /home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/20250708_denovo_c8/output/withScan --remove-prefixes RTINSECONDS= 1/K0=
mgf=/home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/20250708_denovo_c8/output/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf
mkdir -p output/denovo
### primenovo ###
mkdir -p output/denovo/primenovo
source activate PrimeNovo-1
ckpt=/home/user/data3/lit/software/pi-PrimeNovo/model_massive.ckpt
# export PYTHONPATH="/home/user/data3/lit/software/pi-PrimeNovo:$PYTHONPATH"
# ~ 1h 
python -m PrimeNovo.PrimeNovo --mode=denovo --peak_path=$mgf --model=$ckpt --output=output/denovo/primenovo/denovo.tsv
mv denovo.tsv output/denovo/primenovo/
### casanovo ###
mkdir -p output/denovo/casanovo
source activate casanovo
config=/home/user/data3/lit/software/casanovo/config/config.yaml
cpkt=/home/user/data3/lit/software/casanovo/ckpt/casanovo_v4_2_0.ckpt
# ~ 1h 
casanovo sequence -c $config --model $cpkt -o $PWD/output/denovo/casanovo/results.mztab $mgf
mgf_1=/home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/20250708_denovo_c8/output/withScan/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf
casanovo sequence -c $config --model $cpkt -o $PWD/output/denovo/casanovo/results_1.mztab $mgf_1