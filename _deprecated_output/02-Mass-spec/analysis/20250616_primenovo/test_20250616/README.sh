/home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/primenovo_20250616/test_20250616/top1k.ipynb
grep -v '^$' mgf_test_c8/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.top.1k.mgf > mgf_test_c8/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_uncalibrated.top.1k.rmBlank.mgf

grep -v "1/K0" mgf_test_c8/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_uncalibrated.top.1k.rmBlank.mgf > mgf_test_c8/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_uncalibrated.top.1k.rmBlank.rmK0.mgf
cd /home/user/data3/lit/software/pi-PrimeNovo
python -m PrimeNovo.PrimeNovo --mode=eval --peak_path=./bacillus.10k.mgf --model=./model_massive.ckpt
mgf=/home/user/data3/lit/project/sORFs/02-Mass-spec/analysis/primenovo_20250616/test_20250616/mgf_test_c8/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf
python -m PrimeNovo.PrimeNovo --mode=denovo --peak_path=$mgf --model=./model_massive.ckpt

mkdir test
cp $mgf ./test
conda activate base
# 转换下mgf的格式
python /home/user/data2/lit/software/pUniFind/mgf_processor.py -i $PWD/test -o ./c8.processed.mgf -p 8
# 保存上次的结果
mkdir res_20250617
mv denovo.tsv res_20250617/
# 使用转换后的格式跑primenovo
mgf=/home/user/data3/lit/software/pi-PrimeNovo/c8.processed.mgf/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020_uncalibrated.mgf
conda activate PrimeNovo-1
nohup python -m PrimeNovo.PrimeNovo --mode=denovo --peak_path=$mgf --model=./model_massive.ckpt &

mkdir res_20250624
mv denovo.tsv res_20250624/
# 修改容差为20ppm
nohup python -m PrimeNovo.PrimeNovo --mode=denovo --peak_path=$mgf --model=./model_massive.ckpt &
