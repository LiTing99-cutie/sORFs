source activate biotools
# python merge_orf_gpe.py &> ../log/merge_orf_gpe.log
python merge_orf_gpe.v1.py &> ../log/merge_orf_gpe.v1.log
nohup bash flnc_mapping.sh &> ../log/flnc_mapping.log &