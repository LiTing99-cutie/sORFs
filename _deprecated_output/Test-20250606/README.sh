mkdir log
nohup bash /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/Run.20250606.sh &> log/Run.log &


#### 另外一个模块的测试 ####
1.offset.curation.20250606.sh

# 手动check之后选择最终的offset
awk -F'\t' '$4 != ""' all_offset_tab.txt|cut -f 1,4 > all_offset_tab.curated.txt

2.*sh

nohup bash Run.20250716.sh &> log/Run.20250716.log &
