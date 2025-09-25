nohup find . -name "*.gz" | xargs -I {} cp --parents --preserve=links {} \
	/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227 &> cp.20250408.log &
mv cp.20250408.log log/
nohup find . -name "*.sam" | xargs -I {} cp --parents --preserve=links {} \
	/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227 &> log/cp.20250408.1.log &

find . -name "*.gz" -o -name "*.sam"|xargs rm -rf

# 20250410
mv Run.combine.call.orfs.v1.20250312.sh deprecated_script/
rm -rf human_brain_ribo_merge_call_orfs_20250312