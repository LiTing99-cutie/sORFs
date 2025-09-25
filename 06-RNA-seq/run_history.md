# 20250718 合并技术重复
nohup bash merge.bam.20250718.sh &> log/merge.bam.20250718.log &
# 计算reads数量
nohup bash merge.bam.20250719.sh &> log/merge.bam.20250719.log &

# 第一批数据 
Run.20250508.sh

# 第二批数据【完整数据】
Run.20250621.sh

# 使用新的gtf，计算表达量，使用二代的结果应该和三代结果差不多，也就是三代的这些基因应该都表达量大于0
nohup bash Run.20250821.sh &> ./log/Run.20250821.log &

# 使用iso-seq的gtf
nohup bash Run.20250909.isoseq.gtf.sh &> ./log/Run.20250909.isoseq.gtf.log &