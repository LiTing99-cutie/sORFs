bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250909_org_all_data/processed/all.offsetCorrected.merged.sorted.bam
gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/augment_orf_table/sub.gtf
mkdir tmp && cd tmp
less $gtf |awk '$3=="CDS"'|\
awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > ms_sorfs.bed6
sorf_id="PB.10024.19:chr12:-|28|3718:492:2442|canonical|ATG"
grep -F $sorf_id ms_sorfs.bed6 > $sorf_id.bed6
# 大概需要3-4min
bedtools coverage -a $sorf_id.bed6 -b $bam -d -s > $sorf_id.bedgraph

add_pos_frame(){
    bedgraph=$1
    awk '
    BEGIN {
        OFS = "\t";  # 输出字段分隔符设为制表符
    }
    {
        key = $4;    # 按V4列分组
        strand[key] = $6;  # 保存链方向（+/-）
        count[key]++;       # 统计每组的行数
        # 保存整行数据，索引为 (key, 行号)
        rows[key, count[key]] = $0;
    }
    END {
        # 遍历每个分组
        for (key in count) {
            total = count[key];  # 当前组的行数
            is_positive = (strand[key] == "+");  # 判断链方向
            # 遍历组内每行
            for (i = 1; i <= total; i++) {
                # 生成Position：正序或逆序
                pos = is_positive ? i : (total - i + 1);
                # 生成Frame：正向链1-3循环，反向链3-1循环
                frame = is_positive ? ((i - 1) % 3 + 1) : (3 - (i - 1) % 3);
                # 输出原始行 + Position + Frame
                print rows[key, i], pos, frame;
            }
        }
    }' $bedgraph
}

calcu_frame_1_cnt(){
    awk '
    BEGIN {
        OFS = "\t";
        print "Group", "Total_Count", "Frame1_Count";
    }
    {
        key = $4;          # 按第4列分组
        sum_all[key] += $8; # 计算每组第8列的总和
        if ($10 == 1) {
            sum_frame1[key] += $8; # 当第10列=3时，累加第8列的值
        }
    }
    END {
        # 计算并输出每组结果
        for (key in sum_all) {
            printf "%s\t%d\t%d\n", key, sum_all[key], sum_frame1[key];
        }
    }' $1
}
add_pos_frame $sorf_id.bedgraph > target.bedtools.cov.add.pos.frame.bedgraph
calcu_frame_1_cnt target.bedtools.cov.add.pos.frame.bedgraph > target.frame1.txt