i_sorf_id_o_gtf(){
    # 输入sorf id的文件（需要包含一列列名为ORF_id_trans），在目录下生成gpe和gtf文件
    sorfs_id=$1
    gen_genepred_i_sorf_id_script=/home/user/data3/lit/project/sORFs/03-Cross-anna/Uni.gen.genepred.i_sorf_id.py
    conda activate base
    source /home/user/data2/lit/bin/lit_utils.sh
    define_annotation_gencode_v41_human
    python $gen_genepred_i_sorf_id_script \
        $sorfs_id \
        $gpe_15 \
        target.gpe
    genePredToGtf file target.gpe target.gtf
}

i_gtf_id_o_cds_bed6(){
    # 输入gtf，在目录下生成gtf文件
    gtf=$1
    less $gtf |awk '$3=="CDS"'|awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > target.bed6
}

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

sorfs_id_file=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/tmp/sorf_id.20250606.txt
p_site_bam=p_sites.bam
i_sorf_id_o_gtf $sorfs_id_file
i_gtf_id_o_cds_bed6 target.gtf
time bedtools coverage -a target.bed6 -b $p_site_bam -d -s > target.bedtools.cov.bedgraph
add_pos_frame target.bedtools.cov.bedgraph > target.bedtools.cov.add.pos.frame.bedgraph
# 检查
grep ENST00000417816.2-chr10:20785746-20812796 target.bedtools.cov.add.pos.frame.bedgraph

calcu_frame_1_cnt target.bedtools.cov.add.pos.frame.bedgraph > target.frame1.txt