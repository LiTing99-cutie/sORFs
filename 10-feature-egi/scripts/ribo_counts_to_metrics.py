import pandas as pd
import numpy as np
import sys

def calculate_ribo_metrics(input_file, output_file):
    """
    处理Ribo-seq counts数据，计算RPKM、TPM和密码子覆盖度
    
    参数:
    input_file: 输入的counts文件路径
    output_file: 输出文件路径
    """
    
    # 读取数据，跳过注释行
    df_list = []
    with open(input_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                df_list.append(line.strip())
    
    # 创建DataFrame（假设第一行是header）
    import io
    data_str = '\n'.join(df_list)
    df = pd.read_csv(io.StringIO(data_str), sep='\t', header=0)
    
    # 假设列结构: col1=ORF_id, col6=length_bp, col7=RPF_reads, col8=Psites_number
    # 根据awk脚本的列索引调整（awk中$1是第1列，$6是第6列等）
    df.columns = [f'col{i+1}' for i in range(len(df.columns))]
    
    # 提取必要的列
    orf_id = df['col1']
    length_bp = pd.to_numeric(df['col6'], errors='coerce').fillna(0)
    rpf_reads = pd.to_numeric(df['col7'], errors='coerce').fillna(0)
    psites_number = pd.to_numeric(df['col8'], errors='coerce').fillna(0)
    
    # 计算总reads
    rpf_total = rpf_reads.sum()
    ps_total = psites_number.sum()
    
    # 创建结果DataFrame
    results = pd.DataFrame({
        'ORF_id': orf_id,
        'RPF_reads': rpf_reads.astype(int),
        'Psites_number': psites_number.astype(int),
    })
    
    # 计算RPKM
    # RPKM = (reads * 10^9) / (length_bp * total_reads)
    results['RPF_RPKM'] = np.where(
        (length_bp > 0) & (rpf_total > 0),
        (rpf_reads * 1e9) / (length_bp * rpf_total),
        np.nan
    )
    
    results['Psites_RPKM'] = np.where(
        (length_bp > 0) & (ps_total > 0),
        (psites_number * 1e9) / (length_bp * ps_total),
        np.nan
    )
    
    # 计算TPM (Transcripts Per Million)
    # Step 1: 计算RPK (reads per kilobase)
    rpf_rpk = np.where(length_bp > 0, rpf_reads / (length_bp / 1000), 0)
    ps_rpk = np.where(length_bp > 0, psites_number / (length_bp / 1000), 0)
    
    # Step 2: 计算TPM
    rpf_rpk_sum = rpf_rpk.sum()
    ps_rpk_sum = ps_rpk.sum()
    
    results['RPF_TPM'] = np.where(
        rpf_rpk_sum > 0,
        (rpf_rpk / rpf_rpk_sum) * 1e6,
        np.nan
    )
    
    results['Psites_TPM'] = np.where(
        ps_rpk_sum > 0,
        (ps_rpk / ps_rpk_sum) * 1e6,
        np.nan
    )
    
    # 计算每密码子覆盖度 (codon coverage)
    # coverage = (reads * 3) / length_bp
    results['RPF_codon_coverage'] = np.where(
        length_bp > 0,
        (rpf_reads * 3) / length_bp,
        np.nan
    )
    
    results['Psites_codon_coverage'] = np.where(
        length_bp > 0,
        (psites_number * 3) / length_bp,
        np.nan
    )
    
    # 格式化输出，将NaN替换为"NA"
    for col in ['RPF_RPKM', 'Psites_RPKM', 'RPF_TPM', 'Psites_TPM', 
                'RPF_codon_coverage', 'Psites_codon_coverage']:
        results[col] = results[col].apply(lambda x: 'NA' if pd.isna(x) else f'{x:.6f}')
    
    # 保存结果
    results.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    
    # 打印统计信息
    print(f"处理完成！")
    print(f"总RPF reads: {rpf_total:,}")
    print(f"总P-sites: {ps_total:,}")
    print(f"处理了 {len(results)} 个ORF")
    print(f"结果已保存到: {output_file}")
    
    return results

# 使用示例
if __name__ == "__main__":
    # 定义文件路径
    input_file = "../processed/ribo/counts.txt"
    output_file = "../processed/feature_preprare/orf.rpf.psite.1.txt"
    
    # 运行分析
    results_df = calculate_ribo_metrics(input_file, output_file)
    
    # 可选：显示前几行结果
    print("\n结果预览（前5行）：")
    print(results_df.head())