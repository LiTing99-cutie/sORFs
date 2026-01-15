import pandas as pd
import os
import matplotlib.pyplot as plt
import glob

# ================= 配置区域 =================
# 你的结果主目录
BASE_DIR = "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/processed/saturation_analysis"
# 输出统计表的文件名
OUTPUT_CSV = "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/processed/saturation_analysis/saturation_stats.csv"
# 输出图片的文件名
OUTPUT_PLOT = "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/processed/saturation_analysis/saturation_curve.png"

def main():
    print(f"开始读取目录: {BASE_DIR}")
    
    summary_list = []
    
    # 循环读取 10% 到 100%
    for depth in range(10, 101, 10):
        # 构造文件路径
        file_path = os.path.join(BASE_DIR, f"depth_{depth}", "RiboCode_result_collapsed.txt")
        
        if os.path.exists(file_path):
            try:
                # 读取 TSV 文件
                df = pd.read_csv(file_path, sep='\t')
                
                # 1. 统计总数
                total_orfs = len(df)
                
                # 2. 按 ORF_type 分类统计 (annotated, novel, uORF, etc.)
                type_counts = df['ORF_type'].value_counts().to_dict()
                
                # 3. 汇总数据
                row_data = {
                    'Depth_Percent': depth,
                    'Total_ORFs': total_orfs
                }
                # 将分类统计合并进去
                row_data.update(type_counts)
                
                summary_list.append(row_data)
                print(f"Depth {depth}%: Found {total_orfs} ORFs")
                
            except Exception as e:
                print(f"Warning: 读取 depth_{depth} 失败: {e}")
        else:
            print(f"Warning: 文件不存在 {file_path}")

    # 转换为 DataFrame
    if not summary_list:
        print("未找到任何结果文件，请检查路径。")
        return

    result_df = pd.DataFrame(summary_list)
    
    # 整理列顺序：Depth 和 Total 放前面，其他类型自动排
    cols = ['Depth_Percent', 'Total_ORFs'] + [c for c in result_df.columns if c not in ['Depth_Percent', 'Total_ORFs']]
    result_df = result_df[cols].fillna(0) # 空值填0
    
    # 排序
    result_df = result_df.sort_values('Depth_Percent')

    # ================= 保存表格 =================
    result_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n统计表格已保存为: {OUTPUT_CSV}")
    print("-" * 30)
    print(result_df.to_string(index=False))
    print("-" * 30)

    # ================= 绘制曲线图 =================
    plt.figure(figsize=(10, 6))
    
    # 画总数曲线
    plt.plot(result_df['Depth_Percent'], result_df['Total_ORFs'], 
             marker='o', linewidth=2, label='Total ORFs', color='black')
    
    # 画主要类型的曲线 (Annotated vs Novel)
    # 这里的列名取决于 RiboCode 实际输出的类型，通常有 annotated, novel, uORF 等
    interesting_types = ['annotated', 'novel'] 
    colors = ['blue', 'red']
    
    for i, orf_type in enumerate(interesting_types):
        if orf_type in result_df.columns:
            plt.plot(result_df['Depth_Percent'], result_df[orf_type], 
                     marker='s', linestyle='--', label=f'{orf_type}', color=colors[i])

    plt.title('Saturation Analysis of Ribo-seq ORF Identification', fontsize=14)
    plt.xlabel('Sequencing Depth (%)', fontsize=12)
    plt.ylabel('Number of Identified ORFs', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.xticks(range(10, 110, 10))
    
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"饱和度曲线图已保存为: {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()