#!/usr/bin/env python3
"""
Transcript Coverage Analysis Script v1_20250725
分析随着样本数量增加，reads覆盖的转录本数量变化

Author: Assistant
Date: 2025-07-25
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import random
import os
import argparse
from pathlib import Path

def load_count_data(file_path):
    """
    加载featureCounts输出文件
    """
    print(f"Loading data from: {file_path}")
    
    # 读取数据，跳过注释行
    data = pd.read_csv(file_path, sep='\t', comment='#')
    
    # 获取样本列（从第7列开始）
    sample_cols = data.columns[6:]
    
    print(f"Found {len(sample_cols)} samples")
    print(f"Found {len(data)} transcripts")
    
    return data, sample_cols

def calculate_coverage_curves(data, sample_cols, n_permutations=100, min_reads=1):
    """
    计算不同样本数量下的转录本覆盖情况
    
    Parameters:
    - data: 包含count数据的DataFrame
    - sample_cols: 样本列名列表
    - n_permutations: 随机排列的次数
    - min_reads: 最小reads数阈值
    """
    
    # 确保sample_cols是列表格式
    sample_cols = list(sample_cols)
    
    results = []
    
    # 对每个样本数量进行多次随机排列
    for n_samples in range(1, len(sample_cols) + 1):
        print(f"Processing {n_samples} samples...")
        
        n_detected_list = []
        
        for _ in range(n_permutations):
            # 随机选择n_samples个样本
            selected_samples = random.sample(sample_cols, n_samples)
            
            # 计算这些样本中每个转录本的总reads数
            transcript_totals = data[selected_samples].sum(axis=1)
            
            # 统计有reads的转录本数量
            n_detected = (transcript_totals >= min_reads).sum()
            n_detected_list.append(n_detected)
        
        # 计算统计量
        mean_detected = np.mean(n_detected_list)
        std_detected = np.std(n_detected_list)
        min_detected = np.min(n_detected_list)
        max_detected = np.max(n_detected_list)
        
        results.append({
            'n_samples': n_samples,
            'mean_detected': mean_detected,
            'std_detected': std_detected,
            'min_detected': min_detected,
            'max_detected': max_detected,
            'n_permutations': n_permutations
        })
    
    return pd.DataFrame(results)

def plot_coverage_curves(results, output_dir):
    """
    绘制覆盖曲线图
    """
    plt.figure(figsize=(12, 8))
    
    # 主曲线
    plt.plot(results['n_samples'], results['mean_detected'], 
             'o-', linewidth=2, markersize=6, label='Mean detected transcripts')
    
    # 误差条
    plt.fill_between(results['n_samples'], 
                     results['mean_detected'] - results['std_detected'],
                     results['mean_detected'] + results['std_detected'],
                     alpha=0.3, label='±1 SD')
    
    # 最小最大值范围
    plt.fill_between(results['n_samples'], 
                     results['min_detected'],
                     results['max_detected'],
                     alpha=0.1, label='Min-Max range')
    
    plt.xlabel('Number of Samples', fontsize=12)
    plt.ylabel('Number of Detected Transcripts', fontsize=12)
    plt.title('Transcript Coverage vs Number of Samples', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 保存图片
    output_file = os.path.join(output_dir, 'transcript_coverage_curves.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Coverage curve saved to: {output_file}")
    
    plt.show()

def calculate_saturation_analysis(results):
    """
    计算饱和度分析
    """
    # 计算边际增益（每增加一个样本检测到的转录本数量）
    results['marginal_gain'] = results['mean_detected'].diff()
    results['marginal_gain_rate'] = results['marginal_gain'] / results['mean_detected'].shift(1) * 100
    
    # 找到边际增益率低于5%的点（饱和度点）
    saturation_point = None
    for i, rate in enumerate(results['marginal_gain_rate']):
        if not pd.isna(rate) and rate < 5:
            saturation_point = results.iloc[i]
            break
    
    return results, saturation_point

def save_results(results, output_dir):
    """
    保存结果到文件
    """
    # 保存详细结果
    output_file = os.path.join(output_dir, 'transcript_coverage_results.csv')
    results.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    # 保存摘要
    summary_file = os.path.join(output_dir, 'transcript_coverage_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Transcript Coverage Analysis Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total samples analyzed: {len(results)}\n")
        f.write(f"Total permutations per sample size: {results['n_permutations'].iloc[0]}\n\n")
        
        f.write("Coverage Statistics:\n")
        f.write("-" * 20 + "\n")
        for _, row in results.iterrows():
            f.write(f"{row['n_samples']} samples: {row['mean_detected']:.0f} ± {row['std_detected']:.0f} transcripts\n")
        
        # 饱和度分析
        results_with_saturation, saturation_point = calculate_saturation_analysis(results)
        if saturation_point is not None:
            f.write(f"\nSaturation Analysis:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Saturation point: {saturation_point['n_samples']} samples\n")
            f.write(f"Transcripts at saturation: {saturation_point['mean_detected']:.0f}\n")
            f.write(f"Marginal gain rate at saturation: {saturation_point['marginal_gain_rate']:.2f}%\n")
    
    print(f"Summary saved to: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Analyze transcript coverage vs number of samples')
    parser.add_argument('--input', '-i', required=True, 
                       help='Input featureCounts file path')
    parser.add_argument('--output', '-o', default='./transcript_coverage_analysis',
                       help='Output directory')
    parser.add_argument('--permutations', '-p', type=int, default=100,
                       help='Number of permutations per sample size (default: 100)')
    parser.add_argument('--min-reads', '-m', type=int, default=1,
                       help='Minimum reads threshold for transcript detection (default: 1)')
    
    args = parser.parse_args()
    
    # 创建输出目录
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 加载数据
    data, sample_cols = load_count_data(args.input)
    
    # 计算覆盖曲线
    print(f"\nCalculating coverage curves with {args.permutations} permutations...")
    results = calculate_coverage_curves(data, sample_cols, 
                                      n_permutations=args.permutations,
                                      min_reads=args.min_reads)
    
    # 绘制图形
    print("\nGenerating plots...")
    plot_coverage_curves(results, output_dir)
    
    # 保存结果
    print("\nSaving results...")
    save_results(results, output_dir)
    
    print(f"\nAnalysis completed! Results saved to: {output_dir}")

if __name__ == "__main__":
    main() 