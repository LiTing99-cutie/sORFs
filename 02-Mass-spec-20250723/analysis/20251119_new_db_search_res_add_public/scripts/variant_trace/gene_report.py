#!/usr/bin/env python3
"""
生成突变覆盖分析的可视化报告
"""

import pandas as pd
import sys
from typing import Dict, List


def generate_html_report(result_file: str, output_html: str):
    """
    生成HTML格式的可视化报告
    """
    df = pd.read_csv(result_file, sep='\t')
    
    # 计算统计信息
    total_mutations = len(df)
    covered_mutations = df['is_covered'].sum()
    has_peptides = df['has_peptides'].sum()
    orf_found = df['orf_found'].sum()
    
    coverage_rate = (covered_mutations / total_mutations * 100) if total_mutations > 0 else 0
    
    # 按突变类型统计
    mutation_type_stats = []
    if 'mutation_type' in df.columns:
        for mut_type in df['mutation_type'].unique():
            subset = df[df['mutation_type'] == mut_type]
            covered = subset['is_covered'].sum()
            total = len(subset)
            pct = (covered / total * 100) if total > 0 else 0
            mutation_type_stats.append({
                'type': mut_type,
                'covered': covered,
                'total': total,
                'percentage': pct
            })
    
    # 生成HTML
    html_content = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>突变-肽段覆盖分析报告</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-card.success {{
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
        }}
        .stat-card.warning {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        }}
        .stat-card.info {{
            background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        }}
        .stat-value {{
            font-size: 36px;
            font-weight: bold;
            margin: 10px 0;
        }}
        .stat-label {{
            font-size: 14px;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 14px;
        }}
        th {{
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ddd;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .covered {{
            color: #27ae60;
            font-weight: bold;
        }}
        .not-covered {{
            color: #e74c3c;
            font-weight: bold;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 12px;
            font-weight: bold;
        }}
        .badge-success {{
            background-color: #d4edda;
            color: #155724;
        }}
        .badge-danger {{
            background-color: #f8d7da;
            color: #721c24;
        }}
        .badge-warning {{
            background-color: #fff3cd;
            color: #856404;
        }}
        .progress-bar {{
            width: 100%;
            height: 30px;
            background-color: #ecf0f1;
            border-radius: 15px;
            overflow: hidden;
            margin: 20px 0;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #11998e 0%, #38ef7d 100%);
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            transition: width 0.3s ease;
        }}
        .peptide-seq {{
            font-family: 'Courier New', monospace;
            background-color: #f8f9fa;
            padding: 2px 4px;
            border-radius: 3px;
            font-size: 12px;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            text-align: center;
            color: #7f8c8d;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 突变-肽段覆盖分析报告</h1>
        
        <h2>📊 总体统计</h2>
        <div class="stats-grid">
            <div class="stat-card info">
                <div class="stat-label">总突变数</div>
                <div class="stat-value">{total_mutations}</div>
            </div>
            <div class="stat-card success">
                <div class="stat-label">被覆盖的突变</div>
                <div class="stat-value">{covered_mutations}</div>
            </div>
            <div class="stat-card warning">
                <div class="stat-label">未覆盖的突变</div>
                <div class="stat-value">{total_mutations - covered_mutations}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">覆盖率</div>
                <div class="stat-value">{coverage_rate:.1f}%</div>
            </div>
        </div>
        
        <div class="progress-bar">
            <div class="progress-fill" style="width: {coverage_rate}%">
                {coverage_rate:.1f}%
            </div>
        </div>
        
        <h2>📈 按突变类型统计</h2>
        <table>
            <thead>
                <tr>
                    <th>突变类型</th>
                    <th>被覆盖</th>
                    <th>总数</th>
                    <th>覆盖率</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for stat in mutation_type_stats:
        html_content += f"""
                <tr>
                    <td><span class="badge badge-warning">{stat['type']}</span></td>
                    <td>{stat['covered']}</td>
                    <td>{stat['total']}</td>
                    <td><strong>{stat['percentage']:.1f}%</strong></td>
                </tr>
"""
    
    html_content += """
            </tbody>
        </table>
        
        <h2>🔍 详细结果</h2>
        <table>
            <thead>
                <tr>
                    <th>ORF ID</th>
                    <th>突变位置</th>
                    <th>氨基酸变化</th>
                    <th>覆盖状态</th>
                    <th>覆盖肽段</th>
                    <th>肽段位置</th>
                    <th>样本</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for _, row in df.iterrows():
        orf_id_short = row['orf_id'].split(':')[0] if pd.notna(row['orf_id']) else 'N/A'
        
        if row['is_covered']:
            status_html = '<span class="badge badge-success">✓ 已覆盖</span>'
            status_class = 'covered'
        else:
            status_html = '<span class="badge badge-danger">✗ 未覆盖</span>'
            status_class = 'not-covered'
        
        peptides = row['covering_peptides'] if pd.notna(row['covering_peptides']) and row['covering_peptides'] else '-'
        positions = row['peptide_positions'] if pd.notna(row['peptide_positions']) and row['peptide_positions'] else '-'
        samples = row['covering_samples'] if pd.notna(row['covering_samples']) and row['covering_samples'] else '-'
        
        # 格式化肽段序列
        if peptides != '-':
            peptide_list = peptides.split('|')
            peptides_html = '<br>'.join([f'<span class="peptide-seq">{p}</span>' for p in peptide_list])
        else:
            peptides_html = '-'
        
        html_content += f"""
                <tr>
                    <td title="{row['orf_id']}">{orf_id_short}</td>
                    <td><strong>{row['mutation_position']}</strong></td>
                    <td><span class="badge badge-warning">{row['aa_change']}</span></td>
                    <td class="{status_class}">{status_html}</td>
                    <td>{peptides_html}</td>
                    <td>{positions}</td>
                    <td>{samples}</td>
                </tr>
"""
    
    html_content += """
            </tbody>
        </table>
        
        <div class="footer">
            <p>🔬 生成时间: """ + pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S') + """</p>
            <p>突变-肽段覆盖分析工具 v1.0</p>
        </div>
    </div>
</body>
</html>
"""
    
    # 保存HTML文件
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"HTML报告已生成: {output_html}")


def generate_summary_report(result_file: str, output_txt: str):
    """
    生成文本格式的摘要报告
    """
    df = pd.read_csv(result_file, sep='\t')
    
    with open(output_txt, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("突变-肽段覆盖分析摘要报告\n")
        f.write("="*80 + "\n\n")
        
        # 总体统计
        f.write("【总体统计】\n")
        f.write(f"  总突变数: {len(df)}\n")
        f.write(f"  找到ORF序列: {df['orf_found'].sum()}\n")
        f.write(f"  有肽段鉴定: {df['has_peptides'].sum()}\n")
        f.write(f"  被肽段覆盖: {df['is_covered'].sum()}\n")
        
        if len(df) > 0:
            coverage_rate = (df['is_covered'].sum() / len(df)) * 100
            f.write(f"  总覆盖率: {coverage_rate:.2f}%\n")
        f.write("\n")
        
        # 按突变类型统计
        if 'mutation_type' in df.columns:
            f.write("【按突变类型统计】\n")
            for mut_type in df['mutation_type'].unique():
                subset = df[df['mutation_type'] == mut_type]
                covered = subset['is_covered'].sum()
                total = len(subset)
                pct = (covered / total * 100) if total > 0 else 0
                f.write(f"  {mut_type}: {covered}/{total} ({pct:.1f}%)\n")
            f.write("\n")
        
        # 已覆盖的突变
        covered_df = df[df['is_covered'] == True]
        if len(covered_df) > 0:
            f.write("【已覆盖的突变】\n")
            for _, row in covered_df.iterrows():
                f.write(f"  • {row['orf_id'].split(':')[0]}\n")
                f.write(f"    位置: {row['mutation_position']} ({row['aa_change']})\n")
                if pd.notna(row['covering_peptides']):
                    f.write(f"    肽段: {row['covering_peptides']}\n")
                if pd.notna(row['peptide_positions']):
                    f.write(f"    位置: {row['peptide_positions']}\n")
                f.write("\n")
        
        # 未覆盖的突变
        not_covered_df = df[df['is_covered'] == False]
        if len(not_covered_df) > 0:
            f.write("【未覆盖的突变】\n")
            for _, row in not_covered_df.iterrows():
                f.write(f"  • {row['orf_id'].split(':')[0]}\n")
                f.write(f"    位置: {row['mutation_position']} ({row['aa_change']})\n")
                if not row['has_peptides']:
                    f.write(f"    原因: 该ORF没有肽段鉴定\n")
                else:
                    f.write(f"    原因: 肽段未覆盖此位置\n")
                    if pd.notna(row['total_peptides_for_orf']):
                        f.write(f"    该ORF有 {int(row['total_peptides_for_orf'])} 个肽段\n")
                f.write("\n")
        
        f.write("="*80 + "\n")
    
    print(f"文本报告已生成: {output_txt}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_report.py <result_tsv> [output_html] [output_txt]")
        print("\n示例:")
        print("  python generate_report.py mutation_coverage_results.tsv")
        print("  python generate_report.py mutation_coverage_results.tsv report.html summary.txt")
        sys.exit(1)
    
    result_file = sys.argv[1]
    output_html = sys.argv[2] if len(sys.argv) > 2 else result_file.replace('.tsv', '_report.html')
    output_txt = sys.argv[3] if len(sys.argv) > 3 else result_file.replace('.tsv', '_summary.txt')
    
    try:
        generate_html_report(result_file, output_html)
        generate_summary_report(result_file, output_txt)
        print("\n✓ 报告生成完成!")
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)