#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通用文件diff脚本：比较两个文本文件内容的差异（unified diff格式）
用法：
    python diff_workflow.py file1 file2 [--output out.txt]
"""

import difflib
import argparse
import sys

def compare_files(file1_path, file2_path):
    """
    比较两个文件内容的差异，使用unified diff格式输出。
    参数:
        file1_path (str): 第一个文件的路径
        file2_path (str): 第二个文件的路径
    返回:
        list: 包含差异的列表，每项为一行diff结果
    """
    try:
        with open(file1_path, 'r', encoding='utf-8') as file1, open(file2_path, 'r', encoding='utf-8') as file2:
            file1_lines = file1.readlines()
            file2_lines = file2.readlines()
        diff = list(difflib.unified_diff(
            file1_lines, file2_lines,
            fromfile=file1_path, tofile=file2_path,
            lineterm=''
        ))
        return diff
    except FileNotFoundError as e:
        print(f"文件未找到: {e}", file=sys.stderr)
        return []
    except Exception as e:
        print(f"发生错误: {e}", file=sys.stderr)
        return []

def main():
    parser = argparse.ArgumentParser(description="比较两个文本文件内容的差异（unified diff格式）")
    parser.add_argument('file1', help='第一个文件路径')
    parser.add_argument('file2', help='第二个文件路径')
    parser.add_argument('--output', '-o', help='将diff结果输出到指定文件（可选）')
    args = parser.parse_args()

    diff = compare_files(args.file1, args.file2)
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as fout:
            for line in diff:
                fout.write(line + '\n')
        print(f"Diff结果已写入: {args.output}")
    else:
        for line in diff:
            print(line)

if __name__ == '__main__':
    main()