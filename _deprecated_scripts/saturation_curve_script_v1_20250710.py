#!/usr/bin/env python3
"""
Isoform Saturation Curve Analysis Script
Version: v1_20250710
Date: 2025-07-10

This script generates saturation curves for isoform detection by randomly sampling
different numbers of reads and counting the number of isoforms detected.

Usage:
    python saturation_curve_script_v1_20250710.py input_file output_prefix
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import random
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def load_data(file_path):
    """
    Load the isoform data from file
    """
    try:
        # Try to detect separator automatically
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            df = pd.read_csv(file_path, sep='\t')
        else:
            # Try common separators
            for sep in ['\t', ',', ';']:
                try:
                    df = pd.read_csv(file_path, sep=sep)
                    if len(df.columns) > 1:
                        break
                except:
                    continue
            else:
                raise ValueError("Could not determine file format")
        
        print(f"Successfully loaded data with {len(df)} rows and {len(df.columns)} columns")
        return df
    except Exception as e:
        print(f"Error loading file: {e}")
        return None

def prepare_reads_data(df, isoform_col, reads_col):
    """
    Prepare reads data for sampling
    """
    # Create a list where each read is represented by its isoform
    reads_list = []
    for _, row in df.iterrows():
        isoform = row[isoform_col]
        reads_count = row[reads_col]
        
        # Handle NaN or missing values
        if pd.isna(reads_count) or reads_count <= 0:
            continue
            
        # Add the isoform to the list reads_count times
        reads_list.extend([isoform] * int(reads_count))
    
    print(f"Total reads: {len(reads_list)}")
    print(f"Unique isoforms: {len(set(reads_list))}")
    
    return reads_list

def calculate_saturation_curve(reads_list, sample_sizes, n_replicates=10):
    """
    Calculate saturation curve by sampling different numbers of reads
    """
    results = []
    
    for sample_size in sample_sizes:
        sample_isoforms = []
        
        for replicate in range(n_replicates):
            # Randomly sample reads
            if sample_size >= len(reads_list):
                sampled_reads = reads_list
            else:
                sampled_reads = random.sample(reads_list, sample_size)
            
            # Count unique isoforms
            unique_isoforms = len(set(sampled_reads))
            sample_isoforms.append(unique_isoforms)
        
        # Calculate statistics
        mean_isoforms = np.mean(sample_isoforms)
        std_isoforms = np.std(sample_isoforms)
        
        results.append({
            'sample_size': sample_size,
            'mean_isoforms': mean_isoforms,
            'std_isoforms': std_isoforms,
            'min_isoforms': np.min(sample_isoforms),
            'max_isoforms': np.max(sample_isoforms)
        })
    
    return pd.DataFrame(results)

def plot_saturation_curve(results_df, output_prefix):
    """
    Plot the saturation curve
    """
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: Saturation curve
    ax1.plot(results_df['sample_size'], results_df['mean_isoforms'], 
             'o-', linewidth=2, markersize=6, label='Mean isoforms detected')
    
    # Add error bars
    ax1.fill_between(results_df['sample_size'], 
                     results_df['mean_isoforms'] - results_df['std_isoforms'],
                     results_df['mean_isoforms'] + results_df['std_isoforms'],
                     alpha=0.3, label='±1 SD')
    
    ax1.set_xlabel('Number of Reads Sampled', fontsize=12)
    ax1.set_ylabel('Number of Isoforms Detected', fontsize=12)
    ax1.set_title('Isoform Detection Saturation Curve', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Efficiency curve (isoforms per read)
    efficiency = results_df['mean_isoforms'] / results_df['sample_size']
    ax2.plot(results_df['sample_size'], efficiency, 's-', 
             linewidth=2, markersize=6, color='orange', label='Detection efficiency')
    
    ax2.set_xlabel('Number of Reads Sampled', fontsize=12)
    ax2.set_ylabel('Isoforms per Read', fontsize=12)
    ax2.set_title('Detection Efficiency vs Sample Size', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'{output_prefix}_saturation_curve.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_saturation_curve.pdf', bbox_inches='tight')
    print(f"Plots saved as {output_prefix}_saturation_curve.png and .pdf")
    
    plt.show()

def save_results(results_df, output_prefix):
    """
    Save results to file
    """
    output_file = f'{output_prefix}_saturation_results.txt'
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")
    
    # Print summary
    print("\n=== SATURATION ANALYSIS SUMMARY ===")
    print(f"Total reads analyzed: {results_df['sample_size'].max()}")
    print(f"Maximum isoforms detected: {results_df['mean_isoforms'].max():.1f}")
    print(f"Detection efficiency at max: {results_df['mean_isoforms'].max() / results_df['sample_size'].max():.4f} isoforms/read")
    
    # Check if saturation is reached
    last_quarter = results_df.tail(len(results_df)//4)
    if len(last_quarter) > 1:
        slope = np.polyfit(last_quarter['sample_size'], last_quarter['mean_isoforms'], 1)[0]
        print(f"Slope in last quarter: {slope:.6f}")
        if slope < 0.01:
            print("✓ Saturation appears to be reached")
        else:
            print("⚠ Saturation may not be reached yet")

def main():
    parser = argparse.ArgumentParser(description='Generate isoform saturation curves')
    parser.add_argument('input_file', help='Input file with isoform and reads data')
    parser.add_argument('output_prefix', help='Output prefix for results and plots')
    parser.add_argument('--isoform-col', default='isoform', help='Column name for isoform IDs (default: isoform)')
    parser.add_argument('--reads-col', default='FL.BioSample_6', help='Column name for read counts (default: FL.BioSample_6)')
    parser.add_argument('--replicates', type=int, default=10, help='Number of replicates for each sample size (default: 10)')
    parser.add_argument('--max-samples', type=int, default=20, help='Maximum number of sample sizes to test (default: 20)')
    
    args = parser.parse_args()
    
    print("=== Isoform Saturation Curve Analysis ===")
    print(f"Input file: {args.input_file}")
    print(f"Output prefix: {args.output_prefix}")
    print(f"Isoform column: {args.isoform_col}")
    print(f"Reads column: {args.reads_col}")
    
    # Load data
    df = load_data(args.input_file)
    if df is None:
        return
    
    # Check if required columns exist
    if args.isoform_col not in df.columns:
        print(f"Error: Column '{args.isoform_col}' not found in data")
        print(f"Available columns: {list(df.columns)}")
        return
    
    if args.reads_col not in df.columns:
        print(f"Error: Column '{args.reads_col}' not found in data")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Prepare reads data
    reads_list = prepare_reads_data(df, args.isoform_col, args.reads_col)
    if not reads_list:
        print("No valid reads data found")
        return
    
    # Define sample sizes (logarithmically spaced)
    total_reads = len(reads_list)
    sample_sizes = np.logspace(1, np.log10(total_reads), args.max_samples, dtype=int)
    sample_sizes = np.unique(sample_sizes)  # Remove duplicates
    
    print(f"Testing {len(sample_sizes)} sample sizes: {sample_sizes}")
    
    # Calculate saturation curve
    print("Calculating saturation curve...")
    results_df = calculate_saturation_curve(reads_list, sample_sizes, args.replicates)
    
    # Plot and save results
    plot_saturation_curve(results_df, args.output_prefix)
    save_results(results_df, args.output_prefix)
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main() 