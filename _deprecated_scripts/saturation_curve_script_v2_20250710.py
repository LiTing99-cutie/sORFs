#!/usr/bin/env python3
"""
Isoform Saturation Curve Analysis Script (Improved Version)
Version: v2_20250710
Date: 2025-07-10

This script generates saturation curves for isoform detection by randomly sampling
different numbers of reads and counting the number of isoforms detected.
IMPROVED: Better visualization of standard deviation bands.

Usage:
    python saturation_curve_script_v2_20250710.py input_file output_prefix
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
            'max_isoforms': np.max(sample_isoforms),
            'cv_isoforms': std_isoforms / mean_isoforms if mean_isoforms > 0 else 0  # Coefficient of variation
        })
    
    return pd.DataFrame(results)

def plot_saturation_curve(results_df, output_prefix):
    """
    Plot the saturation curve with improved error band visualization
    """
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
    
    # Plot 1: Saturation curve with improved error bands
    # Calculate error band bounds
    y_lower = results_df['mean_isoforms'] - results_df['std_isoforms']
    y_upper = results_df['mean_isoforms'] + results_df['std_isoforms']
    
    # Ensure lower bound doesn't go below 0
    y_lower = np.maximum(y_lower, 0)
    
    # Add error bands first (so they appear behind the line)
    ax1.fill_between(results_df['sample_size'], y_lower, y_upper,
                     alpha=0.4, color='lightcoral', label='±1 SD', zorder=1)
    
    # Plot the main line on top
    ax1.plot(results_df['sample_size'], results_df['mean_isoforms'], 
             'o-', linewidth=2, markersize=6, color='red', 
             label='Mean isoforms detected', zorder=2)
    
    # Add individual error bars for better visibility
    ax1.errorbar(results_df['sample_size'], results_df['mean_isoforms'],
                 yerr=results_df['std_isoforms'], fmt='none', color='darkred',
                 alpha=0.6, capsize=3, capthick=1, zorder=3)
    
    ax1.set_xlabel('Number of Reads Sampled', fontsize=12)
    ax1.set_ylabel('Number of Isoforms Detected', fontsize=12)
    ax1.set_title('Isoform Detection Saturation Curve', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Efficiency curve (isoforms per read)
    efficiency = results_df['mean_isoforms'] / results_df['sample_size']
    efficiency_std = results_df['std_isoforms'] / results_df['sample_size']
    
    # Error bands for efficiency
    eff_lower = efficiency - efficiency_std
    eff_upper = efficiency + efficiency_std
    eff_lower = np.maximum(eff_lower, 0)  # Ensure non-negative
    
    ax2.fill_between(results_df['sample_size'], eff_lower, eff_upper,
                     alpha=0.4, color='lightblue', label='±1 SD', zorder=1)
    
    ax2.plot(results_df['sample_size'], efficiency, 's-', 
             linewidth=2, markersize=6, color='blue', 
             label='Detection efficiency', zorder=2)
    
    ax2.errorbar(results_df['sample_size'], efficiency,
                 yerr=efficiency_std, fmt='none', color='darkblue',
                 alpha=0.6, capsize=3, capthick=1, zorder=3)
    
    ax2.set_xlabel('Number of Reads Sampled', fontsize=12)
    ax2.set_ylabel('Isoforms per Read', fontsize=12)
    ax2.set_title('Detection Efficiency vs Sample Size', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Coefficient of Variation (CV) to show variability
    ax3.plot(results_df['sample_size'], results_df['cv_isoforms'] * 100, 
             '^-', linewidth=2, markersize=6, color='green', 
             label='Coefficient of Variation (%)')
    
    ax3.set_xlabel('Number of Reads Sampled', fontsize=12)
    ax3.set_ylabel('Coefficient of Variation (%)', fontsize=12)
    ax3.set_title('Variability vs Sample Size', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'{output_prefix}_saturation_curve_v2.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_saturation_curve_v2.pdf', bbox_inches='tight')
    print(f"Plots saved as {output_prefix}_saturation_curve_v2.png and .pdf")
    
    plt.show()

def save_results(results_df, output_prefix):
    """
    Save results to file
    """
    output_file = f'{output_prefix}_saturation_results_v2.txt'
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
    parser = argparse.ArgumentParser(description='Generate isoform saturation curves (Improved Version)')
    parser.add_argument('input_file', help='Input file with isoform and reads data')
    parser.add_argument('output_prefix', help='Output prefix for results and plots')
    parser.add_argument('--isoform-col', default='isoform', help='Column name for isoform IDs (default: isoform)')
    parser.add_argument('--reads-col', default='FL.BioSample_6', help='Column name for read counts (default: FL.BioSample_6)')
    parser.add_argument('--replicates', type=int, default=15, help='Number of replicates for each sample size (default: 15)')
    parser.add_argument('--max-samples', type=int, default=20, help='Maximum number of sample sizes to test (default: 20)')
    
    args = parser.parse_args()
    
    print("=== Isoform Saturation Curve Analysis (v2) ===")
    print(f"Input file: {args.input_file}")
    print(f"Output prefix: {args.output_prefix}")
    print(f"Isoform column: {args.isoform_col}")
    print(f"Reads column: {args.reads_col}")
    print(f"Replicates: {args.replicates}")
    
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
    
    # Debug: Print some statistics about std values
    print("\n=== DEBUG: Standard Deviation Analysis ===")
    print(f"Mean std value: {results_df['std_isoforms'].mean():.2f}")
    print(f"Max std value: {results_df['std_isoforms'].max():.2f}")
    print(f"Min std value: {results_df['std_isoforms'].min():.2f}")
    print(f"Std as % of mean (avg): {(results_df['std_isoforms'] / results_df['mean_isoforms'] * 100).mean():.2f}%")
    print(f"Mean CV: {results_df['cv_isoforms'].mean() * 100:.2f}%")
    
    # Check if std values are too small to be visible
    if results_df['std_isoforms'].max() < 1:
        print("⚠ WARNING: Standard deviation values are very small (< 1), error bands may not be visible")
        print("Consider increasing --replicates parameter for more variation")
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main()