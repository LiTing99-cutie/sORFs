#!/usr/bin/env python3
"""
Isoform Saturation Curve Data Generator
Version: v1_20250710
Date: 2025-07-10

This script generates data for isoform saturation curves by randomly sampling
different numbers of reads and counting the number of isoforms detected.
Data is saved for R plotting with ggplot2.

Usage:
    python saturation_curve_data_generator_v1_20250710.py input_file output_prefix
"""

import pandas as pd
import numpy as np
import argparse
import random
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

def prepare_reads_data(df, isoform_col, gene_col, reads_col):
    """
    Prepare reads data for sampling
    """
    # Create lists where each read is represented by its isoform and gene
    reads_list = []
    genes_list = []
    
    for _, row in df.iterrows():
        isoform = row[isoform_col]
        gene = row[gene_col]
        reads_count = row[reads_col]
        
        # Handle NaN or missing values
        if pd.isna(reads_count) or reads_count <= 0:
            continue
            
        # Add the isoform and gene to the lists reads_count times
        reads_list.extend([isoform] * int(reads_count))
        genes_list.extend([gene] * int(reads_count))
    
    print(f"Total reads: {len(reads_list)}")
    print(f"Unique isoforms: {len(set(reads_list))}")
    print(f"Unique genes: {len(set(genes_list))}")
    
    return reads_list, genes_list

def calculate_saturation_curve(reads_list, genes_list, sample_sizes, n_replicates=10):
    """
    Calculate saturation curve by sampling different numbers of reads
    """
    results = []
    
    for sample_size in sample_sizes:
        sample_isoforms = []
        sample_genes = []
        
        for replicate in range(n_replicates):
            # Randomly sample reads
            if sample_size >= len(reads_list):
                sampled_reads = reads_list
                sampled_genes = genes_list
            else:
                indices = random.sample(range(len(reads_list)), sample_size)
                sampled_reads = [reads_list[i] for i in indices]
                sampled_genes = [genes_list[i] for i in indices]
            
            # Count unique isoforms and genes
            unique_isoforms = len(set(sampled_reads))
            unique_genes = len(set(sampled_genes))
            
            sample_isoforms.append(unique_isoforms)
            sample_genes.append(unique_genes)
        
        # Calculate statistics for isoforms
        mean_isoforms = np.mean(sample_isoforms)
        std_isoforms = np.std(sample_isoforms)
        
        # Calculate statistics for genes
        mean_genes = np.mean(sample_genes)
        std_genes = np.std(sample_genes)
        
        results.append({
            'sample_size': sample_size,
            'mean_isoforms': mean_isoforms,
            'std_isoforms': std_isoforms,
            'min_isoforms': np.min(sample_isoforms),
            'max_isoforms': np.max(sample_isoforms),
            'cv_isoforms': std_isoforms / mean_isoforms if mean_isoforms > 0 else 0,
            'mean_genes': mean_genes,
            'std_genes': std_genes,
            'min_genes': np.min(sample_genes),
            'max_genes': np.max(sample_genes),
            'cv_genes': std_genes / mean_genes if mean_genes > 0 else 0,
            'replicates': n_replicates
        })
    
    return pd.DataFrame(results)

def save_data_for_r(results_df, output_prefix):
    """
    Save data in formats suitable for R analysis
    """
    # Save main results
    output_file = f'{output_prefix}_saturation_data.txt'
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Main data saved to {output_file}")
    
    # Create data for R plotting with proper column names
    r_data = results_df.copy()
    r_data.columns = ['sample_size', 'mean_isoforms', 'std_isoforms', 
                      'min_isoforms', 'max_isoforms', 'cv_isoforms',
                      'mean_genes', 'std_genes', 'min_genes', 'max_genes', 'cv_genes', 'replicates']
    
    # Add efficiency columns
    r_data['efficiency_isoforms'] = r_data['mean_isoforms'] / r_data['sample_size']
    r_data['efficiency_isoforms_std'] = r_data['std_isoforms'] / r_data['sample_size']
    r_data['efficiency_genes'] = r_data['mean_genes'] / r_data['sample_size']
    r_data['efficiency_genes_std'] = r_data['std_genes'] / r_data['sample_size']
    
    # Add confidence intervals for isoforms
    r_data['isoforms_upper'] = r_data['mean_isoforms'] + r_data['std_isoforms']
    r_data['isoforms_lower'] = np.maximum(r_data['mean_isoforms'] - r_data['std_isoforms'], 0)
    r_data['efficiency_isoforms_upper'] = r_data['efficiency_isoforms'] + r_data['efficiency_isoforms_std']
    r_data['efficiency_isoforms_lower'] = np.maximum(r_data['efficiency_isoforms'] - r_data['efficiency_isoforms_std'], 0)
    
    # Add confidence intervals for genes
    r_data['genes_upper'] = r_data['mean_genes'] + r_data['std_genes']
    r_data['genes_lower'] = np.maximum(r_data['mean_genes'] - r_data['std_genes'], 0)
    r_data['efficiency_genes_upper'] = r_data['efficiency_genes'] + r_data['efficiency_genes_std']
    r_data['efficiency_genes_lower'] = np.maximum(r_data['efficiency_genes'] - r_data['efficiency_genes_std'], 0)
    
    # Save R-ready data
    r_output_file = f'{output_prefix}_for_r_plotting.txt'
    r_data.to_csv(r_output_file, sep='\t', index=False)
    print(f"R-ready data saved to {r_output_file}")
    
    # Print summary
    print("\n=== SATURATION ANALYSIS SUMMARY ===")
    print(f"Total reads analyzed: {results_df['sample_size'].max()}")
    print(f"Maximum isoforms detected: {results_df['mean_isoforms'].max():.1f}")
    print(f"Detection efficiency at max: {results_df['mean_isoforms'].max() / results_df['sample_size'].max():.4f} isoforms/read")
    
    # Check if saturation is reached
    if len(results_df) >= 2:
        # Use last two points to calculate slope
        last_two = results_df.tail(2)
        slope = (last_two['mean_isoforms'].iloc[-1] - last_two['mean_isoforms'].iloc[-2]) / \
                (last_two['sample_size'].iloc[-1] - last_two['sample_size'].iloc[-2])
        print(f"Slope between last two points: {slope:.6f}")
        if slope < 0.01:
            print("✓ Saturation appears to be reached")
        else:
            print("⚠ Saturation may not be reached yet")
    
    # Debug: Print some statistics about std values
    print("\n=== DEBUG: Standard Deviation Analysis ===")
    print(f"Mean std value: {results_df['std_isoforms'].mean():.2f}")
    print(f"Max std value: {results_df['std_isoforms'].max():.2f}")
    print(f"Min std value: {results_df['std_isoforms'].min():.2f}")
    print(f"Std as % of mean (avg): {(results_df['std_isoforms'] / results_df['mean_isoforms'] * 100).mean():.2f}%")
    
    # Check if std values are too small to be visible
    if results_df['std_isoforms'].max() < 1:
        print("⚠ WARNING: Standard deviation values are very small (< 1), error bands may not be visible")
        print("Consider increasing --replicates parameter for more variation")
    
    return r_data

def main():
    parser = argparse.ArgumentParser(description='Generate isoform saturation curve data for R plotting')
    parser.add_argument('input_file', help='Input file with isoform and reads data')
    parser.add_argument('output_prefix', help='Output prefix for results and plots')
    parser.add_argument('--isoform-col', default='isoform', help='Column name for isoform IDs (default: isoform)')
    parser.add_argument('--gene-col', default='associated_gene', help='Column name for gene IDs (default: associated_gene)')
    parser.add_argument('--reads-col', default='FL.BioSample_6', help='Column name for read counts (default: FL.BioSample_6)')
    parser.add_argument('--replicates', type=int, default=15, help='Number of replicates for each sample size (default: 15)')
    parser.add_argument('--max-samples', type=int, default=20, help='Maximum number of sample sizes to test (default: 20)')
    
    args = parser.parse_args()
    
    print("=== Isoform Saturation Curve Data Generator ===")
    print(f"Input file: {args.input_file}")
    print(f"Output prefix: {args.output_prefix}")
    print(f"Isoform column: {args.isoform_col}")
    print(f"Gene column: {args.gene_col}")
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
    
    if args.gene_col not in df.columns:
        print(f"Error: Column '{args.gene_col}' not found in data")
        print(f"Available columns: {list(df.columns)}")
        return
    
    if args.reads_col not in df.columns:
        print(f"Error: Column '{args.reads_col}' not found in data")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Prepare reads data
    reads_list, genes_list = prepare_reads_data(df, args.isoform_col, args.gene_col, args.reads_col)
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
    results_df = calculate_saturation_curve(reads_list, genes_list, sample_sizes, args.replicates)
    
    # Save data for R plotting
    r_data = save_data_for_r(results_df, args.output_prefix)
    
    print("Data generation completed successfully!")
    print(f"Use the R script to create plots from: {args.output_prefix}_for_r_plotting.txt")

if __name__ == "__main__":
    main() 