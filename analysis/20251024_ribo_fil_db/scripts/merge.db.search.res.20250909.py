#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import re
from pathlib import Path
import argparse

def load_msf(path, source):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    df = df[['Spectrum', 'Peptide']]
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L', regex=False)
    df['Source'] = source
    return df

def load_pfind(path, source):
    df = pd.read_csv(path, sep='\t', header=0, low_memory=False)
    df = df[['File_Name', 'Sequence']]
    df.columns = ['Spectrum', 'Peptide']
    df['Peptide_I_L_equal'] = df['Peptide'].str.replace('I', 'L', regex=False)
    df['Source'] = source
    return df

def format_msf_spectrum(psm):
    split = psm['Spectrum'].str.split('.', n=3, expand=True)
    prefix = split[0]
    scan_number = split[1].str.lstrip('0')
    charge = psm['Charge'] if 'Charge' in psm.columns else split[3]
    psm['Spectrum'] = prefix + '.' + scan_number + '.' + scan_number + '.' + charge.astype(str)
    return psm

def format_pfind_spectrum(df):
    df['Spectrum'] = df['Spectrum'].str.replace(r'_SCANS\d+$', '', regex=True)
    return df

def main():
    ap = argparse.ArgumentParser(description="Merge MSFragger/PFind closed+open PSMs.")
    ap.add_argument('--msf_closed',  type=Path, default=Path('../output/db_search/msfragger/closed/psm.tsv'))
    ap.add_argument('--msf_open',    type=Path, default=Path('../output/db_search/msfragger/open/psm.tsv'))
    ap.add_argument('--pfind_closed',type=Path, default=Path('../output/db_search/pfind/closed/pFind-Filtered.spectra'))
    ap.add_argument('--pfind_open',  type=Path, default=Path('../output/db_search/pfind/open/pFind-Filtered.spectra'))
    ap.add_argument('--out',         type=Path, default=Path('../results/final.csv'))
    args = ap.parse_args()

    # 存在性检查（更友好地报错）
    for p in [args.msf_closed, args.msf_open, args.pfind_closed, args.pfind_open]:
        if not p.exists():
            raise FileNotFoundError(f"找不到文件: {p}")

    # DB search
    msf_closed = format_msf_spectrum(load_msf(args.msf_closed, 'msfragger_closed'))
    msf_open   = format_msf_spectrum(load_msf(args.msf_open,   'msfragger_open'))
    pfind_closed = format_pfind_spectrum(load_pfind(args.pfind_closed, 'pfind_closed'))
    pfind_open   = format_pfind_spectrum(load_pfind(args.pfind_open,   'pfind_open'))

    # msfragger 以 closed 为主
    msf_open_extra = msf_open[~msf_open['Spectrum'].isin(msf_closed['Spectrum'])]
    msfragger_all = pd.concat([msf_closed, msf_open_extra], ignore_index=True)

    # pfind 以 closed 为主
    pfind_open_extra = pfind_open[~pfind_open['Spectrum'].isin(pfind_closed['Spectrum'])]
    pfind_all = pd.concat([pfind_closed, pfind_open_extra], ignore_index=True)

    # 以 msfragger 为主，pfind 补充
    pfind_extra = pfind_all[~pfind_all['Spectrum'].isin(msfragger_all['Spectrum'])]
    final = pd.concat([msfragger_all, pfind_extra], ignore_index=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.out, index=False)
    print(f"Done. Output -> {args.out}")

if __name__ == "__main__":
    main()
