#!/usr/bin/env python3
"""
Generate ORF peptide files and create input list for downstream analysis.

Usage:
    python prepare_orf_files.py <bed_file> <orfs_pep_fa> <nucl_dir> <pep_dir> <output_dir>

Arguments:
    bed_file      : BED file containing ORF coordinates (column 4 = ORF ID)
    orfs_pep_fa   : FASTA file with all ORF peptide sequences
    nucl_dir      : Directory containing nucleotide sequence files (chr*/*.nucl.fa)
    pep_dir       : Directory containing protein sequence files (chr*/*.prot.fa)
    output_dir    : Output directory for results

Example:
    python prepare_orf_files.py \\
        data/orfs.bed \\
        data/all_orfs.fa \\
        results/orfs/nucl \\
        results/orfs/prot \\
        results/output
"""

import os
import sys
import re
from pathlib import Path

def sanitize_filename(s: str) -> str:
    """Produce a filesystem-safe name from an arbitrary string."""
    return re.sub(r"[^A-Za-z0-9._-]+", "__", s)

def generate_id_mapping(bed_file, output_dir):
    """Extract unique IDs from BED file and create safe ID mapping."""
    print("Step 1: Generating ID mapping...")
    
    # Extract unique IDs from column 4
    ids = set()
    with open(bed_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    ids.add(fields[3])
    
    if not ids:
        print("  ERROR: No IDs found in BED file!")
        sys.exit(1)
    
    # Write IDs and create mapping
    id_file = output_dir / "orf.id.txt"
    mapping_file = output_dir / "id_mapping.txt"
    
    with open(id_file, 'w') as f_id, open(mapping_file, 'w') as f_map:
        for orig_id in sorted(ids):
            f_id.write(f"{orig_id}\n")
            safe_id = sanitize_filename(orig_id)
            f_map.write(f"{orig_id}\t{safe_id}\n")
    
    print(f"  - Found {len(ids)} unique ORF IDs")
    print(f"  - ID list saved to: {id_file}")
    print(f"  - Mapping saved to: {mapping_file}")
    return mapping_file

def extract_sequences(mapping_file, orfs_pep_fa, output_dir):
    """Extract sequences from FASTA based on ID mapping."""
    print("\nStep 2: Extracting sequences from FASTA...")
    
    # Read mapping
    id_map = {}
    with open(mapping_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                orig_id, safe_id = parts
                id_map[orig_id] = safe_id
    
    print(f"  - Loaded {len(id_map)} ID mappings")
    
    # Process FASTA
    output_files = {}
    current_id = None
    count = 0
    
    with open(orfs_pep_fa) as f:
        for line in f:
            if line.startswith('>'):
                # Extract sequence ID (first word after >)
                seq_id = line[1:].split()[0]
                
                if seq_id in id_map:
                    current_id = id_map[seq_id]
                    if current_id not in output_files:
                        out_file = output_dir / f"{current_id}.ORF_pep.fa"
                        output_files[current_id] = open(out_file, 'w')
                        count += 1
                    output_files[current_id].write(line)
                else:
                    current_id = None
            elif current_id:
                output_files[current_id].write(line)
    
    # Close all files
    for f in output_files.values():
        f.close()
    
    print(f"  - Extracted {count} ORF sequences")
    print(f"  - Output directory: {output_dir}")
    
    # Report unmapped IDs
    unmapped = set(id_map.values()) - set(output_files.keys())
    if unmapped:
        print(f"  WARNING: {len(unmapped)} IDs not found in FASTA")
        if len(unmapped) <= 10:
            for uid in sorted(unmapped):
                print(f"    - {uid}")
    
    return count

def create_input_list(nucl_dir, pep_dir, orf_pep_dir, output_file):
    """Create input list file with paths to all required files."""
    print("\nStep 3: Creating input list...")
    
    nucl_dir = Path(nucl_dir)
    pep_dir = Path(pep_dir)
    orf_pep_dir = Path(orf_pep_dir)
    
    if not nucl_dir.exists():
        print(f"  ERROR: nucl_dir does not exist: {nucl_dir}")
        sys.exit(1)
    if not pep_dir.exists():
        print(f"  ERROR: pep_dir does not exist: {pep_dir}")
        sys.exit(1)
    
    # Pre-build lookup dictionaries (一次性构建索引)
    print("  - Building file index...")
    
    # Index all prot files by base name
    prot_index = {}
    for prot_fa in pep_dir.rglob("*.prot.fa"):
        id_base = prot_fa.stem.replace('.prot', '')
        prot_index[id_base] = prot_fa
    
    # Index all ORF_pep files by base name
    orf_index = {}
    for orf_file in orf_pep_dir.glob("*.ORF_pep.fa"):
        id_base = orf_file.stem.replace('.ORF_pep', '')
        orf_index[id_base] = orf_file
    
    print(f"    Found {len(prot_index)} prot files")
    print(f"    Found {len(orf_index)} ORF_pep files")
    
    if output_file.exists():
        output_file.unlink()
    
    count = 0
    missing_prot = []
    missing_orf = []
    
    with open(output_file, 'w') as f:
        # Find all nucl.fa files
        nucl_files = list(nucl_dir.rglob("*.nucl.fa"))
        
        if not nucl_files:
            print(f"  WARNING: No .nucl.fa files found in {nucl_dir}")
        
        for nucl_fa in sorted(nucl_files):
            id_base = nucl_fa.stem.replace('.nucl', '')
            
            # Look up files in pre-built indexes (O(1) 查找)
            prot_fa = prot_index.get(id_base)
            orf_file = orf_index.get(id_base)
            
            # Check if all files exist
            if prot_fa is None:
                missing_prot.append(id_base)
                continue
            if orf_file is None:
                missing_orf.append(id_base)
                continue
            
            # Write to output
            f.write(f"{prot_fa}\t{orf_file}\t{nucl_fa}\n")
            count += 1
    
    print(f"  - Created {count} entries")
    print(f"  - Output file: {output_file}")
    
    # Report missing files
    if missing_prot:
        print(f"\n  WARNING: {len(missing_prot)} .prot.fa files not found")
        for m in missing_prot[:5]:
            print(f"    - {m}")
        if len(missing_prot) > 5:
            print(f"    ... and {len(missing_prot) - 5} more")
    
    if missing_orf:
        print(f"\n  WARNING: {len(missing_orf)} .ORF_pep.fa files not found")
        for m in missing_orf[:5]:
            print(f"    - {m}")
        if len(missing_orf) > 5:
            print(f"    ... and {len(missing_orf) - 5} more")
    
    return count

def main():
    # Check arguments
    if len(sys.argv) != 6:
        print(__doc__)
        sys.exit(1)
    
    # Parse arguments
    bed_file = Path(sys.argv[1])
    orfs_pep_fa = Path(sys.argv[2])
    nucl_dir = Path(sys.argv[3])
    pep_dir = Path(sys.argv[4])
    out_dir = Path(sys.argv[5])
    
    # Validate input files
    if not bed_file.exists():
        print(f"ERROR: BED file not found: {bed_file}")
        sys.exit(1)
    if not orfs_pep_fa.exists():
        print(f"ERROR: FASTA file not found: {orfs_pep_fa}")
        sys.exit(1)
    
    # Create output directories
    orf_pep_dir = out_dir / "orf_pep"
    orf_pep_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("ORF Peptide File Generator")
    print("=" * 70)
    print(f"\nInput files:")
    print(f"  BED file      : {bed_file}")
    print(f"  ORFs FASTA    : {orfs_pep_fa}")
    print(f"  Nucl dir      : {nucl_dir}")
    print(f"  Prot dir      : {pep_dir}")
    print(f"  Output dir    : {out_dir}")
    print(f"  ORF pep dir   : {orf_pep_dir}")
    print()
    
    try:
        # Step 1: Generate ID mapping
        mapping_file = generate_id_mapping(bed_file, orf_pep_dir)
        
        # Step 2: Extract sequences
        extract_sequences(mapping_file, orfs_pep_fa, orf_pep_dir)
        
        # Step 3: Create input list
        output_file = out_dir / "input_list.txt"
        create_input_list(nucl_dir, pep_dir, orf_pep_dir, output_file)
        
        print("\n" + "=" * 70)
        print("✓ Processing complete!")
        print("=" * 70)
        print(f"\nOutput files:")
        print(f"  ORF peptides  : {orf_pep_dir}/")
        print(f"  Input list    : {output_file}")
        print(f"  ID mapping    : {mapping_file}")
        print()
        
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()