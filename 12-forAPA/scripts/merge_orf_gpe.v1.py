import pandas as pd
import sys

# ================= 路径配置 =================
FILE1_PATH = '/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251119_new_db_search_res_add_public/results/augment_orf_table/augmented.tsv'
FILE2_PATH = '/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251119_new_db_search_res_add_public/results/augment_orf_table/sub.genepred'
FILE3_PATH = '/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gpe'
OUTPUT_PATH = '../processed/gencode.v41.basic.annotation.with_new_orfs.gpe'

# ================= 辅助函数 =================

def get_introns_strictly_within_cds(exon_starts, exon_ends, cds_start, cds_end):
    """
    提取完全位于 CDS 范围内的内含子。
    用于验证 CDS 内部的剪接结构是否一致。
    """
    introns = set()
    for i in range(len(exon_starts) - 1):
        intron_s = exon_ends[i]
        intron_e = exon_starts[i+1]
        if cds_start <= intron_s and cds_end >= intron_e:
            introns.add((intron_s, intron_e))
    return introns

def is_coord_inside_exons(coord, exon_starts, exon_ends, is_end_coord=False):
    """
    检查坐标点是否落在任意一个外显子区间内。
    """
    for s, e in zip(exon_starts, exon_ends):
        if is_end_coord:
            if s < coord <= e: return True
        else:
            if s <= coord < e: return True
    return False

def calculate_frames(exon_starts, exon_ends, cds_start, cds_end, strand):
    """重算 Frames"""
    cds_lengths = []
    has_cds = []
    for s, e in zip(exon_starts, exon_ends):
        o_s = max(s, cds_start)
        o_e = min(e, cds_end)
        if o_s < o_e:
            cds_lengths.append(o_e - o_s)
            has_cds.append(True)
        else:
            cds_lengths.append(0)
            has_cds.append(False)

    final_frames = [-1] * len(exon_starts)
    if strand == '+':
        curr = 0
        for i in range(len(exon_starts)):
            if has_cds[i]:
                final_frames[i] = curr
                curr = (curr + cds_lengths[i]) % 3
    else:
        curr = 0
        for i in range(len(exon_starts) - 1, -1, -1):
            if has_cds[i]:
                final_frames[i] = curr
                curr = (curr + cds_lengths[i]) % 3

    return ",".join(map(str, final_frames)) + ","

# ================= 主程序 =================

def main():
    print("Step 1: Loading Augmented ORF Table (File 1)...")
    try:
        df1 = pd.read_csv(FILE1_PATH, sep='\t')
    except Exception as e:
        print(f"Error: {e}")
        return

    # 1.1 筛选
    mask = (df1['ORF_type'] == 'noncoding') & \
           ((df1['Gene_type'] == 'lncRNA') | (df1['Gene_type'].str.contains('pseudogene', case=False, na=False)))
    filtered = df1[mask].sort_values('ORF_length', ascending=False).drop_duplicates(subset=['Isoform_id'])
    
    # 1.2 建立映射 (Symbol -> List of ORFs)
    gene_map = {}
    target_ids = set()
    for _, row in filtered.iterrows():
        sym = str(row['Geneid']).strip()
        oid = str(row['ORF_id']).strip()
        if sym not in gene_map: gene_map[sym] = []
        gene_map[sym].append({'id': oid, 'len': row['ORF_length']})
        target_ids.add(oid)
        
    print(f"  Filtered ORFs: {len(filtered)}")
    print(f"  Mapped {len(target_ids)} unique ORFs to {len(gene_map)} genes.")

    # ---------------------------------------------------------
    print("\nStep 2: Loading Coordinates (File 2)...")
    orf_coords = {}
    
    with open(FILE2_PATH, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 10: continue
            
            curr_id = parts[0].strip()
            
            if curr_id in target_ids:
                try:
                    c_s = int(parts[5]) 
                    c_e = int(parts[6])
                    e_s = [int(x) for x in parts[8].split(',') if x]
                    e_e = [int(x) for x in parts[9].split(',') if x]
                    
                    req_introns = get_introns_strictly_within_cds(e_s, e_e, c_s, c_e)
                    
                    orf_coords[curr_id] = {
                        'chrom': parts[1],
                        'strand': parts[2],
                        'cdsStart': c_s,
                        'cdsEnd': c_e,
                        'req_introns': req_introns
                    }
                except ValueError: continue
                    
    print(f"  Loaded coordinates for {len(orf_coords)} ORFs.")

    # ---------------------------------------------------------
    print("\nStep 3: Merging into GENCODE...")
    
    cnt_total = 0
    cnt_updated = 0
    cnt_skip = 0
    
    # === 新增：追踪被整合的ORF ===
    integrated_orfs = set()
    integrated_genes = set()
    
    with open(FILE3_PATH, 'r') as f_in, open(OUTPUT_PATH, 'w') as f_out:
        for line in f_in:
            cnt_total += 1
            parts = line.strip().split('\t')
            if len(parts) < 15:
                f_out.write(line); continue
            
            symbol = parts[11].strip()
            
            if symbol in gene_map:
                try:
                    g_chrom = parts[1]
                    g_strand = parts[2]
                    g_ex_starts = [int(x) for x in parts[8].split(',') if x]
                    g_ex_ends = [int(x) for x in parts[9].split(',') if x]
                except:
                    f_out.write(line); continue
                
                valid_candidates = []
                
                for cand in gene_map[symbol]:
                    oid = cand['id']
                    if oid not in orf_coords: continue
                    s_data = orf_coords[oid]
                    
                    if s_data['chrom'] != g_chrom or s_data['strand'] != g_strand: continue
                    
                    if not is_coord_inside_exons(s_data['cdsStart'], g_ex_starts, g_ex_ends, is_end_coord=False): continue
                    if not is_coord_inside_exons(s_data['cdsEnd'], g_ex_starts, g_ex_ends, is_end_coord=True): continue
                        
                    g_introns = get_introns_strictly_within_cds(g_ex_starts, g_ex_ends, s_data['cdsStart'], s_data['cdsEnd'])
                    
                    if s_data['req_introns'] == g_introns:
                        valid_candidates.append(cand)

                if valid_candidates:
                    best = max(valid_candidates, key=lambda x: x['len'])
                    best_data = orf_coords[best['id']]
                    
                    parts[6] = str(best_data['cdsStart'])
                    parts[7] = str(best_data['cdsEnd'])
                    parts[12] = 'cmpl'
                    parts[13] = 'cmpl'
                    new_frames = calculate_frames(g_ex_starts, g_ex_ends, best_data['cdsStart'], best_data['cdsEnd'], g_strand)
                    if len(parts) > 14: parts[14] = new_frames
                    else: parts.append(new_frames)
                    
                    f_out.write('\t'.join(parts) + '\n')
                    cnt_updated += 1
                    
                    # === 新增：记录被整合的ORF和Gene ===
                    integrated_orfs.add(best['id'])
                    integrated_genes.add(symbol)
                else:
                    cnt_skip += 1
                    f_out.write(line)
            else:
                f_out.write(line)

    # === 新增：详细统计输出 ===
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Input Statistics:")
    print(f"  - Candidate ORFs (after filtering): {len(target_ids)}")
    print(f"  - Candidate Genes: {len(gene_map)}")
    print(f"  - ORFs with valid coordinates: {len(orf_coords)}")
    print()
    print(f"GENCODE Integration:")
    print(f"  - Total transcripts scanned: {cnt_total}")
    print(f"  - Transcripts updated: {cnt_updated}")
    print(f"  - Transcripts skipped (structure mismatch): {cnt_skip}")
    print()
    print(f"ORF Integration:")
    print(f"  - Unique ORFs integrated: {len(integrated_orfs)}")
    print(f"  - Unique Genes with ORF annotation: {len(integrated_genes)}")
    print(f"  - ORFs not integrated: {len(target_ids) - len(integrated_orfs)}")
    print()
    print(f"Output saved to: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()