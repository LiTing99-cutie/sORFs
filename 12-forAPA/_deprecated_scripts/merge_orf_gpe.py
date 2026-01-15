import pandas as pd
import sys

# ================= 路径配置 =================
FILE1_PATH = '/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251119_new_db_search_res_add_public/results/augment_orf_table/augmented.tsv'
FILE2_PATH = '/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251119_new_db_search_res_add_public/results/augment_orf_table/sub.genepred'
FILE3_PATH = '/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gpe'
OUTPUT_PATH = 'gencode.v41.basic.annotation.with_new_orfs.gpe'

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
        # 只有当内含子被 CDS 范围完全覆盖时（即 CDS 跨越了这个内含子）
        if cds_start <= intron_s and cds_end >= intron_e:
            introns.add((intron_s, intron_e))
    return introns

def is_coord_inside_exons(coord, exon_starts, exon_ends, is_end_coord=False):
    """
    检查坐标点是否落在任意一个外显子区间内。
    """
    for s, e in zip(exon_starts, exon_ends):
        if is_end_coord:
            # genePred End 是不包含的(exclusive)，所以允许 coord == exon_end
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
        
    print(f"  Mapped {len(target_ids)} unique ORFs to {len(gene_map)} genes.")

    # ---------------------------------------------------------
    print("Step 2: Loading Coordinates (File 2) [CORRECTED INDICES]...")
    orf_coords = {}
    
    with open(FILE2_PATH, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 10: continue
            
            curr_id = parts[0].strip() # name (Index 0)
            
            if curr_id in target_ids:
                try:
                    # === 这里是本次修正的关键 ===
                    # 0: name
                    # 1: chrom
                    # 2: strand
                    # 3: txStart
                    # 4: txEnd
                    # 5: cdsStart  <--- 修正
                    # 6: cdsEnd    <--- 修正
                    # 7: exonCount
                    # 8: exonStarts
                    # 9: exonEnds
                    
                    c_s = int(parts[5]) 
                    c_e = int(parts[6])
                    e_s = [int(x) for x in parts[8].split(',') if x]
                    e_e = [int(x) for x in parts[9].split(',') if x]
                    
                    # 预计算所需的内含子结构
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
    print("Step 3: Merging into GENCODE (Option B)...")
    
    cnt_total = 0
    cnt_updated = 0
    cnt_skip = 0
    
    with open(FILE3_PATH, 'r') as f_in, open(OUTPUT_PATH, 'w') as f_out:
        for line in f_in:
            cnt_total += 1
            parts = line.strip().split('\t')
            if len(parts) < 15:
                f_out.write(line); continue
            
            # File 3: Gene Symbol 在 Index 11
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
                    
                    # 1. 染色体与链必须一致
                    if s_data['chrom'] != g_chrom or s_data['strand'] != g_strand: continue
                    
                    # 2. 边界检查 (Option B: 只要 CDS 落在外显子内即可，忽略 UTR)
                    if not is_coord_inside_exons(s_data['cdsStart'], g_ex_starts, g_ex_ends, is_end_coord=False): continue
                    if not is_coord_inside_exons(s_data['cdsEnd'], g_ex_starts, g_ex_ends, is_end_coord=True): continue
                        
                    # 3. 内含子结构检查
                    # 获取 GENCODE 在 CDS 范围内的内含子
                    g_introns = get_introns_strictly_within_cds(g_ex_starts, g_ex_ends, s_data['cdsStart'], s_data['cdsEnd'])
                    
                    # 必须一致
                    if s_data['req_introns'] == g_introns:
                        valid_candidates.append(cand)

                if valid_candidates:
                    # 选最长
                    best = max(valid_candidates, key=lambda x: x['len'])
                    best_data = orf_coords[best['id']]
                    
                    # 更新 CDS 坐标 (Index 6, 7)
                    parts[6] = str(best_data['cdsStart'])
                    parts[7] = str(best_data['cdsEnd'])
                    # 更新状态 (Index 12, 13)
                    parts[12] = 'cmpl'
                    parts[13] = 'cmpl'
                    # 更新 Frames (Index 14)
                    new_frames = calculate_frames(g_ex_starts, g_ex_ends, best_data['cdsStart'], best_data['cdsEnd'], g_strand)
                    if len(parts) > 14: parts[14] = new_frames
                    else: parts.append(new_frames)
                    
                    f_out.write('\t'.join(parts) + '\n')
                    cnt_updated += 1
                else:
                    cnt_skip += 1
                    f_out.write(line)
            else:
                f_out.write(line)

    print("-" * 30)
    print(f"Total GENCODE transcripts scanned: {cnt_total}")
    print(f"Successfully Updated: {cnt_updated}")
    print(f"Skipped (Structure Mismatch): {cnt_skip}")
    print(f"Output saved to: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()