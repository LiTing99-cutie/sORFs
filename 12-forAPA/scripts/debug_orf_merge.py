import pandas as pd
import sys
import os

# ================= 路径配置 =================
# 请确保路径绝对正确
FILE1_PATH = '/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251119_new_db_search_res_add_public/results/augment_orf_table/augmented.tsv'
FILE2_PATH = '/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251119_new_db_search_res_add_public/results/augment_orf_table/sub.genepred'
FILE3_PATH = '/home/user/data2/lit/project/ZNF271/02-APA-1/annotation/gencode.v41.basic.annotation.gpe'

# 输出的调试文件路径
DEBUG_LOADED_PATH = 'debug_loaded_orfs.tsv'
DEBUG_FAILURE_PATH = 'debug_match_failures.tsv'

# ================= 辅助逻辑 =================

def get_introns_strictly_within_cds(exon_starts, exon_ends, cds_start, cds_end):
    introns = set()
    for i in range(len(exon_starts) - 1):
        intron_s = exon_ends[i]
        intron_e = exon_starts[i+1]
        # 内含子完全位于 CDS 范围内
        if cds_start <= intron_s and cds_end >= intron_e:
            introns.add(f"{intron_s}-{intron_e}") # 用字符串存储方便打印
    return introns

def is_coord_inside_exons(coord, exon_starts, exon_ends, is_end=False):
    for s, e in zip(exon_starts, exon_ends):
        if is_end:
            if s < coord <= e: return True
        else:
            if s <= coord < e: return True
    return False

# ================= 主程序 =================

def main():
    print("=== STEP 1: LOADING DATA & GENERATING 'debug_loaded_orfs.tsv' ===")
    
    # 1. Load File 1
    df1 = pd.read_csv(FILE1_PATH, sep='\t')
    mask = (df1['ORF_type'] == 'noncoding') & \
           ((df1['Gene_type'] == 'lncRNA') | (df1['Gene_type'].str.contains('pseudogene', case=False, na=False)))
    filtered = df1[mask].sort_values('ORF_length', ascending=False).drop_duplicates(subset=['Isoform_id'])
    
    gene_map = {}
    target_ids = set()
    
    # Debug List 1
    loaded_log = []
    
    for _, row in filtered.iterrows():
        sym = str(row['Geneid']).strip()
        oid = str(row['ORF_id']).strip()
        
        if sym not in gene_map: gene_map[sym] = []
        gene_map[sym].append({'id': oid, 'len': row['ORF_length']})
        target_ids.add(oid)
        
        # Log basic info
        loaded_log.append({
            'Gene_Symbol': sym,
            'ORF_ID': oid,
            'Length': row['ORF_length'],
            'Source': 'File 1'
        })

    # 2. Load File 2
    orf_coords = {}
    with open(FILE2_PATH, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 10: continue
            curr_id = parts[0].strip()
            
            if curr_id in target_ids:
                try:
                    c_s, c_e = int(parts[6]), int(parts[7])
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
                    
                    # Update Log
                    # 找到对应的 log entry 补充坐标信息 (效率较低但为了debug无所谓)
                    for entry in loaded_log:
                        if entry['ORF_ID'] == curr_id:
                            entry['Chrom'] = parts[1]
                            entry['CDS_Start'] = c_s
                            entry['CDS_End'] = c_e
                            entry['Found_In_File2'] = 'YES'
                except:
                    continue

    # 保存第一个中间文件
    pd.DataFrame(loaded_log).to_csv(DEBUG_LOADED_PATH, sep='\t', index=False)
    print(f"-> Saved loaded ORFs check to: {DEBUG_LOADED_PATH}")
    print(f"-> Please check this file: Are 'Chrom' and 'CDS_Start' populated? If empty, ID mismatch between File 1 & 2.")

    # ---------------------------------------------------------
    print("\n=== STEP 2: SCANNING GENCODE & GENERATING 'debug_match_failures.tsv' ===")
    
    failure_log = []
    success_count = 0
    
    with open(FILE3_PATH, 'r') as f:
        for line_idx, line in enumerate(f):
            parts = line.strip().split('\t')
            if len(parts) < 15: continue
            
            # Index 11 Check
            symbol = parts[11].strip()
            
            # 只记录那些我们在 File 1 里有的基因
            if symbol in gene_map:
                transcript_id = parts[0]
                g_chrom = parts[1]
                g_strand = parts[2]
                try:
                    g_ex_starts = [int(x) for x in parts[8].split(',') if x]
                    g_ex_ends = [int(x) for x in parts[9].split(',') if x]
                except: continue

                # 遍历所有候选 ORF
                for cand in gene_map[symbol]:
                    oid = cand['id']
                    
                    # 如果 File 2 没找到坐标，记录
                    if oid not in orf_coords:
                        failure_log.append({
                            'Gene': symbol, 'Transcript': transcript_id, 'ORF': oid,
                            'Reason': 'Coordinates_Missing', 'Details': 'ID not found in File 2'
                        })
                        continue
                    
                    s_data = orf_coords[oid]
                    
                    # 1. Chrom Check
                    if s_data['chrom'] != g_chrom:
                        failure_log.append({
                            'Gene': symbol, 'Transcript': transcript_id, 'ORF': oid,
                            'Reason': 'Chrom_Mismatch',
                            'Details': f"ORF={s_data['chrom']} vs GENCODE={g_chrom}"
                        })
                        continue
                        
                    # 2. Strand Check
                    if s_data['strand'] != g_strand:
                        failure_log.append({
                            'Gene': symbol, 'Transcript': transcript_id, 'ORF': oid,
                            'Reason': 'Strand_Mismatch',
                            'Details': f"ORF={s_data['strand']} vs GENCODE={g_strand}"
                        })
                        continue
                        
                    # 3. Start Check
                    if not is_coord_inside_exons(s_data['cdsStart'], g_ex_starts, g_ex_ends, is_end=False):
                        # 详细记录 GENCODE 外显子范围，方便肉眼看
                        exons_str = ",".join([f"{s}-{e}" for s,e in zip(g_ex_starts, g_ex_ends)])
                        failure_log.append({
                            'Gene': symbol, 'Transcript': transcript_id, 'ORF': oid,
                            'Reason': 'Start_Out_Of_Bounds',
                            'Details': f"Start {s_data['cdsStart']} not in {exons_str}"
                        })
                        continue

                    # 4. End Check
                    if not is_coord_inside_exons(s_data['cdsEnd'], g_ex_starts, g_ex_ends, is_end=True):
                        exons_str = ",".join([f"{s}-{e}" for s,e in zip(g_ex_starts, g_ex_ends)])
                        failure_log.append({
                            'Gene': symbol, 'Transcript': transcript_id, 'ORF': oid,
                            'Reason': 'End_Out_Of_Bounds',
                            'Details': f"End {s_data['cdsEnd']} not in {exons_str}"
                        })
                        continue
                        
                    # 5. Intron Check
                    g_introns = get_introns_strictly_within_cds(g_ex_starts, g_ex_ends, s_data['cdsStart'], s_data['cdsEnd'])
                    if s_data['req_introns'] != g_introns:
                        failure_log.append({
                            'Gene': symbol, 'Transcript': transcript_id, 'ORF': oid,
                            'Reason': 'Intron_Mismatch',
                            'Details': f"ORF_Needs={s_data['req_introns']} vs GENCODE_Has={g_introns}"
                        })
                        continue

                    # SUCCESS CASE
                    success_count += 1
                    # (Diagnosis mode doesn't write output GPE, just counting)
            
    # 保存失败日志
    df_fail = pd.DataFrame(failure_log)
    if not df_fail.empty:
        df_fail.to_csv(DEBUG_FAILURE_PATH, sep='\t', index=False)
        print(f"-> Saved failure details to: {DEBUG_FAILURE_PATH}")
    else:
        print("-> Strange... No failures logged? (Did any gene symbols match?)")

    print(f"\nTotal Success Matches Identified: {success_count}")
    print("Please open 'debug_match_failures.tsv' and filter for 'ZNF271P'.")

if __name__ == "__main__":
    main()