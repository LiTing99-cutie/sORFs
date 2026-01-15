#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, csv, random

def load_expr_table(path):
    expr = {}
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        iso_col = 0
        tpm_col = next((i for i,h in enumerate(header) if "FL_TPM" in h), len(header)-1)
        for row in reader:
            if not row: continue
            iso = row[iso_col].strip()
            if not iso: continue
            try: expr[iso] = float(row[tpm_col])
            except: expr[iso] = 0.0
    return expr

def parse_id_info(s):
    s = s.strip()
    if s.startswith(("sp|","tr|")):
        return s, True, "ATG", None
    parts = s.split("|")
    if len(parts) >= 5:
        orf_type = parts[-2].strip().lower()
        start_codon = parts[-1].strip().strip(",").upper() or "NA"
        iso = parts[0].split(":",1)[0]  # PB.x.x
        return s, (orf_type=="canonical"), start_codon, iso
    return s, False, "NA", None

def pick_representative(ids, expr_map, rng):
    infos = [parse_id_info(x) for x in ids]
    cand = [inf for inf in infos if inf[1] is True] or infos
    def score(inf):
        _, _, start_codon, iso = inf
        tpm = expr_map.get(iso, 0.0) if iso else 0.0
        start_rank = 0 if start_codon == "ATG" else 1
        return (-tpm, start_rank, rng.random())  # 仅用于打破完全并列
    return min(cand, key=score)

def main():
    ap = argparse.ArgumentParser(description="代表性ID：canonical最高；否则按表达量→ATG→随机。并输出逗号分隔后的第一个ID。")
    ap.add_argument("-i","--input", required=True, help="两列：count\\tids（逗号分隔）")
    ap.add_argument("-e","--expr",  required=True, help="isoform表达表（含FL_TPM列）")
    ap.add_argument("-o","--output",required=True, help="输出TSV")
    ap.add_argument("--seed", type=int, default=20250910, help="随机种子（仅用于打破代表性并列）")
    args = ap.parse_args()

    expr_map = load_expr_table(args.expr)
    rng = random.Random(args.seed)

    with open(args.output,"w",newline="") as w, open(args.input) as f:
        writer = csv.writer(w, delimiter="\t")
        writer.writerow([
            "group_size","chosen_id","first_id_in_line",
            "all_ids","chosen_is_canonical","chosen_start_codon","chosen_expr"
        ])
        for line in f:
            line=line.strip()
            if not line: continue
            cols=line.split("\t",1)
            if len(cols)<2: continue
            try: n=int(cols[0])
            except ValueError: continue
            ids=[x.strip() for x in cols[1].split(",") if x.strip()]
            if not ids: continue

            # 代表性ID
            full_id, is_canon, start_codon, iso = pick_representative(ids, expr_map, rng)
            tpm = expr_map.get(iso, 0.0) if iso else 0.0

            # 逗号分隔后的第一个ID
            first_id = ids[0]

            writer.writerow([
                n, full_id, first_id,
                ", ".join(ids), "Yes" if is_canon else "No",
                start_codon, f"{tpm:.6f}"
            ])

if __name__ == "__main__":
    main()
