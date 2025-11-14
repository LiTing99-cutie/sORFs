#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, logging
import pandas as pd
from pathlib import Path
import numpy as np

TRUE, FALSE = "TRUE", "FALSE"

def setup_logger(verbosity:int):
    lvl = logging.WARNING if verbosity==0 else (logging.INFO if verbosity==1 else logging.DEBUG)
    logging.basicConfig(
        level=lvl,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger("merge_orf_flags")

def main():
    ap = argparse.ArgumentParser(description="orfs_3_ways 合并 + 统计（含进度日志）")
    ap.add_argument("--orfs-3-ways", required=True, help="列：ORF_id,Seq,Software,ORF_id_custom（默认TSV）")
    ap.add_argument("--sep3ways", default="\t", help="orfs_3_ways 分隔符（默认 \\t）")
    ap.add_argument("--chunksize", type=int, default=2_000_000, help="orfs_3_ways 按行分块大小")
    ap.add_argument("--in-with-custom", required=True, help="主表（含 ORF_id_custom）")
    ap.add_argument("--out-with-flags", required=True, help="输出表（新增 Is_*）")
    ap.add_argument("--out-stats", default=None, help="可选：导出统计（.tsv 或 .xlsx）")
    ap.add_argument("-v","--verbose", action="count", default=1, help="-v 普通日志，-vv 调试日志，默认 INFO")
    ap.add_argument("--tqdm", action="store_true", help="使用进度条（需安装 tqdm）")
    args = ap.parse_args()
    log = setup_logger(args.verbose)

    # 优化1: 读取主表时只加载必要列，使用更快的数据类型
    log.info("读取主表：%s", args.in_with_custom)
    main_df = pd.read_csv(
        args.in_with_custom, 
        sep="\t", 
        dtype=str, 
        quoting=3, 
        compression="infer",
        engine='c'  # 使用C引擎加速
    )
    if "ORF_id_custom" not in main_df.columns:
        raise SystemExit("[ERR] 主表缺少 ORF_id_custom 列")
    n_main = len(main_df)
    log.info("主表行数：%d", n_main)

    # 优化2: 预转换为category类型以节省内存和加速比较
    main_df["ORF_id_custom"] = main_df["ORF_id_custom"].astype('category')

    # 优化3: 分块读取 orfs_3_ways，使用numpy array加速集合构建
    log.info("分块读取 orfs_3_ways：%s （chunksize=%d）", args.orfs_3_ways, args.chunksize)
    map_sw = {"price":"PRICE", "ribocode":"RiboCode", "riborf":"RibORF"}
    
    # 使用列表累积，最后一次性转set（比逐步update更快）
    price_list, rc_list, ro_list = [], [], []
    rows_price = rows_rc = rows_ro = 0
    
    try:
        pbar = None
        use_pbar = False
        if args.tqdm:
            try:
                from tqdm import tqdm
                pbar = tqdm(unit="rows", desc="reading orfs_3_ways")
                use_pbar = True
            except Exception:
                pass

        reader = pd.read_csv(
            args.orfs_3_ways,
            sep=args.sep3ways,
            header=None,
            usecols=[2,3],  # Software, ORF_id_custom
            names=["Software","ORF_id_custom"],
            dtype=str,
            engine=("c" if args.sep3ways == "\t" else "python"),  # 优化4: 优先使用C引擎
            compression="infer",
            quoting=3,
            chunksize=args.chunksize
        )
        total_rows = 0
        for i,chunk in enumerate(reader, 1):
            chunk_len = len(chunk)
            total_rows += chunk_len
            if use_pbar:
                pbar.update(chunk_len)
            
            # 优化5: 向量化操作，减少中间变量
            sw_lower = chunk["Software"].str.strip().str.lower()
            orf_ids = chunk["ORF_id_custom"].values  # 使用numpy array
            
            # 优化6: 使用布尔索引直接获取，避免多次map操作
            price_mask = sw_lower == "price"
            rc_mask = sw_lower == "ribocode"
            ro_mask = sw_lower == "riborf"
            
            # 计数（使用numpy sum更快）
            rows_price += int(price_mask.sum())
            rows_rc    += int(rc_mask.sum())
            rows_ro    += int(ro_mask.sum())
            
            # 累积到列表（比逐个update set快）
            if price_mask.any():
                price_list.append(orf_ids[price_mask])
            if rc_mask.any():
                rc_list.append(orf_ids[rc_mask])
            if ro_mask.any():
                ro_list.append(orf_ids[ro_mask])
            
            if i % 20 == 0:
                log.info("已处理 %d 行（块 %d）", total_rows, i)
        
        if use_pbar:
            pbar.close()
        log.info("orfs_3_ways 总行数：%d", total_rows)
    except pd.errors.EmptyDataError:
        log.warning("orfs_3_ways 为空；将全部标记为 FALSE")

    # 优化7: 一次性构建集合（比逐步update快）
    price_set = set(np.concatenate(price_list)) if price_list else set()
    rc_set = set(np.concatenate(rc_list)) if rc_list else set()
    ro_set = set(np.concatenate(ro_list)) if ro_list else set()

    log.info("PRICE: 行=%d, 唯一ORF=%d", rows_price, len(price_set))
    log.info("RiboCode: 行=%d, 唯一ORF=%d", rows_rc, len(rc_set))
    log.info("RibORF: 行=%d, 唯一ORF=%d", rows_ro, len(ro_set))

    # 优化8: 直接使用where生成字符串，避免map操作
    s = main_df["ORF_id_custom"]
    main_df["Is_PRICE"]    = np.where(s.isin(price_set), TRUE, FALSE)
    main_df["Is_RiboCode"] = np.where(s.isin(rc_set), TRUE, FALSE)
    main_df["Is_RibORF"]   = np.where(s.isin(ro_set), TRUE, FALSE)

    # 写出
    Path(args.out_with_flags).parent.mkdir(parents=True, exist_ok=True)
    main_df.to_csv(args.out_with_flags, sep="\t", index=False)
    log.info("已写出：%s（%d 行）", args.out_with_flags, len(main_df))

    # 优化9: 使用value_counts一次性统计
    def tf(df, col):
        vc = df[col].value_counts()
        return int(vc.get(TRUE,0)), int(vc.get(FALSE,0))
    t1,f1 = tf(main_df,"Is_PRICE")
    t2,f2 = tf(main_df,"Is_RiboCode")
    t3,f3 = tf(main_df,"Is_RibORF")

    print("\n[orfs_3_ways 软件计数]")
    print("Software\trows\tunique_orfs")
    print(f"PRICE\t{rows_price}\t{len(price_set)}")
    print(f"RiboCode\t{rows_rc}\t{len(rc_set)}")
    print(f"RibORF\t{rows_ro}\t{len(ro_set)}")

    print("\n[主表 Is_* 计数]")
    print("flag\tTRUE\tFALSE")
    print(f"Is_PRICE\t{t1}\t{f1}")
    print(f"Is_RiboCode\t{t2}\t{f2}")
    print(f"Is_RibORF\t{t3}\t{f3}")

    # 可选导出统计
    if args.out_stats:
        if args.out_stats.endswith((".xlsx",".xls")):
            with pd.ExcelWriter(args.out_stats, engine='openpyxl') as xw:
                pd.DataFrame(
                    {"Software":["PRICE","RiboCode","RibORF"],
                     "rows":[rows_price, rows_rc, rows_ro],
                     "unique_orfs":[len(price_set), len(rc_set), len(ro_set)]}
                ).to_excel(xw, sheet_name="software_stats", index=False)
                pd.DataFrame(
                    {"flag":["Is_PRICE","Is_RiboCode","Is_RibORF"],
                     "TRUE":[t1,t2,t3],
                     "FALSE":[f1,f2,f3]}
                ).to_excel(xw, sheet_name="is_flag_stats", index=False)
        else:
            base = Path(args.out_stats)
            pd.DataFrame(
                {"Software":["PRICE","RiboCode","RibORF"],
                 "rows":[rows_price, rows_rc, rows_ro],
                 "unique_orfs":[len(price_set), len(rc_set), len(ro_set)]}
            ).to_csv(str(base).replace(".tsv",".software.tsv"), sep="\t", index=False)
            pd.DataFrame(
                {"flag":["Is_PRICE","Is_RiboCode","Is_RibORF"],
                 "TRUE":[t1,t2,t3],
                 "FALSE":[f1,f2,f3]}
            ).to_csv(str(base).replace(".tsv",".is_flags.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()