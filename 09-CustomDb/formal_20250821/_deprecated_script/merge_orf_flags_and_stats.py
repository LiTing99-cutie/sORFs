#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, logging
import pandas as pd
from pathlib import Path

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
    ap.add_argument("--chunksize", type=int, default=1_000_000, help="orfs_3_ways 按行分块大小")
    ap.add_argument("--in-with-custom", required=True, help="主表（含 ORF_id_custom）")
    ap.add_argument("--out-with-flags", required=True, help="输出表（新增 Is_*）")
    ap.add_argument("--out-stats", default=None, help="可选：导出统计（.tsv 或 .xlsx）")
    ap.add_argument("-v","--verbose", action="count", default=1, help="-v 普通日志，-vv 调试日志，默认 INFO")
    ap.add_argument("--tqdm", action="store_true", help="使用进度条（需安装 tqdm）")
    args = ap.parse_args()
    log = setup_logger(args.verbose)

    # 读取主表
    log.info("读取主表：%s", args.in_with_custom)
    main_df = pd.read_csv(args.in_with_custom, sep="\t", dtype=str, quoting=3, compression="infer")
    if "ORF_id_custom" not in main_df.columns:
        raise SystemExit("[ERR] 主表缺少 ORF_id_custom 列")
    n_main = len(main_df)
    log.info("主表行数：%d", n_main)

    # 分块读取 orfs_3_ways，构建三个集合 + 行数统计
    log.info("分块读取 orfs_3_ways：%s （chunksize=%d）", args.orfs_3_ways, args.chunksize)
    map_sw = {"price":"PRICE", "ribocode":"RiboCode", "riborf":"RibORF"}
    price_set, rc_set, ro_set = set(), set(), set()
    rows_price = rows_rc = rows_ro = 0
    try:
        pbar = None
        if args.tqdm:
            try:
                from tqdm import tqdm
                pbar = tqdm(unit="rows", desc="reading orfs_3_ways")
            except Exception:
                pbar = None

        reader = pd.read_csv(
            args.orfs_3_ways,
            sep=args.sep3ways,
            header=None,
            usecols=[2,3],  # Software, ORF_id_custom
            names=["Software","ORF_id_custom"],
            dtype=str,
            engine=("python" if args.sep3ways != "\t" else None),
            compression="infer",
            quoting=3,
            chunksize=args.chunksize
        )
        total_rows = 0
        for i,chunk in enumerate(reader, 1):
            total_rows += len(chunk)
            if pbar: pbar.update(len(chunk))
            # 规范软件名
            sw = chunk["Software"].str.strip().str.lower().map(map_sw)
            key = chunk["ORF_id_custom"].astype(str)
            # 计数
            rows_price += int((sw=="PRICE").sum())
            rows_rc    += int((sw=="RiboCode").sum())
            rows_ro    += int((sw=="RibORF").sum())
            # 集合
            price_set.update(key[sw=="PRICE"])
            rc_set.update(key[sw=="RiboCode"])
            ro_set.update(key[sw=="RibORF"])
            if i % 20 == 0:
                log.info("已处理 %d 行（块 %d）", total_rows, i)
        if pbar: pbar.close()
        log.info("orfs_3_ways 总行数：%d", total_rows)
    except pd.errors.EmptyDataError:
        log.warning("orfs_3_ways 为空；将全部标记为 FALSE")

    log.info("PRICE: 行=%d, 唯一ORF=%d", rows_price, len(price_set))
    log.info("RiboCode: 行=%d, 唯一ORF=%d", rows_rc, len(rc_set))
    log.info("RibORF: 行=%d, 唯一ORF=%d", rows_ro, len(ro_set))

    # 生成标记（向量化 isin）
    s = main_df["ORF_id_custom"]
    main_df["Is_PRICE"]    = s.isin(price_set).map({True: TRUE, False: FALSE})
    main_df["Is_RiboCode"] = s.isin(rc_set).map({True: TRUE, False: FALSE})
    main_df["Is_RibORF"]   = s.isin(ro_set).map({True: TRUE, False: FALSE})

    # 写出
    Path(args.out_with_flags).parent.mkdir(parents=True, exist_ok=True)
    main_df.to_csv(args.out_with_flags, sep="\t", index=False)
    log.info("已写出：%s（%d 行）", args.out_with_flags, len(main_df))

    # 统计并打印
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
            with pd.ExcelWriter(args.out_stats) as xw:
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
