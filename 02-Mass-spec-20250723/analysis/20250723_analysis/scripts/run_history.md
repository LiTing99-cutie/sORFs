
# 20250723
## 该脚本处理≤150个氨基酸的人类UniProt蛋白质，并将其与Ensembl数据库进行匹配，以进行全面的注释分析。
nohup bash process_uniprot_small_proteins.v1.20250723.sh &> ../log/process_uniprot_small_proteins.20250723.log &
## 得到经典小肽的其他信息
cano.sep.tab.info.basic.20250723.ipynb
- 输出../processed/tab_info/cano.sep.tab.info.txt

# 20250724
## 得到质谱证据支持的小肽的被支持的PSM数量以及肽段数量
export.ms.info.20250724.ipynb
- 输出../processed/ms_res/ms_info.txt

## 结合蛋白质注释等级信息以及表达信息，查看注释蛋白被质谱鉴定的比例
cano.sep.add.ms.add.expr.iden.rate.20250724.ipynb
- 输出../processed/ms_res/cano.sep.tab.info.add.ms.txt

## 得到被质谱证据支持的非经典小肽的list的其他信息
noncano.sep.tab.info.basic.20250724.ipynb
- 输出../processed/tab_info/noncano.sep.ms.supported.tab.info.txt

## 非经典小肽类别绘图
2.3a.noncano.sep.type.plot.20250726.ipynb

## 整合正对照以及负对照
3.1.combine.sep.20250728.R.ipynb

## 在kozak序列上进行比较
3.2.sep.start.codon.kozak.20250726.R.ipynb

## 在肽段覆盖上进行比较
3.3.cano.nocano.pep.cov.20250728.R.ipynb