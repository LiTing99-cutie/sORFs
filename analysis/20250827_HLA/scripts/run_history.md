nohup bash run.fragpipe.hla.1.20250827.sh &> ../log/hla.1.log &
# HLA-I型数据；修改数据库为100w+数据库，肽段大小为8-12aa；修改数据库为100w-数据库，肽段大小为8-12aa
nohup bash -c '
  bash run.fragpipe.hla.2025908.db.1.sh &> ../log/hla.db.1.log || true
  bash run.fragpipe.hla.2025908.db.2.sh &> ../log/hla.db.2.log
' >/dev/null 2>&1 &
# HLA-II型数据；修改数据库为100w+数据库，肽段大小为8-25aa；修改数据库为100w-数据库，肽段大小为8-25aa
nohup bash -c '
  bash run.fragpipe.hla.2.2025908.db.1.sh &> ../log/hla.2.db.1.log || true
  bash run.fragpipe.hla.2.2025908.db.2.sh &> ../log/hla.2.db.2.log
' >/dev/null 2>&1 &

nohup bash run.fragpipe.hla.2.2025908.db.2.sh &> ../log/hla.2.db.2.log &
nohup bash run.fragpipe.hla.2.2025908.db.1.sh &> ../log/hla.2.db.1.log &
nohup bash run.fragpipe.hla.2.2025908.db.3.sh &> ../log/hla.2.db.3.log &
nohup bash run.fragpipe.hla.2.2025908.db.4.sh &> ../log/hla.2.db.4.log &
