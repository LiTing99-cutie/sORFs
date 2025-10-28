
## 数据集
对于单个样本（使用C8柱富集小肽，使用trypsin进行酶切）进行LC-MS/MS实验

## 分析目的
以数据库搜库的结果作为golden standard,benchmark不同的de novo测序方法
- casanovo
- primenovo
- novor
- punifind
设置母离子容差为20ppm

数据库搜库
分别使用了msfragger的closed search和open search进行搜库
以及pfind的closed search以及open search进行搜库
比较msfragger和pfind的结果

得到最终谱图肽段对应的结果文件:
- 谱图ID
- 肽段
- 来源的鉴定方法