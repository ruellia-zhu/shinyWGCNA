# shinyWGCNA

shinyWGCNA is a Shiny-based graphical interface for WGCNA analysis on Linux and HPC servers. It wraps common WGCNA steps into an interactive workflow so users can upload expression data, clean/filter genes, choose a soft-threshold power, build gene co-expression modules, connect modules with traits, and export hub gene or Cytoscape-ready results.

## V3.3 更新内容

- 表达矩阵和 trait matrix 上传支持 `.csv`、tab 分隔 `.txt/.tsv`、`.xlsx`。
- 上传文件会自动判断文本分隔符或 xlsx 格式并读取。
- Data import 页面中 normalized expression 的显示统一为 `FPKM/TPM/CPM`。
- 应用版本号更新为 `V3.3`。
- 删除 `Before Use` 页面中的旧版 `NOTE` 说明区块。
- 新增 `AGENT.md`，用于记录后续维护和更新约定。

## 运行方式

在安装 R 和 Shiny 的服务器或本地环境中运行：

```r
shiny::runApp("shinyWGCNA.R")
```

也可以在 RStudio Server 中打开 `shinyWGCNA.R` 后点击 Run App。

脚本会按需安装缺失的 R 包。首次启动可能较慢，尤其是在需要安装 Bioconductor 或 GitHub 依赖时。建议在服务器管理员允许的 R library 路径中运行，并提前确认网络可以访问 CRAN、Bioconductor 和 GitHub。

## 输入文件格式

### Expression Matrix

在 `Data import and cleaning` 页面上传表达矩阵。

支持格式：

- comma-delimited `.csv`
- tab-delimited `.txt` 或 `.tsv`
- Excel `.xlsx`

文件要求：

- 第一列为 gene ID、probe ID 或 feature ID。
- 第二列开始为样本表达量。
- 表达量列必须为数值。
- 至少包含 1 个基因列和 2 个样本表达量列。
- `.xlsx` 默认读取第一个 worksheet。

示例：

```text
gene_id  sample_1  sample_2  sample_3
GeneA    12.4      10.8      15.1
GeneB    0.5       0.7       0.6
GeneC    23.1      19.5      20.4
```

### Trait Matrix

在 `Module-trait` 页面上传 trait matrix。

支持格式：

- comma-delimited `.csv`
- tab-delimited `.txt` 或 `.tsv`
- Excel `.xlsx`

文件要求：

- 第一列为 sample ID。
- 第二列开始为 trait 数值列。
- sample ID 必须能匹配 expression matrix 中的样本名。
- trait 列必须为数值。
- sample ID 不允许重复。
- `.xlsx` 默认读取第一个 worksheet。

示例：

```text
sample_id  treatment  time
sample_1   0          1
sample_2   1          1
sample_3   1          2
```

## 分析流程

### 1. Before Use

设置 WGCNA 可使用的服务器线程数。

- 默认线程数为 20。
- 点击 `Apply thread setting` 后，后续 WGCNA 步骤会使用当前设置。
- 如果服务器多人共用，建议根据机器负载适当降低线程数。

### 2. Data import and cleaning

上传 expression matrix，并设置数据类型和过滤参数。

- `Format` 可选择 `count` 或 `FPKM/TPM/CPM`。
- `Normalized method` 在 count 数据下使用 `VST`。
- `Normalized method` 在 `FPKM/TPM/CPM` 数据下可选择原始值或 `log10(FPKM/TPM/CPM)`。
- 第一轮过滤用于去除低表达基因。
- 第二轮过滤可按 MAD 或 Var 保留高变化基因。

完成导入后，页面会提示检测到的文件类型、分隔符或编码，并提供输入预览和过滤后数据预览。

### 3. Soft-threshold

根据清理后的表达矩阵选择 WGCNA soft-threshold power。

- 可以使用推荐 power。
- 也可以手动设置 power。
- 页面提供 scale-free topology 相关图表和表格下载。

### 4. Module-net

构建共表达网络并识别模块。

- 设置最小模块大小、合并阈值和 block size。
- 输出模块聚类图、模块数量表、TOMplot 和 gene-to-module 表格。

### 5. Module-trait

上传 trait matrix，并计算模块与表型/性状的相关性。

- expression matrix 的样本名必须能和 trait matrix 第一列 sample ID 对齐。
- trait matrix 后续列必须为数值。
- 页面输出 module-trait heatmap。

### 6. Hub Gene

根据模块成员关系和 trait 相关性筛选 hub genes。

- 可按 kME 和 GS 阈值筛选。
- 可下载 hub gene 表格。
- 可导出 Cytoscape 所需的 edge 和 node 文件。

## 输出文件

应用中的图表和表格下载按钮按分析流程分布在各页面中。常见输出包括：

- sample clustering tree
- soft-threshold plots and tables
- module dendrogram
- TOMplot
- gene-to-module table
- module-trait heatmap
- kME / GS hub gene table
- Cytoscape edge and node tables

## 常见注意事项

- `.xlsx` 文件默认读取第一个 worksheet；如果数据不在第一个 sheet，请先整理工作簿。
- 表达矩阵和 trait matrix 的样本名必须一致，否则 Module-trait 分析会失败。
- 表达矩阵中的空行、空列会在读取时自动移除，但表达量中的空值会阻止继续分析。
- 如果 count 数据过滤后剩余基因过少，请降低 expression cutoff 或调整 sample percentage。
- 如果上传文本文件解析错误，请确认文件确实使用逗号或 tab 分隔。
- 当前仓库环境可能没有 `R` 或 `Rscript`，最终启动测试应在可用 R/Shiny 的环境中完成。

## 维护说明

后续维护约定见 `AGENT.md`。
