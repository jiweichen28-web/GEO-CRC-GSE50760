# GSE50760 — 结直肠癌 Bulk RNA-seq 分析

> 作者：陈季威 | 浙江大学海洋学院
> 标签：`bulk RNA-seq` `GEO` `差异表达` `富集` `GSEA` `WGCNA` `预后模型` `R`

从一个 GEO 数据集走完"差异 → 富集 → 共表达网络 → 预后建模"的完整流程。前半程（01–04）在 GSE50760 上做表达层面的分析；后半程（05–07）把挖到的基因带到 TCGA-COAD 上建生存预后模型。

## 数据来源

| 字段 | 内容 |
|------|------|
| GEO | [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760)，*Homo sapiens*，Illumina HiSeq 2000 |
| 设计 | 54 样本 = 18 例患者 × 3 组织（Normal / Primary CRC / Metastasis） |
| 外部队列 | TCGA-COAD（UCSC Xena GDC hub），434 病人带生存随访 |

## 分析流程

| 脚本 | 内容 |
|------|------|
| 01 | 数据下载与分组 |
| 02 | QC + 标准化 + PCA |
| 03 | 差异表达（limma）+ 火山图 + 热图 |
| 04 | GO/KEGG 富集（ORA）+ GSEA |
| 05 | WGCNA 共表达网络 → 模块-性状相关 → hub 基因 |
| 06 | 下载 TCGA-COAD 表达谱 + 生存数据 |
| 07 | 预后模型：单因素 Cox → LASSO-Cox → 风险评分 → KM / timeROC → 列线图 |

## 主要结果

**差异表达（limma，adj.P<0.05，|logFC|>1）**

| 对比 | 上调 | 下调 | 合计 |
|------|------|------|------|
| Primary CRC vs Normal | 316 | 521 | 837 |
| Metastasis vs Normal | 724 | 924 | 1648 |
| Metastasis vs CRC | 241 | 141 | 382 |

三组都是下调多于上调，提示 CRC 主要丢失正常肠上皮功能基因。最显著下调 OTOP2、BEST4（肠上皮标志），上调含 MMP1、WNT2（侵袭相关）。

**富集**：CRC 上调富集到细胞外基质重塑、PI3K-Akt；转移癌上调出现补体凝血、胆固醇代谢等肝转移微环境特征。GSEA 全局排序进一步抓到核糖体生物发生、DNA 复制等增殖信号。两组抑制端一致指向"刷状缘 + 离子转运"，即正常肠上皮分化功能的系统性丢失。

**WGCNA**：以肿瘤进展程度（Normal=0 / CRC=1 / Metastasis=2）为性状，9 个模块中 **blue 模块**与进展相关性最强（cor=−0.755，p=4.1e-11，随进展整体下调），含 463 个基因，从中筛出 hub 基因。

**预后模型**（建在 TCGA-COAD，434 病人 / 95 死亡 / 中位随访 22.1 月）：

- 候选 = blue 模块 hub 基因 ∩ TCGA 表达谱
- 单因素 Cox 留 43 个 → LASSO-Cox 选中 **20 基因签名** → 风险评分按中位数分高/低危
- KM 高/低危分离 **log-rank p=4.5e-8**
- timeROC **AUC 1/3/5 年 = 0.71 / 0.69 / 0.76**
- 多因素 Cox：风险评分、年龄、分期均独立预后；列线图 **C-index=0.765**

> 这套"GEO 挖基因 → TCGA 建预后模型"是肿瘤生信的常见范式。GSE50760 本身只有组织类型、无生存随访，所以预后必须借 TCGA 完成。

## 输出

- `results/figures/` — 每图 PDF（矢量，投稿用）+ TIF（300 ppi）
- `results/tables/` — 各步骤 CSV（差异基因、富集、模块成员、Cox/LASSO 系数、风险评分等）
- 各图含义与读图方法见 `notes-lessons/No.1_lessons_GEO_analysis`

## 环境

```r
# Bioconductor
GEOquery, limma, clusterProfiler, org.Hs.eg.db, enrichplot, WGCNA, AnnotationDbi
# CRAN
tidyverse, ggplot2, ggrepel, pheatmap, data.table
survival, survminer, glmnet, timeROC, rms
```

数据不纳入版本库（见 `.gitignore`）：GEO 数据按 accession 从 NCBI 下载，TCGA 数据由脚本 06 从 UCSC Xena 获取。

## 参考

Kim SK et al. A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients. *Mol Oncol* 2014;8(8):1653-66.
