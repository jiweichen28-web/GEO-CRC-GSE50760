# GSE50760 — Colorectal Cancer Bulk RNA-seq Analysis

> 🧑‍💻 作者：陈季威 | 浙江大学海洋学院 博士在读
> 📅 开始日期：2026-06-11
> 🏷️ 标签：`bulk RNA-seq` `GEO` `差异表达` `富集分析` `GSEA` `结直肠癌` `R`

---

## 数据来源

| 字段 | 内容 |
|------|------|
| GEO Accession | [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760) |
| 物种 | *Homo sapiens* |
| 平台 | Illumina HiSeq 2000 |
| 样本数 | 54 samples，18 例患者，每例 3 个组织 |

| 分组 | 样本数 | 说明 |
|------|--------|------|
| Normal | 18 | 正常结肠组织 |
| Primary CRC | 18 | 原发结直肠癌 |
| Metastasis | 18 | 肝转移癌 |

---

## 分析内容

- [x] 01 数据下载与分组整理
- [x] 02 质控 + 标准化 + PCA
- [x] 03 差异表达分析（limma）+ 火山图 + 热图
- [x] 04 GO / KEGG 富集分析 + GSEA

---

## 主要结果

**差异基因数（limma，BH adj.P<0.05，|logFC|>1）**

| 对比 | 上调 | 下调 | 合计 |
|------|------|------|------|
| Primary CRC vs Normal | 316 | 521 | 837 |
| Metastasis vs Normal  | 724 | 924 | 1648 |
| Metastasis vs CRC     | 241 | 141 | 382 |

各组下调基因数均多于上调，提示 CRC 特征以正常肠上皮功能基因大规模丢失为主。
最显著下调基因：OTOP2、BEST4（正常肠上皮标志物）；上调基因含 MMP1、WNT2（肿瘤侵袭相关）。

**图件输出**

QC：

- `results/figures/02_QC_boxplot.pdf` / `.tif`
- `results/figures/02_QC_PCA.pdf` / `.tif`

差异表达：

- `results/figures/03_DEG_volcano_CRC_vs_Normal.pdf` / `.tif`
- `results/figures/03_DEG_volcano_Meta_vs_Normal.pdf` / `.tif`
- `results/figures/03_DEG_volcano_Meta_vs_CRC.pdf` / `.tif`
- `results/figures/03_DEG_heatmap_CRC_vs_Normal.pdf` / `.tif`
- `results/figures/03_DEG_heatmap_Meta_vs_Normal.pdf` / `.tif`
- `results/figures/03_DEG_heatmap_Meta_vs_CRC.pdf` / `.tif`

富集分析：

- `results/figures/04_GO_ORA_CRC_Up.pdf` / `.tif`
- `results/figures/04_GO_ORA_CRC_Down.pdf` / `.tif`
- `results/figures/04_GO_ORA_Meta_Up.pdf` / `.tif`
- `results/figures/04_GO_ORA_Meta_Down.pdf` / `.tif`
- `results/figures/04_KEGG_ORA_CRC_Up.pdf` / `.tif`
- `results/figures/04_KEGG_ORA_CRC_Down.pdf` / `.tif`
- `results/figures/04_KEGG_ORA_Meta_Up.pdf` / `.tif`
- `results/figures/04_KEGG_ORA_Meta_Down.pdf` / `.tif`
- `results/figures/04_GSEA_GO_CRC_vs_Normal.pdf` / `.tif`
- `results/figures/04_GSEA_GO_Meta_vs_Normal.pdf` / `.tif`

所有 PDF 均由 `cairo_pdf()` 输出；TIF 为 300 ppi。
- `results/tables/DEG_CRC_vs_Normal.csv` / `DEG_Meta_vs_Normal.csv` / `DEG_Meta_vs_CRC.csv`
- `results/tables/ORA_GO_CRC_Up.csv` / `ORA_GO_CRC_Down.csv` / `ORA_GO_Meta_Up.csv` / `ORA_GO_Meta_Down.csv`
- `results/tables/ORA_KEGG_CRC_Up.csv` / `ORA_KEGG_CRC_Down.csv` / `ORA_KEGG_Meta_Up.csv` / `ORA_KEGG_Meta_Down.csv`
- `results/tables/GSEA_GO_CRC.csv` / `GSEA_GO_Meta.csv`

**富集分析摘要**

ORA GO（ont="ALL"，BP/CC/MF 一并输出）— CRC 上调基因（316个）富集 292 条通路，主题为细胞外基质重塑（ECM organization）、胶原代谢、上皮形态发生；CRC 下调富集 242 条，集中在解毒/消化/刷状缘功能。转移癌上调富集 981 条（有机酸生物合成、体液免疫、补体激活、外源物代谢），下调 559 条（离子转运、白细胞迁移调控）。

ORA KEGG — CRC 上调富集 28 条（PI3K-Akt 信号、IL-17 信号、黏着斑/ECM-receptor 互作），下调 18 条；转移癌上调 32 条（补体与凝血级联、胆固醇代谢、细胞色素 P450 药物代谢），下调 15 条（矿物质吸收、胆汁分泌、胰腺分泌）——呈现明显的肝转移微环境特征。

GSEA GO（ont="ALL"）— CRC 显著富集 942 条（激活 574 / 抑制 368），激活以核糖体生物发生、DNA 复制、rRNA 加工为主（增殖特征），抑制以刷状缘、铜离子响应、离子转运为主。转移癌富集 790 条（激活 579 / 抑制 211），激活以急性期反应、血浆微粒、血小板 α 颗粒为主（肝脏/炎症-凝血特征），抑制同样落在刷状缘与氯离子转运。

---

## 环境依赖

```r
# Bioconductor
GEOquery, limma, clusterProfiler, org.Hs.eg.db, enrichplot

# CRAN
tidyverse, ggplot2, ggrepel, pheatmap, RColorBrewer
```

---

## 学习日志

| 日期 | 脚本 | 内容 | 遇到的问题 |
|------|------|------|------------|
| 2026-06-11 | 01 | 搭建项目框架，下载数据 | Series Matrix 为空（RNA-seq 补充文件下载）；full_join 内存爆炸（改用 cbind） |
| 2026-06-11 | 02 | QC + 标准化 + PCA | FPKM 不需要 quantile normalization |
| 2026-06-12 | 03 | limma 差异分析 + 火山图 + 热图 | topTable 行名为数字索引（fit\$genes=NULL，需回查 rownames）；LW 系数混用；pdf() 改 cairo_pdf()；补齐 3 组各自 Top50 热图（6页） |
| 2026-06-12 | 04 | GO/KEGG ORA + GSEA | enrichplot barplot() 新版 fill 参数已废弃，改用 dotplot()；bitr 5-10% 映射失败属正常；GSEA ties 警告无害；enrichGO/gseGO 用 ont="ALL" 一次覆盖 BP/CC/MF 并按 ONTOLOGY 分面；KEGG geneID 用 setReadable() 转回 symbol；补齐 CRC+Meta 两对比、GO/KEGG Up/Down、GSEA（10页） |

---

## 参考文献

Kim SK et al. A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients. *Mol Oncol* 2014;8(8):1653-66.
