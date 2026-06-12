# 01_download_data.R — GSE50760 数据整理（FPKM 文件版本）
library(GEOquery); library(tidyverse); library(stringr)

data_dir <- file.path("data", "")

# ── 1. 下载 Series Matrix（只用于触发缓存，实际数据来自 FPKM 文件）
getGEO("GSE50760", destdir=data_dir, getGPL=FALSE)

# ── 2. 读取 54 个 FPKM 文件并合并（cbind，不用 join）────
fpkm_files <- sort(list.files(data_dir, pattern="FPKM\\.txt\\.gz$", full.names=TRUE))
cat("找到 FPKM 文件:", length(fpkm_files), "个\n")

# 验证基因顺序一致（抽查首、次、末三个文件）
g1 <- read_tsv(fpkm_files[1],  col_types="cd", show_col_types=FALSE)[[1]]
g2 <- read_tsv(fpkm_files[2],  col_types="cd", show_col_types=FALSE)[[1]]
g3 <- read_tsv(fpkm_files[54], col_types="cd", show_col_types=FALSE)[[1]]
stopifnot("基因顺序不一致，不能直接 cbind" = identical(g1,g2) && identical(g1,g3))
cat("基因顺序验证通过 ✓\n")

# 逐文件读取第二列（FPKM值），cbind 拼成矩阵
vals     <- vapply(fpkm_files,
                   function(f) read_tsv(f, col_types="cd", show_col_types=FALSE)[[2]],
                   numeric(length(g1)))
rownames(vals) <- g1
colnames(vals) <- str_extract(basename(fpkm_files), "AMC_[0-9]+\\.[123]")
cat("表达矩阵维度:", dim(vals), "\n")  # 应为 23505 54

# ── 3. log2(FPKM+1) 转换 ────────────────────────────────
expr_log2 <- log2(vals + 1)
cat("log2 值域:", round(range(expr_log2, na.rm=TRUE), 2), "\n")

# ── 4. 从文件名构建样本信息表 ────────────────────────────
sample_ids <- colnames(expr_log2)
pheno <- data.frame(
  sample_id = sample_ids,
  group     = case_when(
    str_detect(sample_ids, "\\.1$") ~ "Primary_CRC",
    str_detect(sample_ids, "\\.2$") ~ "Normal",
    str_detect(sample_ids, "\\.3$") ~ "Metastasis"),
  row.names = sample_ids)
cat("分组统计:\n"); print(table(pheno$group))

# ── 5. 保存 ──────────────────────────────────────────────
saveRDS(expr_log2, paste0(data_dir, "expr_matrix_log2.rds"))
saveRDS(pheno,     paste0(data_dir, "sample_info.rds"))
write.csv(pheno,   paste0(data_dir, "sample_info.csv"))
cat("✓ 数据已保存\n")
