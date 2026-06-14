# 06_tcga_download.R — 下载 TCGA-COAD 表达谱+生存数据，Ensembl→symbol，备预后建模
# 数据源 UCSC Xena GDC hub（公开、可复现）：STAR fpkm-uq + survival
# 产物 data/tcga_coad_expr.rds（symbol×肿瘤样本，log2）, data/tcga_coad_surv.rds
library(data.table); library(dplyr); library(org.Hs.eg.db); library(AnnotationDbi)
.a<-commandArgs(FALSE); .f<-sub("^--file=","",.a[grepl("^--file=",.a)])
if(length(.f)&&nzchar(.f[1])) setwd(dirname(dirname(normalizePath(.f[1],winslash="/"))))
dir.create("data/tcga", showWarnings=FALSE, recursive=TRUE)
base <- "https://gdc.xenahubs.net/download/"
# 文件名按 Xena 当前 STAR 流程（htseq 旧名已下线）
files <- c(
  expr = "TCGA-COAD.star_fpkm-uq.tsv.gz",
  surv = "TCGA-COAD.survival.tsv.gz",
  clin = "TCGA-COAD.clinical.tsv.gz")

# 下载（断点续传 + 重试 3 次；已存在且非空则跳过）
# 表达谱 128MB，必须放宽 timeout（默认60s 会半途中断）
options(timeout=1200)
fetch <- function(fn) {
  dest <- file.path("data/tcga", fn)
  if (file.exists(dest) && file.info(dest)$size > 1e3) {
    cat("已存在跳过:", fn, "\n"); return(dest) }
  url <- paste0(base, fn)
  for (i in 1:4) {
    if (file.exists(dest)) unlink(dest)   # 清掉上次的半截文件再重下
    ok <- tryCatch({
      download.file(url, dest, mode="wb", quiet=TRUE, method="libcurl")
      TRUE
    }, error=function(e){cat("第",i,"次失败:",conditionMessage(e),"\n");FALSE})
    if (ok && file.exists(dest) && file.info(dest)$size > 1e3) {
      cat("下载完成:", fn, "(", round(file.info(dest)$size/1e6,2), "MB)\n")
      return(dest) }
    Sys.sleep(5)
  }
  stop("下载失败: ", fn)
}
fp <- sapply(files, fetch)

# 1) 表达矩阵：行=Ensembl_ID(带版本号) 列=TCGA样本条码，值=log2(fpkm-uq+1)
expr <- fread(fp["expr"], data.table=FALSE)
ens <- sub("\\..*$", "", expr[[1]])   # 去 Ensembl 版本号
expr[[1]] <- NULL
expr <- as.matrix(expr)
# 去版本号后 PAR/Y 区基因会产生重复 Ensembl ID，按行均值保留最高者
ord  <- order(rowMeans(expr), decreasing=TRUE)
expr <- expr[ord, ]; ens <- ens[ord]
keep <- !duplicated(ens); expr <- expr[keep, ]; rownames(expr) <- ens[keep]
cat("TCGA 表达矩阵:", nrow(expr), "基因 x", ncol(expr), "样本\n")

# 2) Ensembl→symbol；一symbol多探针取平均表达最高者
map <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(expr),
         keytype="ENSEMBL", columns="SYMBOL")
map <- map[!is.na(map$SYMBOL) & !duplicated(map$ENSEMBL), ]
expr <- expr[map$ENSEMBL, ]; rownames(expr) <- map$SYMBOL
# 同名 symbol 去重：保留行均值最大的
ord  <- order(rowMeans(expr), decreasing=TRUE)
expr <- expr[ord, ]; expr <- expr[!duplicated(rownames(expr)), ]
cat("symbol 去重后:", nrow(expr), "基因\n")

# 3) 仅留原发肿瘤样本（条码第14-15位"01"=primary tumor），并裁条码到病人级
sampType <- substr(colnames(expr), 14, 15)
expr <- expr[, sampType=="01"]
colnames(expr) <- substr(colnames(expr), 1, 12)
expr <- expr[, !duplicated(colnames(expr))]
cat("原发肿瘤样本:", ncol(expr), "\n")

# 4) 生存数据：sample/OS.time/OS/_PATIENT；裁到病人级，月为单位
surv <- fread(fp["surv"], data.table=FALSE)
surv$patient <- substr(surv$sample, 1, 12)
surv <- surv[substr(surv$sample,14,15)=="01", ]
surv <- surv[!duplicated(surv$patient), ]
surv <- data.frame(patient=surv$patient,
                   OS=as.integer(surv$OS),
                   OS.time=as.numeric(surv$OS.time)/30.44,  # 天→月
                   row.names=surv$patient)
surv <- surv[!is.na(surv$OS) & !is.na(surv$OS.time) & surv$OS.time>0, ]

# 5) 临床协变量（年龄/性别/分期）供列线图
clin <- tryCatch({
  cl <- fread(fp["clin"], data.table=FALSE)
  cl$patient <- substr(cl$sample, 1, 12)
  pick <- function(nm) if(nm %in% names(cl)) cl[[nm]] else NA
  data.frame(patient=cl$patient,
    age=suppressWarnings(as.numeric(pick("age_at_index.demographic"))),
    gender=pick("gender.demographic"),
    stage=pick("ajcc_pathologic_stage.diagnoses"),
    stringsAsFactors=FALSE) %>% filter(!duplicated(patient))
}, error=function(e){cat("临床解析跳过:",conditionMessage(e),"\n");NULL})

# 6) 对齐三者交集
common <- Reduce(intersect, list(colnames(expr), rownames(surv),
            if(!is.null(clin)) clin$patient else rownames(surv)))
cat("表达∩生存", if(!is.null(clin))"∩临床" else "", "病人:", length(common), "\n")
expr <- expr[, common]; surv <- surv[common, ]
if (!is.null(clin)) { rownames(clin)<-clin$patient; clin<-clin[common,] }

saveRDS(expr, "data/tcga_coad_expr.rds")
saveRDS(surv, "data/tcga_coad_surv.rds")
if (!is.null(clin)) saveRDS(clin, "data/tcga_coad_clin.rds")
cat("事件数(死亡):", sum(surv$OS), " 中位随访:",
    round(median(surv$OS.time),1), "月\n")
cat("✓ 06 TCGA 下载完成\n")
