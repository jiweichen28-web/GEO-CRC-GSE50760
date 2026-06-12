# 04_enrichment_analysis.R — ORA (GO/KEGG) + GSEA，覆盖 CRC 和 Meta 两个对比
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
library(ggplot2); library(dplyr)
FONT_SIZE <- 8; PT <- FONT_SIZE/2.845; LW <- 0.5/2.1333
COL7 <- c("#00468A","#EC0000","#42B540","#0099B4","#925E9F","#FCAE91","#AC002A")

theme_pub_bw <- function() {
  theme_bw(base_size=FONT_SIZE) +
    theme(
      text          =element_text(size=FONT_SIZE, color="black", face="plain"),
      plot.title    =element_text(size=FONT_SIZE, color="black", face="plain", hjust=0.5),
      plot.subtitle =element_text(size=FONT_SIZE, color="black", face="plain"),
      plot.caption  =element_text(size=FONT_SIZE, color="black", face="plain"),
      axis.title    =element_text(size=FONT_SIZE, color="black", face="plain"),
      axis.text     =element_text(size=FONT_SIZE, color="black", face="plain"),
      legend.title  =element_text(size=FONT_SIZE, color="black", face="plain"),
      legend.text   =element_text(size=FONT_SIZE, color="black", face="plain"),
      strip.text    =element_text(size=FONT_SIZE, color="black", face="plain"),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
      panel.border  =element_rect(color="black", fill=NA, linewidth=LW),
      axis.line     =element_blank(),
      axis.ticks    =element_line(color="black", linewidth=LW),
      panel.background=element_rect(fill="transparent", color=NA),
      plot.background =element_rect(fill="transparent", color=NA),
      legend.background=element_rect(fill="transparent", color=NA),
      legend.key    =element_rect(fill="transparent", color=NA)
    )
}

# 读取 DEG 结果
deg_CRC  <- readRDS(file.path("data", "deg_CRC_vs_Normal.rds"))
deg_Meta <- readRDS(file.path("data", "deg_Meta_vs_Normal.rds"))
out <- file.path("results", "tables", "")

# Symbol → Entrez ID 转换（bitr 会丢弃无法映射的基因，5-10% 失败属正常）
to_entrez <- function(symbols) {
  bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db,
       drop=TRUE)$ENTREZID
}

# 提取各方向基因列表
up_CRC   <- deg_CRC  %>% filter(direction=="Up")   %>% pull(gene)
down_CRC <- deg_CRC  %>% filter(direction=="Down") %>% pull(gene)
up_meta  <- deg_Meta %>% filter(direction=="Up")   %>% pull(gene)
down_meta<- deg_Meta %>% filter(direction=="Down") %>% pull(gene)
cat("CRC  Up:", length(up_CRC), " Down:", length(down_CRC), "\n")
cat("Meta Up:", length(up_meta), " Down:", length(down_meta), "\n")

# ORA — enrichGO(ont="ALL" 一次覆盖 BP/CC/MF) + enrichKEGG
# ont="ALL" 返回结果含 ONTOLOGY 列，dotplot 可按本体分面
run_ora <- function(genes, label) {
  eid <- to_entrez(genes)
  go  <- enrichGO(eid, OrgDb=org.Hs.eg.db, ont="ALL",
                  pAdjustMethod="BH", pvalueCutoff=0.05,
                  qvalueCutoff=0.2, readable=TRUE)
  kg  <- enrichKEGG(eid, organism="hsa",
                    pAdjustMethod="BH", pvalueCutoff=0.05)
  # KEGG 结果 geneID 默认是 entrez，转为可读 symbol
  if (!is.null(kg) && nrow(kg@result) > 0)
    kg <- setReadable(kg, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
  go_n <- if (is.null(go)) 0 else sum(go@result$p.adjust < 0.05)
  kg_n <- if (is.null(kg)) 0 else sum(kg@result$p.adjust < 0.05)
  cat(label, "| GO terms:", go_n, " KEGG pathways:", kg_n, "\n")
  list(go=go, kegg=kg)
}
ora_CRC_up   <- run_ora(up_CRC,    "CRC  Up  ")
ora_CRC_dn   <- run_ora(down_CRC,  "CRC  Down")
ora_Meta_up  <- run_ora(up_meta,   "Meta Up  ")
ora_Meta_dn  <- run_ora(down_meta, "Meta Down")

# GSEA — 全基因按 logFC 排序，gseGO(ont="ALL") 覆盖 BP/CC/MF
make_gsea_input <- function(deg) {
  m <- bitr(deg$gene, fromType="SYMBOL", toType="ENTREZID",
            OrgDb=org.Hs.eg.db, drop=TRUE)
  d <- deg %>% filter(gene %in% m$SYMBOL) %>%
    left_join(m, by=c("gene"="SYMBOL"))
  sort(setNames(d$logFC, d$ENTREZID), decreasing=TRUE)
}
gl_CRC  <- make_gsea_input(deg_CRC)
gl_Meta <- make_gsea_input(deg_Meta)
set.seed(42)
gsea_CRC  <- gseGO(gl_CRC,  OrgDb=org.Hs.eg.db, ont="ALL",
                   minGSSize=10, maxGSSize=500,
                   pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE)
gsea_Meta <- gseGO(gl_Meta, OrgDb=org.Hs.eg.db, ont="ALL",
                   minGSSize=10, maxGSSize=500,
                   pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE)
cat("GSEA CRC  sig terms:", sum(gsea_CRC@result$p.adjust<0.05), "\n")
cat("GSEA Meta sig terms:", sum(gsea_Meta@result$p.adjust<0.05), "\n")

# 保存结果表（GO ORA 4 个方向 + KEGG 4 个方向 + GSEA 2 个对比）
save_csv <- function(obj, fn) {
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0)
    write.csv(as.data.frame(obj), paste0(out, fn), row.names=FALSE)
}
save_csv(ora_CRC_up$go,    "ORA_GO_CRC_Up.csv")
save_csv(ora_CRC_dn$go,    "ORA_GO_CRC_Down.csv")
save_csv(ora_Meta_up$go,   "ORA_GO_Meta_Up.csv")
save_csv(ora_Meta_dn$go,   "ORA_GO_Meta_Down.csv")
save_csv(ora_CRC_up$kegg,  "ORA_KEGG_CRC_Up.csv")
save_csv(ora_CRC_dn$kegg,  "ORA_KEGG_CRC_Down.csv")
save_csv(ora_Meta_up$kegg, "ORA_KEGG_Meta_Up.csv")
save_csv(ora_Meta_dn$kegg, "ORA_KEGG_Meta_Down.csv")
save_csv(gsea_CRC,         "GSEA_GO_CRC.csv")
save_csv(gsea_Meta,        "GSEA_GO_Meta.csv")

# ORA dotplot：ont="ALL" 用 split="ONTOLOGY" 保证 BP/CC/MF 各取 top-N，按本体分面
# color 映射 p.adjust，up 用红系、down 用蓝系（COL7 体系）
ora_dot <- function(eo, title, low_col) {
  if (is.null(eo) || sum(eo@result$p.adjust<0.05)==0)
    return(ggplot() + labs(title=paste0(title," (无显著通路)")) + theme_pub_bw())
  dotplot(eo, showCategory=5, split="ONTOLOGY", color="p.adjust", label_format=45) +
    facet_grid(ONTOLOGY~., scales="free") +
    scale_color_gradient(low=low_col, high=COL7[6]) +
    labs(title=title) + theme_pub_bw() +
    # 通路名较长，y轴标签缩至0.75倍防溢出
    theme(axis.text.y=element_text(size=FONT_SIZE*0.75, color="black", face="plain"),
          legend.position="right")
}

# KEGG dotplot（barplot 在新版 enrichplot 中 fill 参数已废弃，统一用 dotplot）
kegg_dot <- function(eo, title, low_col) {
  if (is.null(eo) || sum(eo@result$p.adjust<0.05)==0)
    return(ggplot() + labs(title=paste0(title," (无显著通路)")) + theme_pub_bw())
  dotplot(eo, showCategory=15, color="p.adjust", label_format=45) +
    scale_color_gradient(low=low_col, high=COL7[6]) +
    labs(title=title) + theme_pub_bw() +
    theme(axis.text.y=element_text(size=FONT_SIZE*0.75, color="black", face="plain"),
          legend.position="right")
}

# GSEA dotplot：split="ONTOLOGY" 保证 BP/CC/MF 各取 top-N，NES 着色（红激活/蓝抑制）
gsea_dot <- function(eo, title) {
  if (is.null(eo) || sum(eo@result$p.adjust<0.05)==0)
    return(ggplot() + labs(title=paste0(title," (无显著通路)")) + theme_pub_bw())
  dotplot(eo, showCategory=5, split="ONTOLOGY", color="NES", label_format=45) +
    facet_grid(ONTOLOGY~., scales="free") +
    scale_color_gradient2(low=COL7[1], mid="white", high=COL7[2], midpoint=0) +
    labs(title=title) + theme_pub_bw() +
    theme(axis.text.y=element_text(size=FONT_SIZE*0.75, color="black", face="plain"),
          legend.position="right")
}

# 单图单 PDF 与 300 ppi TIF
fig_dir <- file.path("results", "figures", "")

save_fig <- function(p, name) {
  pdf_file <- paste0(fig_dir, name, ".pdf")
  tif_file <- paste0(fig_dir, name, ".tif")
  cairo_pdf(pdf_file, width=9, height=8, pointsize=8,
            bg="transparent", onefile=FALSE)
  print(p); dev.off()
  tiff(tif_file, width=9, height=8, units="in", res=300,
       pointsize=8, compression="lzw", bg="transparent", type="cairo")
  print(p); dev.off()
}

plots <- list(
  "04_GO_ORA_CRC_Up"=
    ora_dot(ora_CRC_up$go, "GO ORA: CRC Up-regulated", COL7[2]),
  "04_GO_ORA_CRC_Down"=
    ora_dot(ora_CRC_dn$go, "GO ORA: CRC Down-regulated", COL7[1]),
  "04_GO_ORA_Meta_Up"=
    ora_dot(ora_Meta_up$go, "GO ORA: Metastasis Up-regulated", COL7[2]),
  "04_GO_ORA_Meta_Down"=
    ora_dot(ora_Meta_dn$go, "GO ORA: Metastasis Down-regulated", COL7[1]),
  "04_KEGG_ORA_CRC_Up"=
    kegg_dot(ora_CRC_up$kegg, "KEGG ORA: CRC Up-regulated", COL7[2]),
  "04_KEGG_ORA_CRC_Down"=
    kegg_dot(ora_CRC_dn$kegg, "KEGG ORA: CRC Down-regulated", COL7[1]),
  "04_KEGG_ORA_Meta_Up"=
    kegg_dot(ora_Meta_up$kegg, "KEGG ORA: Metastasis Up-regulated", COL7[2]),
  "04_KEGG_ORA_Meta_Down"=
    kegg_dot(ora_Meta_dn$kegg, "KEGG ORA: Metastasis Down-regulated",
             COL7[1]),
  "04_GSEA_GO_CRC_vs_Normal"=
    gsea_dot(gsea_CRC, "GSEA GO: CRC vs Normal"),
  "04_GSEA_GO_Meta_vs_Normal"=
    gsea_dot(gsea_Meta, "GSEA GO: Metastasis vs Normal")
)

for (nm in names(plots)) save_fig(plots[[nm]], nm)

cat("✓ 完成\n")
