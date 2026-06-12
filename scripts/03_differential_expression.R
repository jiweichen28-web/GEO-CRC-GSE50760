# 03_differential_expression.R — limma 差异分析、火山图、热图
library(limma); library(dplyr); library(ggplot2); library(ggrepel); library(pheatmap)
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

expr_norm <- readRDS(file.path("data", "expr_norm.rds"))
group     <- readRDS(file.path("data", "group_factor.rds"))

# limma 拟合
design <- model.matrix(~0 + group); colnames(design) <- levels(group)
fit    <- lmFit(expr_norm, design)
fit2   <- contrasts.fit(fit, makeContrasts(
  CRC_vs_Normal  = Primary_CRC - Normal,
  Meta_vs_Normal = Metastasis  - Normal,
  Meta_vs_CRC    = Metastasis  - Primary_CRC,
  levels=design))
fit2 <- eBayes(fit2)

# 提取 DEG（adj.P.Val<0.05，|logFC|>1）
# topTable 当 fit$genes=NULL 时行名为行索引数字，需回查 expr_norm 行名换回 symbol
get_DEG <- function(coef, lfc=1, p=0.05) {
  tt <- topTable(fit2, coef=coef, number=Inf, adjust.method="BH")
  tt$gene <- rownames(expr_norm)[as.integer(rownames(tt))]
  tt %>% mutate(direction=case_when(logFC>lfc & adj.P.Val<p ~ "Up",
                                    logFC< -lfc & adj.P.Val<p ~ "Down",
                                    TRUE ~ "NS"))
}
deg_CRC  <- get_DEG("CRC_vs_Normal")
deg_Meta <- get_DEG("Meta_vs_Normal")
deg_MC   <- get_DEG("Meta_vs_CRC")

# 打印各对比上调/下调基因数
for (nm in c("CRC_vs_Normal","Meta_vs_Normal","Meta_vs_CRC")) {
  d <- list(CRC_vs_Normal=deg_CRC, Meta_vs_Normal=deg_Meta, Meta_vs_CRC=deg_MC)[[nm]]
  cat(nm, "| Up:", sum(d$direction=="Up"),
      " Down:", sum(d$direction=="Down"), "\n")
}

# 保存 CSV 表格
out <- file.path("results", "tables", "")
write.csv(deg_CRC,  paste0(out,"DEG_CRC_vs_Normal.csv"),  row.names=FALSE)
write.csv(deg_Meta, paste0(out,"DEG_Meta_vs_Normal.csv"), row.names=FALSE)
write.csv(deg_MC,   paste0(out,"DEG_Meta_vs_CRC.csv"),    row.names=FALSE)

# 火山图函数（shape=21：fill 映射方向，黑色边框）
volcano <- function(deg, title) {
  top10 <- deg %>% filter(direction!="NS") %>% arrange(adj.P.Val) %>% head(10)
  ggplot(deg, aes(logFC, -log10(adj.P.Val), fill=direction)) +
    geom_point(shape=21, size=0.8, stroke=LW*2, color="black", alpha=0.7) +
    geom_text_repel(data=top10, aes(label=gene), size=PT,
                    color="black", fontface="plain", max.overlaps=20) +
    geom_vline(xintercept=c(-1,1), linetype="dashed",
               color="grey60", linewidth=LW) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               color="grey60", linewidth=LW) +
    # Up=正红, Down=深蓝, NS=浅灰（COL7体系）
    scale_fill_manual(values=c(Up=COL7[2], Down=COL7[1], NS="#BDC3C7")) +
    labs(title=title, x="log2FC", y="-log10(adj.P)", fill=NULL) +
    theme_pub_bw()
}

# 热图函数：每个对比取各自 Top50 DEG，仅取该对比涉及的两组样本
# scale="row" 行内 z-score；颜色 COL7 蓝-白-红；分组注释用 COL7 对应组色
draw_heatmap <- function(deg, grp_levels, ann_cols, title) {
  top50 <- deg %>% filter(direction!="NS") %>%
    arrange(adj.P.Val) %>% head(50) %>% pull(gene)
  sel   <- which(group %in% grp_levels)
  mat_h <- expr_norm[top50, sel]
  ann   <- data.frame(Group=group[sel]); rownames(ann) <- colnames(mat_h)
  pheatmap(mat_h, annotation_col=ann, annotation_colors=list(Group=ann_cols),
           scale="row", show_colnames=FALSE, fontsize=8, fontsize_row=6,
           color=colorRampPalette(c(COL7[1],"white",COL7[2]))(100),
           main=title, clustering_method="ward.D2")
}

# 单图单 PDF 与 300 ppi TIF
fig_dir <- file.path("results", "figures", "")

save_fig <- function(name, expr) {
  pdf_file <- paste0(fig_dir, name, ".pdf")
  tif_file <- paste0(fig_dir, name, ".tif")
  cairo_pdf(pdf_file, width=7, height=6, pointsize=8,
            bg="transparent", onefile=FALSE)
  eval(substitute(expr), parent.frame()); dev.off()
  tiff(tif_file, width=7, height=6, units="in", res=300,
       pointsize=8, compression="lzw", bg="transparent", type="cairo")
  eval(substitute(expr), parent.frame()); dev.off()
}

save_fig("03_DEG_volcano_CRC_vs_Normal",
         print(volcano(deg_CRC, "Primary CRC vs Normal")))
save_fig("03_DEG_volcano_Meta_vs_Normal",
         print(volcano(deg_Meta, "Metastasis vs Normal")))
save_fig("03_DEG_volcano_Meta_vs_CRC",
         print(volcano(deg_MC, "Metastasis vs CRC")))

save_fig("03_DEG_heatmap_CRC_vs_Normal",
         draw_heatmap(deg_CRC, c("Normal","Primary_CRC"),
                      c(Normal=COL7[1], Primary_CRC=COL7[2]),
                      "Top 50 DEGs: Primary CRC vs Normal"))
save_fig("03_DEG_heatmap_Meta_vs_Normal",
         draw_heatmap(deg_Meta, c("Normal","Metastasis"),
                      c(Normal=COL7[1], Metastasis=COL7[3]),
                      "Top 50 DEGs: Metastasis vs Normal"))
save_fig("03_DEG_heatmap_Meta_vs_CRC",
         draw_heatmap(deg_MC, c("Primary_CRC","Metastasis"),
                      c(Primary_CRC=COL7[2], Metastasis=COL7[3]),
                      "Top 50 DEGs: Metastasis vs Primary CRC"))

# 保存 DEG 结果供 04 使用
saveRDS(deg_CRC,  file.path("data", "deg_CRC_vs_Normal.rds"))
saveRDS(deg_Meta, file.path("data", "deg_Meta_vs_Normal.rds"))
saveRDS(deg_MC,   file.path("data", "deg_Meta_vs_CRC.rds"))
cat("✓ 完成\n")
