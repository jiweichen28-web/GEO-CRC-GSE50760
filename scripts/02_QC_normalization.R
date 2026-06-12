# 02_QC_normalization.R — 质控、过滤、PCA
library(limma); library(ggplot2); library(ggrepel); library(RColorBrewer)
FONT_SIZE <- 8; PT <- FONT_SIZE/2.845; LW <- 0.5/2.1333
COL7 <- c("#00468A","#EC0000","#42B540","#0099B4","#925E9F","#FCAE91","#AC002A")

theme_pub_bw <- function() {
  theme_bw(base_size=FONT_SIZE) +
    theme(
      text          =element_text(size=FONT_SIZE, color="black", face="plain"),
      plot.title    =element_text(size=FONT_SIZE, color="black", face="plain", hjust=0.5),
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

expr_mat <- readRDS(file.path("data", "expr_matrix_log2.rds"))
pheno    <- readRDS(file.path("data", "sample_info.rds"))
stopifnot(all(colnames(expr_mat) == rownames(pheno)))
group <- factor(pheno$group, levels=c("Normal","Primary_CRC","Metastasis"))
# 三组配色：深蓝=Normal，正红=Primary_CRC，草绿=Metastasis
col3  <- c("Normal"=COL7[1], "Primary_CRC"=COL7[2], "Metastasis"=COL7[3])
cat("样本数:", ncol(expr_mat), " 基因数:", nrow(expr_mat), "\n")

# 过滤低表达基因
keep          <- rowMeans(expr_mat > 1) > (1/3)
expr_filtered <- expr_mat[keep, ]
cat("过滤后:", nrow(expr_filtered), "基因（过滤掉",
    nrow(expr_mat)-nrow(expr_filtered), "个）\n")
expr_norm <- expr_filtered

# PCA
pca_res <- prcomp(t(expr_norm), scale.=TRUE)
pca_df  <- data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2],
                      group=group, sample=rownames(pca_res$x))
var_exp <- round(summary(pca_res)$importance[2,1:2]*100, 1)

# shape=21：黑色细边框+填充色透明度，重叠点可分辨
p_pca <- ggplot(pca_df, aes(PC1, PC2, fill=group, label=sample)) +
  geom_point(shape=21, size=2, stroke=LW*2, color="black", alpha=0.7) +
  geom_text_repel(size=PT, max.overlaps=15,
                  color="black", face="plain") +
  scale_fill_manual(values=col3) +
  labs(title="PCA of GSE50760",
       x=paste0("PC1 (",var_exp[1],"%)"),
       y=paste0("PC2 (",var_exp[2],"%)"),
       fill="Group") +
  theme_pub_bw() + theme(legend.position="right")

# 箱线图
pbox <- ggplot(
  data.frame(
    value=as.vector(expr_mat),
    sample=rep(colnames(expr_mat), each=nrow(expr_mat)),
    group=rep(as.character(group), each=nrow(expr_mat))),
  aes(x=sample, y=value, fill=group)) +
  geom_boxplot(outlier.size=0.3, outlier.alpha=0.3, linewidth=LW) +
  scale_fill_manual(values=col3) +
  labs(title="Sample distribution (log2 FPKM)",
       x=NULL, y="log2(FPKM+1)", fill="Group") +
  theme_pub_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="right")

# 单图单 PDF 与 300 ppi TIF
fig_dir <- file.path("results", "figures", "")

save_pdf <- function(p, file) {
  cairo_pdf(file, width=10, height=5, pointsize=8,
            bg="transparent", onefile=FALSE)
  print(p); dev.off()
}

save_tif <- function(p, file) {
  tiff(file, width=10, height=5, units="in", res=300,
       pointsize=8, compression="lzw", bg="transparent", type="cairo")
  print(p); dev.off()
}

save_pdf(pbox, paste0(fig_dir, "02_QC_boxplot.pdf"))
save_tif(pbox, paste0(fig_dir, "02_QC_boxplot.tif"))
save_pdf(p_pca, paste0(fig_dir, "02_QC_PCA.pdf"))
save_tif(p_pca, paste0(fig_dir, "02_QC_PCA.tif"))

saveRDS(expr_norm, file.path("data", "expr_norm.rds"))
saveRDS(group,     file.path("data", "group_factor.rds"))
cat("✓ 完成\n")
