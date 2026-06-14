# 05_WGCNA.R — 加权基因共表达网络：软阈值→模块→模块-性状相关→hub基因→模块富集
# 接 02 的 expr_norm.rds / group_factor.rds，产出 hub 基因供 07 预后建模使用
library(WGCNA); library(ggplot2); library(dplyr); library(clusterProfiler)
library(org.Hs.eg.db)
options(stringsAsFactors=FALSE); disableWGCNAThreads()
# 自定位项目根：scripts/ 上溯一级，任意目录均可运行
.a<-commandArgs(FALSE); .f<-sub("^--file=","",.a[grepl("^--file=",.a)])
if(length(.f)&&nzchar(.f[1])) setwd(dirname(dirname(normalizePath(.f[1],winslash="/"))))
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
fig_dir <- "results/figures/"; tab_dir <- "results/tables/"

# ggplot 单图 PDF+TIF（与 02-04 一致）
save_fig <- function(p, name, width=7, height=6) {
  cairo_pdf(paste0(fig_dir,name,".pdf"), width=width, height=height,
            pointsize=8, bg="transparent", onefile=FALSE); print(p); dev.off()
  tiff(paste0(fig_dir,name,".tif"), width=width, height=height, units="in",
       res=300, pointsize=8, compression="lzw", bg="transparent", type="cairo")
  print(p); dev.off()
}
# base R 图（WGCNA 自带绘图为 base R）单图 PDF+TIF
# 绘图代码以函数传入（FUN），避免 substitute/eval.parent 不绘制到设备的惰性求值坑
save_base <- function(name, FUN, width=7, height=6) {
  cairo_pdf(paste0(fig_dir,name,".pdf"), width=width, height=height,
            pointsize=8, bg="white", onefile=FALSE)
  FUN(); dev.off()
  tiff(paste0(fig_dir,name,".tif"), width=width, height=height, units="in",
       res=300, pointsize=8, compression="lzw", bg="white", type="cairo")
  FUN(); dev.off()
}

# 数据：expr_norm 13362基因×54样本（log2 FPKM），group 三组因子
expr_norm <- readRDS(file.path("data","expr_norm.rds"))
group     <- readRDS(file.path("data","group_factor.rds"))
cat("输入:", nrow(expr_norm), "基因 x", ncol(expr_norm), "样本\n")

# 取 MAD 变异最大的 5000 基因（WGCNA 惯例，降噪+提速）
nTop <- 5000
mads <- apply(expr_norm, 1, mad)
sel  <- names(sort(mads, decreasing=TRUE))[1:min(nTop, sum(mads>0))]
datExpr <- t(expr_norm[sel, ])   # 行=样本 列=基因
cat("WGCNA 输入基因:", ncol(datExpr), " 样本:", nrow(datExpr), "\n")
gsg <- goodSamplesGenes(datExpr, verbose=0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# 1) 软阈值挑选（signed 网络，RsquaredCut=0.85）
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector=powers, networkType="signed",
                         RsquaredCut=0.85, verbose=0)
power <- sft$powerEstimate
if (is.na(power)) power <- 12   # 54样本 signed 网络经验值兜底
cat("选定软阈值 power =", power, "\n")
sfttab <- sft$fitIndices
sfttab$SFT.R.sq.signed <- -sign(sfttab$slope)*sfttab$SFT.R.sq
write.csv(sfttab, paste0(tab_dir,"05_WGCNA_softthreshold.csv"), row.names=FALSE)

# 软阈值诊断图：scale-free R² 与平均连通度（base R 双panel）
save_base("05_WGCNA_softpower", function() {
  par(mfrow=c(1,2), mar=c(4,4,2,1), cex=0.8, font.main=1, font.lab=1, font.axis=1)
  plot(sfttab$Power, sfttab$SFT.R.sq.signed, type="n", ylim=c(-0.4,1),
       xlab="Soft threshold (power)", ylab="Scale-free topology R^2",
       main="Scale independence")
  text(sfttab$Power, sfttab$SFT.R.sq.signed, labels=powers, col=COL7[2], cex=0.7)
  abline(h=0.85, col=COL7[1], lty=2)
  plot(sfttab$Power, sfttab$mean.k., type="n",
       xlab="Soft threshold (power)", ylab="Mean connectivity",
       main="Mean connectivity")
  text(sfttab$Power, sfttab$mean.k., labels=powers, col=COL7[2], cex=0.7)
}, width=8, height=4)

# 2) 一步法构建网络 + 模块识别（signed，所有基因一个 block）
cor <- WGCNA::cor   # 防止 cor 被其它包覆盖导致 blockwiseModules 报错
net_cache <- file.path("data","wgcna_net.rds")
if (file.exists(net_cache)) {
  net <- readRDS(net_cache); cat("复用已缓存网络:", net_cache, "\n")
} else {
  net <- blockwiseModules(datExpr, power=power, networkType="signed",
                        TOMType="signed", minModuleSize=30, deepSplit=2,
                        mergeCutHeight=0.25, numericLabels=TRUE,
                        maxBlockSize=6000, saveTOMs=FALSE,
                        pamRespectsDendro=FALSE, verbose=0, seed=42)
  saveRDS(net, net_cache)
}
cor <- stats::cor
moduleColors <- labels2colors(net$colors)
nMod <- length(unique(moduleColors[moduleColors!="grey"]))
cat("识别模块数（不含grey）:", nMod, "\n")
print(table(moduleColors))

# 模块树状图
save_base("05_WGCNA_dendrogram", function() {
  plotDendroAndColors(net$dendrograms[[1]],
    moduleColors[net$blockGenes[[1]]], "Module",
    dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
    main="Gene dendrogram and module colors", cex.main=0.9,
    cex.axis=0.8, cex.lab=0.8)
}, width=8, height=5)

# 3) 模块-性状相关：进展(0/1/2) + 三组各自二元
prog <- as.integer(group)-1
datTraits <- data.frame(Progression=prog,
  Normal=as.integer(group=="Normal"),
  Primary_CRC=as.integer(group=="Primary_CRC"),
  Metastasis=as.integer(group=="Metastasis"))
rownames(datTraits) <- rownames(datExpr)

MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
moduleTraitCor <- cor(MEs, datTraits, use="p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
write.csv(data.frame(module=rownames(moduleTraitCor), moduleTraitCor,
          check.names=FALSE), paste0(tab_dir,"05_WGCNA_module_trait_cor.csv"),
          row.names=FALSE)

# 模块-性状相关热图（每格：相关系数\n(p值)）
textMat <- paste0(signif(moduleTraitCor,2),"\n(",signif(moduleTraitP,1),")")
dim(textMat) <- dim(moduleTraitCor)
save_base("05_WGCNA_module_trait_heatmap", function() {
  par(mar=c(5,8,2,1), cex=0.7, font.main=1, font.lab=1, font.axis=1)
  labeledHeatmap(Matrix=moduleTraitCor, xLabels=colnames(datTraits),
    yLabels=rownames(moduleTraitCor), ySymbols=rownames(moduleTraitCor),
    colorLabels=FALSE,
    colors=blueWhiteRed(50), textMatrix=textMat, setStdMargins=FALSE,
    cex.text=0.6, zlim=c(-1,1), main="Module-trait relationships")
}, width=6, height=max(4, 0.5*ncol(MEs)+2))

# 4) 目标模块 = 与 Progression 绝对相关最高的非 grey 模块
mt <- data.frame(module=rownames(moduleTraitCor),
                 cor=moduleTraitCor[,"Progression"],
                 p=moduleTraitP[,"Progression"])
mt <- mt[mt$module!="MEgrey", ]
mt <- mt[order(-abs(mt$cor)), ]
targetME <- mt$module[1]; targetColor <- sub("^ME","",targetME)
cat("\n目标模块:", targetColor, " cor(Progression)=",
    round(mt$cor[1],3), " p=", signif(mt$p[1],2), "\n")

# 5) hub 基因：模块成员度 MM 与基因显著性 GS
GS <- as.numeric(cor(datExpr, prog, use="p"))
names(GS) <- colnames(datExpr)
MM <- cor(datExpr, MEs, use="p")
modGenes <- colnames(datExpr)[moduleColors==targetColor]
mmCol <- MM[, targetME]
hubTab <- data.frame(gene=modGenes, module=targetColor,
            MM=mmCol[modGenes], GS=GS[modGenes]) %>%
  mutate(absMM=abs(MM), absGS=abs(GS)) %>% arrange(desc(absMM))
write.csv(hubTab, paste0(tab_dir,"05_WGCNA_",targetColor,"_members.csv"),
          row.names=FALSE)
# hub 判定：|MM|>0.8 且 |GS|>0.5；不足则放宽取 |MM| top30
hub <- hubTab %>% filter(absMM>0.8, absGS>0.5) %>% pull(gene)
if (length(hub) < 10) hub <- head(hubTab$gene, 30)
cat("hub 基因数:", length(hub), "\n")
saveRDS(hub, file.path("data","wgcna_hub_genes.rds"))
write.csv(hubTab %>% filter(gene %in% hub),
          paste0(tab_dir,"05_WGCNA_hub_genes.csv"), row.names=FALSE)

# MM vs GS 散点（目标模块），hub 高亮+标注
mmgs <- hubTab %>% mutate(isHub=gene %in% hub)
p_mmgs <- ggplot(mmgs, aes(absMM, absGS)) +
  geom_point(shape=21, size=1.6, stroke=LW*2, color="black",
             aes(fill=isHub), alpha=0.7) +
  ggrepel::geom_text_repel(data=head(mmgs[mmgs$isHub,],15),
    aes(label=gene), size=PT, color="black", fontface="plain",
    max.overlaps=20) +
  geom_vline(xintercept=0.8, linetype="dashed", color="grey60", linewidth=LW) +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey60", linewidth=LW) +
  scale_fill_manual(values=c(`FALSE`="#BDC3C7",`TRUE`=COL7[2]), guide="none") +
  labs(title=paste0("Module ",targetColor,": MM vs GS"),
       x="Module membership (|MM|)", y="Gene significance (|GS| to progression)") +
  theme_pub_bw()
save_fig(p_mmgs, "05_WGCNA_MM_vs_GS", width=6, height=5)

# 模块基因数条形图
modSize <- as.data.frame(table(moduleColors)) %>%
  filter(moduleColors!="grey") %>% arrange(desc(Freq))
modSize$moduleColors <- factor(modSize$moduleColors, levels=modSize$moduleColors)
p_size <- ggplot(modSize, aes(moduleColors, Freq)) +
  geom_col(fill=COL7[1], width=0.7) +
  geom_text(aes(label=Freq), vjust=-0.3, size=PT, color="black") +
  labs(title="Module sizes", x="Module", y="Number of genes") +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  theme_pub_bw() + theme(axis.text.x=element_text(angle=45, hjust=1,
    size=FONT_SIZE, color="black", face="plain"))
save_fig(p_size, "05_WGCNA_module_sizes", width=6, height=4)

# 6) 目标模块 GO 富集（try 包裹，失败不影响主产物）
try({
  eid <- bitr(modGenes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  ego <- enrichGO(eid, org.Hs.eg.db, ont="BP", pAdjustMethod="BH",
                  pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego))>0) {
    write.csv(as.data.frame(ego),
      paste0(tab_dir,"05_WGCNA_",targetColor,"_GO_BP.csv"), row.names=FALSE)
    p_go <- dotplot(ego, showCategory=15, color="p.adjust", label_format=45) +
      scale_color_gradient(low=COL7[2], high=COL7[6]) +
      labs(title=paste0("GO BP: module ",targetColor)) + theme_pub_bw() +
      theme(axis.text.y=element_text(size=FONT_SIZE*0.75, color="black",
        face="plain"), legend.position="right")
    save_fig(p_go, paste0("05_WGCNA_",targetColor,"_GO_BP"), width=8, height=7)
    cat("目标模块 GO BP 显著条目:", sum(ego@result$p.adjust<0.05), "\n")
  }
}, silent=FALSE)

cat("✓ 05 WGCNA 完成；hub 基因已存 data/wgcna_hub_genes.rds\n")
