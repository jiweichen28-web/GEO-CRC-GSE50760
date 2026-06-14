# 07_prognosis.R — TCGA-COAD 预后模型：单因素Cox→LASSO-Cox→风险评分→KM/timeROC→列线图
# 输入：05 的 hub 基因 + 06 的 TCGA 表达谱/生存/临床。基因来自 GEO(GSE50760)，模型建在 TCGA
library(survival); library(survminer); library(glmnet); library(timeROC)
library(rms); library(dplyr); library(ggplot2)
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
save_fig <- function(p, name, width=7, height=6) {
  cairo_pdf(paste0(fig_dir,name,".pdf"), width=width, height=height,
            pointsize=8, bg="transparent", onefile=FALSE); print(p); dev.off()
  tiff(paste0(fig_dir,name,".tif"), width=width, height=height, units="in",
       res=300, pointsize=8, compression="lzw", bg="transparent", type="cairo")
  print(p); dev.off()
}
save_base <- function(name, FUN, width=7, height=6) {
  cairo_pdf(paste0(fig_dir,name,".pdf"), width=width, height=height,
            pointsize=8, bg="white", onefile=FALSE)
  FUN(); dev.off()
  tiff(paste0(fig_dir,name,".tif"), width=width, height=height, units="in",
       res=300, pointsize=8, compression="lzw", bg="white", type="cairo")
  FUN(); dev.off()
}

# 数据加载
expr <- readRDS("data/tcga_coad_expr.rds")    # symbol × 病人，log2
surv <- readRDS("data/tcga_coad_surv.rds")     # OS, OS.time(月)
hub  <- readRDS("data/wgcna_hub_genes.rds")    # GEO 来源候选基因
clin <- tryCatch(readRDS("data/tcga_coad_clin.rds"), error=function(e) NULL)
cand <- intersect(hub, rownames(expr))
cat("候选基因(hub∩TCGA表达):", length(cand), "/", length(hub), "\n")
cat("样本数:", ncol(expr), " 事件:", sum(surv$OS), "\n")

# 构建建模数据框：行=病人，列=OS/OS.time + 候选基因表达
X <- t(expr[cand, rownames(surv)])
dat <- data.frame(OS.time=surv$OS.time, OS=surv$OS, X, check.names=FALSE)

# 1) 单因素 Cox 初筛（p<0.05）
uni <- lapply(cand, function(g){
  f <- tryCatch(coxph(Surv(OS.time,OS) ~ dat[[g]], data=dat),
                error=function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- summary(f)
  data.frame(gene=g, HR=s$coef[,"exp(coef)"], z=s$coef[,"z"],
             p=s$coef[,"Pr(>|z|)"],
             lower=s$conf.int[,"lower .95"], upper=s$conf.int[,"upper .95"])
})
uni <- bind_rows(uni) %>% arrange(p)
write.csv(uni, paste0(tab_dir,"07_univariate_cox.csv"), row.names=FALSE)
sig <- uni %>% filter(p<0.05) %>% pull(gene)
cat("单因素 Cox 显著(p<0.05):", length(sig), "个\n")
# 不足 5 个则放宽到 p<0.1，再不足取 z 值 top10 保证 LASSO 有输入
if (length(sig) < 5) sig <- uni %>% filter(p<0.1) %>% pull(gene)
if (length(sig) < 5) sig <- head(uni$gene, 10)
cat("进入 LASSO 的基因:", length(sig), "个\n")

# 单因素 Cox 森林图（top 显著基因）
fdf <- uni %>% filter(gene %in% sig) %>% arrange(HR) %>%
  mutate(gene=factor(gene, levels=gene),
         risk=ifelse(HR>1,"Risk","Protective"))
p_forest <- ggplot(fdf, aes(HR, gene, color=risk)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60", linewidth=LW) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.25, linewidth=LW) +
  geom_point(size=1.6) +
  scale_color_manual(values=c(Risk=COL7[2], Protective=COL7[1])) +
  scale_x_log10() +
  labs(title="Univariate Cox (TCGA-COAD)", x="Hazard ratio (log scale)",
       y=NULL, color=NULL) + theme_pub_bw()
save_fig(p_forest, "07_univariate_forest", width=6, height=max(3,0.25*nrow(fdf)+1.5))

# 2) LASSO-Cox（10折CV，lambda.min）
set.seed(42)
xmat <- as.matrix(dat[, sig, drop=FALSE])
ysur <- Surv(dat$OS.time, dat$OS)
cvf  <- cv.glmnet(xmat, ysur, family="cox", alpha=1, nfolds=10)
fit  <- glmnet(xmat, ysur, family="cox", alpha=1, lambda=cvf$lambda.min)
coef_lasso <- coef(fit); sel <- rownames(coef_lasso)[as.numeric(coef_lasso)!=0]
coefs <- setNames(as.numeric(coef_lasso)[as.numeric(coef_lasso)!=0], sel)
cat("LASSO 选中基因:", length(sel), "->", paste(sel,collapse=", "), "\n")
write.csv(data.frame(gene=sel, coef=coefs),
          paste0(tab_dir,"07_lasso_coefficients.csv"), row.names=FALSE)

# LASSO CV 曲线
save_base("07_lasso_cv", function() {
  par(mar=c(4,4,3,1), cex=0.8, font.main=1, font.lab=1, font.axis=1)
  plot(cvf); title("LASSO-Cox CV", line=2.5, cex.main=0.9)
}, width=6, height=5)
# LASSO 系数路径
save_base("07_lasso_path", function() {
  par(mar=c(4,4,3,1), cex=0.8, font.main=1, font.lab=1, font.axis=1)
  plot(glmnet(xmat, ysur, family="cox", alpha=1), xvar="lambda")
  abline(v=log(cvf$lambda.min), lty=2, col=COL7[2])
}, width=6, height=5)

# 3) 风险评分 = Σ coef×表达；按中位数分高/低危
risk_score <- as.numeric(xmat[, sel, drop=FALSE] %*% coefs)
dat$riskScore <- risk_score
dat$riskGroup <- factor(ifelse(risk_score > median(risk_score),
                               "High","Low"), levels=c("Low","High"))
write.csv(data.frame(patient=rownames(dat), riskScore=dat$riskScore,
          riskGroup=dat$riskGroup, OS=dat$OS, OS.time=dat$OS.time),
          paste0(tab_dir,"07_risk_score.csv"), row.names=FALSE)

# 4) KM 生存曲线（高 vs 低危）
sfit <- survfit(Surv(OS.time,OS) ~ riskGroup, data=dat)
km <- ggsurvplot(sfit, data=dat, pval=TRUE, risk.table=TRUE,
        palette=c(COL7[1],COL7[2]), legend.labs=c("Low risk","High risk"),
        legend.title="", xlab="Time (months)", ylab="Overall survival",
        ggtheme=theme_pub_bw(), risk.table.height=0.28,
        font.main=8, font.x=8, font.y=8, font.tickslab=8, font.legend=8,
        pval.size=3)
cairo_pdf(paste0(fig_dir,"07_KM_riskgroup.pdf"), width=6, height=6,
          pointsize=8, bg="white", onefile=FALSE); print(km); dev.off()
tiff(paste0(fig_dir,"07_KM_riskgroup.tif"), width=6, height=6, units="in",
     res=300, pointsize=8, compression="lzw", bg="white", type="cairo")
print(km); dev.off()
cat("KM log-rank p =", signif(surv_pvalue(sfit,data=dat)$pval,3), "\n")

# 5) timeROC：1/3/5 年 AUC
tROC <- timeROC(T=dat$OS.time, delta=dat$OS, marker=dat$riskScore,
                cause=1, times=c(12,36,60), iid=TRUE)
auc <- tROC$AUC; cat("AUC 1/3/5y:", paste(round(auc,3),collapse=" / "), "\n")
rocdf <- do.call(rbind, lapply(seq_along(c(12,36,60)), function(i){
  data.frame(FP=tROC$FP[,i], TP=tROC$TP[,i],
             yr=paste0(c(1,3,5)[i],"-year (AUC=",sprintf("%.2f",auc[i]),")"))
}))
p_roc <- ggplot(rocdf, aes(FP, TP, color=yr)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60", linewidth=LW) +
  geom_path(linewidth=LW*2) +
  scale_color_manual(values=COL7[c(1,2,3)]) +
  labs(title="Time-dependent ROC (TCGA-COAD)", x="1 - Specificity",
       y="Sensitivity", color=NULL) +
  coord_equal() + theme_pub_bw() + theme(legend.position=c(0.65,0.2))
save_fig(p_roc, "07_timeROC", width=5, height=5)

# 风险因子三联图：评分排序 + 生存状态散点 + 基因热图
ord <- order(dat$riskScore); dord <- dat[ord, ]
rk <- data.frame(idx=seq_len(nrow(dord)), score=dord$riskScore,
                 group=dord$riskGroup, time=dord$OS.time,
                 status=factor(dord$OS,labels=c("Alive","Dead")))
p_rank <- ggplot(rk, aes(idx, score, color=group)) + geom_point(size=0.6) +
  geom_vline(xintercept=sum(dord$riskGroup=="Low"), linetype="dashed",
             color="grey60", linewidth=LW) +
  scale_color_manual(values=c(Low=COL7[1],High=COL7[2])) +
  labs(title="Risk score", x=NULL, y="Risk score", color=NULL) + theme_pub_bw()
p_scat <- ggplot(rk, aes(idx, time, color=status)) + geom_point(size=0.6) +
  geom_vline(xintercept=sum(dord$riskGroup=="Low"), linetype="dashed",
             color="grey60", linewidth=LW) +
  scale_color_manual(values=c(Alive=COL7[1],Dead=COL7[2])) +
  labs(x="Patients (ranked by risk)", y="Survival (months)", color=NULL) +
  theme_pub_bw()
save_fig(p_rank, "07_riskplot_score", width=6, height=3)
save_fig(p_scat, "07_riskplot_scatter", width=6, height=3)

# 6) 多因素 Cox + 列线图（风险评分 + 临床协变量）
ndf <- data.frame(OS.time=dat$OS.time, OS=dat$OS, riskScore=dat$riskScore)
if (!is.null(clin)) {
  cc <- clin[rownames(dat), ]
  ndf$age <- suppressWarnings(as.numeric(cc$age))
  st <- toupper(as.character(cc$stage))
  ndf$stage <- ifelse(grepl("IV", st), "III-IV",
                ifelse(grepl("III", st), "III-IV", "I-II"))
  ndf$stage[is.na(st) | st=="NA"] <- NA
  ndf$stage <- factor(ndf$stage, levels=c("I-II","III-IV"))
}
ndf2 <- ndf[complete.cases(ndf), ]
cat("列线图建模样本:", nrow(ndf2), "（去除协变量缺失）\n")
dd <- datadist(ndf2); options(datadist="dd")
vars <- intersect(c("riskScore","age","stage"), names(ndf2))
fml <- as.formula(paste("Surv(OS.time,OS) ~", paste(vars, collapse=" + ")))
cph <- cph(fml, data=ndf2, surv=TRUE, x=TRUE, y=TRUE)
# 多因素 Cox 表
mcox <- summary(coxph(fml, data=ndf2))
write.csv(data.frame(var=rownames(mcox$coef), HR=mcox$coef[,"exp(coef)"],
  p=mcox$coef[,"Pr(>|z|)"], lower=mcox$conf.int[,"lower .95"],
  upper=mcox$conf.int[,"upper .95"]),
  paste0(tab_dir,"07_multivariate_cox.csv"), row.names=FALSE)

surv_fun <- Survival(cph)
nom <- nomogram(cph, fun=list(function(x) surv_fun(12,x),
                              function(x) surv_fun(36,x),
                              function(x) surv_fun(60,x)),
                funlabel=c("1-year OS","3-year OS","5-year OS"),
                fun.at=seq(0.1,0.9,0.2))
save_base("07_nomogram", function() { plot(nom, cex.axis=0.7, cex.var=0.8) },
          width=9, height=6)

# 校准曲线（3 年，可改 u）
cph3 <- cph(fml, data=ndf2, surv=TRUE, x=TRUE, y=TRUE, time.inc=36)
cal <- calibrate(cph3, cmethod="KM", method="boot", u=36,
                 m=max(20, floor(nrow(ndf2)/4)), B=200)
save_base("07_calibration_3y", function() {
  par(mar=c(4,4,2,1), cex=0.8, font.main=1, font.lab=1, font.axis=1)
  plot(cal, xlab="Nomogram-predicted 3-year OS",
       ylab="Observed OS", subtitles=FALSE)
  title("Calibration (3-year)", cex.main=0.9)
}, width=5, height=5)

# C-index
cidx <- rcorrcens(Surv(OS.time,OS) ~ predict(cph), data=ndf2)
cat("C-index =", round(1-cidx[1],3), "\n")
saveRDS(dat, "data/tcga_risk_model.rds")
cat("✓ 07 预后模型完成\n")
