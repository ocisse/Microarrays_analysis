module load R/3.2.1-goolf-1.7.20
R

# load the saved session
library("affy")
library(limma)
library(xps)
library(oligo)
library(pd.mogene.2.0.st)
library(mas5)
library(gcrma)
library(GEOquery)
library(simpleaffy)
library(RColorBrewer)
library(AnnotationForge)
library(mouse.db0)
library(org.Mm.eg.db)
library(annotate)
library(mogene20sttranscriptcluster.db)

# --- load data
targets_ccmd <- readTargets("affy_targets.txt")
data <- read.celfiles(filenames=targets_ccmd$FileName)
est <-rma(data)
annotation(est) <- "mogene20sttranscriptcluster.db"

# ---  QC
pdf("boxUnnormalized.pdf")
boxplot(data, which=c("all"), transfo=log2)
dev.off()

pdf("boxNormalized.pdf")
boxplot(est, transfo=identity)
dev.off()

pdf("hist_density_vs_log.unnormalized.pdf")
hist(data, which=c("all"), transfo=log2)
dev.off()

pdf("hist_density_vs_log.normalized.pdf")
hist(est, transfo=identity)
dev.off()

fit <- fitProbeLevelModel(data)
pdf("fit.pdf")
image(fit)
dev.off()

pdf("NUSE.pdf")
NUSE(fit)
dev.off()

pdf("RLE.pdf")
RLE(fit)
dev.off()

eset  <- exprs(est)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
pdf("clusters.pdf")
plot(clusters)
dev.off()

pc = prcomp(t(exprs(est)), scale.=TRUE)
pdf("2ndScan_PCA.pdf")
#par(mfrow = c(1, 2))
color=c('green','green','green','green',
        'red','red','red','red',
        'blue','blue','blue','blue',
        'darkorange','darkorange','darkorange','darkorange',
        'gray4','gray4','gray4','gray4',
        'magenta4','magenta4','magenta4','magenta4',
        'seagreen4','seagreen4','seagreen4','seagreen4',
        'darkviolet','darkviolet','darkviolet','darkviolet')

plot(pc$x[,1:2 ], xlab = "PC1", ylab = "PC2",
#       pch = unclass(as.factor(pData(est)[, 1])),
#       col = unclass(as.factor(pData(est)[,2])),
        main = "2nd Scan PCA",
        col=color)
dev.off()
summary(pc)


# --- FILETERING DATA
library(genefilter)
celfiles.filtered <- nsFilter(est, require.entrez=FALSE, remove.dupEntrez=FALSE)
# What got removed and why?
celfiles.filtered$filter.log

# Finding differentially expressed probesets
samples <- targets_ccmd$Name
# check the results of this
samples
# convert into factors
samples <- as.factor(samples)
# check factors have been assigned
design <- model.matrix(~0 + samples)
colnames(design) <- c("KO_d120", "KO_d35", "KO_d45", "KO_d75","WT_d120", "WT_d35", "WT_d45", "WT_d75")

library(limma)

# fit the linear model to the filtered expression set
fit <- lmFit(exprs(celfiles.filtered$eset), design) # stop here
contrast.matrix_wt_d75_vs_ko <-  makeContrasts( WT_d75_KO_d45 = WT_d75 - KO_d45, WT_d75_KO_d75 = WT_d75 - KO_d75, WT_d75_KO_d120 = WT_d75 - KO_d120, levels=design)
WT_d75_vs_KO_fits <-contrasts.fit(fit,contrast.matrix_wt_d75_vs_ko)
WT_d75_vs_KO_ebFit <- eBayes(WT_d75_vs_KO_fits)
probeset.list_WT_d75_vs_KO <- topTable(WT_d75_vs_KO_ebFit, coef=2, p.value=0.001, adjust.method="fdr", number=10000, lfc=1.5)

# annotation
gene.symbols_WT_d75_vs_KO <- getSYMBOL(rownames(probeset.list_WT_d75_vs_KO), "mogene20sttranscriptcluster.db")
results_WT_d75_vs_KO <- cbind(probeset.list_WT_d75_vs_KO, gene.symbols_WT_d75_vs_KO)
write.table(results_WT_d75_vs_KO, "results_WT_d75_vs_KO.txt", sep="\t", quote=FALSE)
ids_KO <- rownames(probeset.list_KO)
ann_KO <- select(mogene20sttranscriptcluster.db, ids_KO, c("ENTREZID","GENENAME"), "PROBEID")
write.table(ann_KO,"results.probes_AND_Names_KO.txt", sep="\t", quote=FALSE)

ids_WT_d75_vs_KO <- rownames(probeset.list_WT_d75_vs_KO)
ann_WT_d75_vs_KO <- select(mogene20sttranscriptcluster.db, ids_WT_d75_vs_KO , c("ENTREZID","GENENAME"), "PROBEID")
write.table(ann_WT_d75_vs_KO,"results_WT_d75_vs_KO.ANN.txt", sep="\t", quote=FALSE)

# merge
perl merged_result_AND_annotations.pl results_WT_d75_vs_KO.txt results_WT_d75_vs_KO.ANN.txt


## --- Paired Two group analysis analysis + HM

# Heatmap
library("gplots")

wt <- est[, est$varLabels %in% c("WT_d75","IL-21R_KO_d75")]
f_wt <- factor(as.character(as.factor(wt$varLabels)))
design_wt <- model.matrix(~f_wt)
fit_wt <- eBayes(lmFit(wt,design_wt))
topTable(fit_wt, coef=2, p.value=0.001, adjust.method="fdr", number=10000, lfc=1.5)) 

selected  <- p.adjust(fit_wt$p.value[, 2], method ="fdr", n = length(fit_wt$p.value[, 2])) < 0.001
wt_2fc_Pval <- wt [selected, ]

color.map <- function(varLabels) { if (varLabels=="WT_d35") "#FF0000" else "#0000FF" }
micecolors <- unlist(lapply(wt_2fc_Pval$varLabels, color.map))

pdf("wt_d35_vs_wt_d45_2fc_Pval_0.001.pdf")
heatmap.2(exprs(wt_2fc_Pval), col=redgreen(75), scale="row", ColSideColors=micecolors,
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
dev.off()
