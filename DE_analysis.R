library(SummarizedExperiment)
library(ggplot2)
library(edgeR)

#Introduce data
se <- readRDS( "seLUAD.rds")
se

#Process dge list
dge_luad <- DGEList(counts = assays(se)$counts, genes = as.data.frame(mcols(se)), group = se$type)

#CPM scaling
CPM <- t(t(dge_luad$counts)/(dge_luad$samples$lib.size/1e+06))
assays(se)$logCPM <- cpm(dge_luad, log = TRUE, prior.count = 0.25)
assays(se)$logCPM[1:3, 1:7]

#Expression level by sample
library(geneplotter)
#multidensity(assays(se)$logCPM, xlab = "log2 CPM", legend = NULL, 
#             main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

# Distribution of expression across genes
avgexp <- rowMeans(assays(se)$logCPM)
hist(avgexp, xlab = expression(log[2] * "CPM"), main = "", las = 1, col = "gray")

#Filtering genes
cpmcutoff <- round(10/min(dge_luad$sample$lib.size/1e+06), digits = 1)
cpmcutoff

nsamplescutoff <- min(table(se$type))
nsamplescutoff

mask <- rowSums(cpm(dge_luad) > cpmcutoff) >= nsamplescutoff

se.filt <- se[mask, ]
dge_luad.filt <- dge_luad[mask, ]
dim(se.filt)

h <- hist(avgexp, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "", 
          las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(se.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))

### MA plots

# Vuvuzela plot
plotSmear(dge_luad.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

# Vuvuzela plot refiltered
mask <- rowMeans(assays(se)$logCPM) > 1
se.filt <- se[mask, ]
dge_luad.filt <- dge_luad[mask, ]
dim(se.filt)

plotSmear(dge_luad.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

#Vuvuzela plot normalized (TMM)

dge.filt <- calcNormFactors(dge.filt, normalize.method="quantile")

par(mfrow = c(1, 2))
plotSmear(dge_luad, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)
plotSmear(dge_luad.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

# Expression by sample (two first samples)
par(mfrow=c(1, 2), mar=c(4, 5, 3, 1))
for (i in 1:2) {
  A <- rowMeans(assays(lclse.filt)$logCPM) ; M <- assays(lclse.filt)$logCPM[, i] - A
  smoothScatter(A, M, main=colnames(lclse.filt)[i], las=1, cex.axis=1.2, cex.lab=1.5, cex.main=2)
  abline(h=0, col="blue", lwd=2) ; lo <- lowess(M ~ A) ; lines(lo$x, lo$y, col="red", lwd=2)
}

## MDS
plotMDS(dge_luad.filt, col = c("red", "blue")[as.integer(dge_luad.filt$samples$group)], cex = 6, pch = 22)
legend("topleft", c("Normal", "Tumor"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)

# faltaria fer nomes els punts sense nom
# faltaria fer per els altres 3 casos comprobats amb el GGplot