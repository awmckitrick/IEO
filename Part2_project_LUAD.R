#################
## Charge modules
#################
library(knitr)
library(SummarizedExperiment)
library(edgeR)
library(ggplot2)
library(plyr)
library(sva)
library(grid)
library(geneplotter)
library(limma)
library(corpcor)

##################
## Initialize data
##################

# Charge data
se <- readRDS( "seLUAD.rds")
#Charge sample variables
dge_luad <- DGEList(counts = assays(se)$counts, genes = as.data.frame(mcols(se)), group = se$type)

# Logaritmize
CPM <- t(t(dge_luad$counts)/(dge_luad$samples$lib.size/1e+06))
assays(se)$logCPM <- cpm(dge_luad, log = TRUE, prior.count = 0.25)
assays(se)$logCPM[1:3, 1:7]

# Apply filters
cpmcutoff <- round(10/min(dge_luad$sample$lib.size/1e+06), digits = 1)
cpmcutoff

nsamplescutoff <- min(table(se$type))
nsamplescutoff

mask <- rowSums(cpm(dge_luad) > cpmcutoff) >= nsamplescutoff

se.filt <- se[mask, ]
dge_luad.filt <- dge_luad[mask, ]
dim(se.filt)

# Normalize
dge_luad.filt <- calcNormFactors(dge_luad.filt, normalize.method="quantile")

# Batch effect variables and add them to the summarized experiment
tss <- substr(colnames(se.filt), 6, 7)
names(tss) <- colnames(se.filt)
se.filt$tss <- tss

center <- substr(colnames(se.filt), 27, 28)
names(center) <- colnames(se.filt)
se.filt$center <- center

plate <- substr(colnames(se.filt), 22, 25)
names(plate) <- colnames(se.filt)
se.filt$plate <- plate

portionanalyte <- substr(colnames(se.filt), 18, 20)
names(portionanalyte) <- colnames(se.filt)
se.filt$portionanalyte <- portionanalyte

samplevial <- substr(colnames(se.filt), 14, 16)
names(samplevial) <- colnames(se.filt)
se.filt$samplevial <- samplevial

gender <- unname(se.filt$gender)
names(gender) <- colnames(se.filt)
se.filt$gender <- gender

race <- unname(se.filt$race)
histo <- unname(se.filt$histologic_diagnosis.1)

# Removing batch effect
#logCPM <- cpm(dge_luad.filt, log=TRUE, prior.count = 3)
#batch <- as.integer(factor(samplevial))
#logCPM.batch <- logCPM
#logCPM <-removeBatchEffect(logCPM, batch, design = mod)
#assays(se.filt)$logCPM <- logCPM
#s <- fast.svd(t(scale(t(logCPM), center = TRUE, scale = TRUE)))
#pcSds <- s$d
#pcSds[1] <- 0
#svdexp <- s$u %*% diag(pcSds) %*% t(s$v)
#colnames(svdexp) <- colnames(se.filt)
#d <- as.dist(1-cor(svdexp, method="spearman"))
#sampleClustering <- hclust(d)
#sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
#names(outcome) <- colnames(se.filt)

# Filter by logCPM
mask <- rowMeans(assays(se.filt)$logCPM) > 1
sum(mask)

se.filt <- se.filt[mask, ]
dim(se.filt)

dge_luad.filt <- dge_luad.filt[mask, ]
dim(dge_luad.filt)

##########################
## DIfferential expression
##########################

# 1. Create SVA models
mod <- model.matrix(~ se.filt$type + se.filt$samplevial, colData(se.filt))
mod0 <- model.matrix(~ se.filt$samplevial, colData(se.filt))
pv <- f.pvalue(assays(se.filt)$logCPM, mod, mod0)
sum(p.adjust(pv, method = "bonferroni") < 0.01)

sv <- sva(assays(se.filt)$logCPM, mod, mod0)
sv$n

modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(se.filt)$logCPM, mod, mod0)
sum(p.adjust(pvsv, method="bonferroni") < 0.01)

par(mfrow = c(1,2))
v <- voom(dge_luad.filt, mod, plot=TRUE)
vsv <- voom(dge_luad.filt, modsv, plot=TRUE)

#2. Fit linear model
par(mfrow = c(1,1))
fit5 <- lmFit(v, mod)
fit5sv <- lmFit(vsv,modsv)

#3. T-statistics
fit5 <- eBayes(fit5)
fit5sv <- eBayes(fit5sv)

#4. Extend of differencial expression
FDRcutoff <- 0.05
res5 <- decideTests(fit5, p.value = FDRcutoff)
summary(res5)
res5sv <- decideTests(fit5sv, p.value = FDRcutoff)
summary(res5sv)

#5. Metadata and fetch table
genesmd <- data.frame(
  chr = as.character(seqnames(rowRanges(se.filt))),
  symbol = rowData(se.filt)[, 1],
  stringsAsFactors = FALSE) # We have to say explicitly metadata is not a factor
fit5$genes <- genesmd
tt5 <- topTable(fit5, coef = 2, n = Inf)
head(tt5)

genesmd <- data.frame(
  chr = as.character(seqnames(rowRanges(se.filt))),
  symbol = rowData(se.filt)[, 1],
  stringsAsFactors = FALSE) # We have to say explicitly metadata is not a factor
fit5sv$genes <- genesmd
tt5sv <- topTable(fit5sv, coef = 2, n = Inf)
head(tt5sv)

#6. Chromosome distribution
sort(table(tt5$chr[tt5$adj.P.Val < FDRcutoff]), decreasing = TRUE)
sort(table(tt5sv$chr[tt5sv$adj.P.Val < FDRcutoff]), decreasing = TRUE)

#8. Diagnostic plots
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(tt5sv$P.Value, xlab = "Raw P-values", main = "", las = 1)
qqt(fit5$t[, 2], df = fit5$df.prior + fit5$df.residual, main = "", pch = ".", cex = 3)
abline(0, 1, lwd = 2)

pvsv <- f.pvalue(assays(se.filt)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="bonferroni") < 0.01)

par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(tt5sv$P.Value, xlab = "Raw P-values", main = "", las = 1)
qqt(fit5sv$t[, 2], df = fit5sv$df.prior + fit5sv$df.residual, main = "", pch = ".", cex = 3)
abline(0, 1, lwd = 2)

pvsv <- f.pvalue(assays(se.filt)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="bonferroni") < 0.01)


# 9. Volcano plot
par(mfrow = c(1,1))
volcanoplot(fit5, coef = 2,
            highlight = 10,
            names = fit5$genes$symbol,
            main = "Model fit5",
            las = 1)
volcanoplot(fit5sv, coef = 2,
            highlight = 10,
            names = fit5sv$genes$symbol,
            main = "Model fit5 sva",
            las = 1)

# Looking for replicated data
table(table(se.filt$tissue_source_site))

################################
## Gene enrichment (Fisher test)
################################

# 0. Prepaing gene universe
geneUniverse <- rownames(se)
DEgenes <- rownames(tt5)[tt5$adj.P.Val < 0.01 & tt5$logFC > 2]

# 1. Prepare parameter used
library(GOstats)
conditional(params) <- TRUE
params <- new("GOHyperGParams",
              geneIds=DEgenes,
              universeGeneIds=geneUniverse,
              annotation="org.Hs.eg.db",
              ontology="BP",
              pvalueCutoff=0.05,
              testDirection="over")

# 2. Functional enrichment analysis
hgOver <- hyperGTest(params)
hgOver

#3. Store and visualize results
htmlReport(hgOver, file = "gotests_overtumor.html")#Not working for version problems
browseURL("gotests_overtumor.html")

summary(hgOver)

###### Intepreting table
#Check
goresults <- summary(hgOver)
head(goresults)

#Filter to size between 3-300 and Mmore than 3 counts
goresults <- goresults[goresults$Size >= 3 & goresults$Size <= 300 & goresults$Count >= 3, ]
goresults <- goresults[order(goresults$OddsRatio, decreasing=TRUE), ]
head(goresults)

# Extract genes that enrich GO terms and paste them to results
geneIDs <- geneIdsByCategory(hgOver)[goresults$GOBPID]
geneSYMs <- sapply(geneIDs, function(id) select(org.Hs.eg.db, columns = "SYMBOL", key = id, 
                                                keytype = "ENTREZID")$SYMBOL)
geneSYMs <- sapply(geneSYMs, paste, collapse = ", ")
goresults <- cbind(goresults, Genes = geneSYMs)
rownames(goresults) <- 1:nrow(goresults)

library(xtable)
xtab <- xtable(goresults, align = "l|c|r|r|r|r|r|p{3cm}|p{3cm}|")
print(xtab, file = "goresults_overtumor.html", type = "html")
browseURL("goresults_overtumor.html")

#######Provisional checking volcano plot
par(mfrow = c(1,1))
volcanoplot(fit5sv, coef = 2,
            highlight = 7,
            col = ifelse(fit5sv$genes$symbol %in% c("EPAS1", "ABCA12", "HOMER2", "HOMER1" ),"red","blue"),
            names = fit5sv$genes$symbol,
            cex = ifelse(fit5sv$genes$symbol %in% c("EPAS1", "ABCA12", "HOMER2", "HOMER1"),1,0.1),
            main = "Model fit5 sva",
            las = 1)
#######
## GSVA
#######

# Charge gene setss
library(GSVAdata)
data("c2BroadSets")
head(names(c2BroadSets))
gsc <- GeneSetCollection(c(c2BroadSets))

# 0. Map identifiers
gsc <- mapIdentifiers(gsc, AnnoOrEntrezIdentifier(metadata(se.filt)$annotation))
gsc

# 1. Incidence matrix: which genes belong to each dataset
Im <- incidence(gsc)
dim(Im)
Im[1:6, 1:10]

# Discard genes not part of our data
Im <- Im[, colnames(Im) %in% rownames(se.filt)]
dim(Im)

# Discard genes not present in any of the datasets
se.filt <- se.filt[colnames(Im), ]
dim(se.filt)
dge_luad.filt <- dge_luad.filt[colnames(Im), ]
dim(dge_luad.filt)

# Limit to gene sets with at least 5 genes, to give more robustness to our statistics
Im <- Im[rowSums(Im) >= 5, ]
dim(Im)

#2. Make sure the DE analysis is already done

#3 T-values plot
par(mfrow = c(1,1))
qq <- qqnorm(tt5sv$t)
abline(0, 1)
chroutliers <- tt5sv$chr[abs(tt5sv$t) > 10]
text(qq$x[abs(qq$y) > 10], qq$y[abs(qq$y) > 10], chroutliers, pos = 4)# So wrong

#7. Let's try normal: Calculate Z-scores and qqplot
tGSgenes <- tt5sv[match(colnames(Im), rownames(tt5sv)), "t"]
zS <- sqrt(rowSums(Im)) * (as.vector(Im %*% tGSgenes)/rowSums(Im))
qqnorm(zS)
abline(0, 1)

#8. Sorting by z-scores to obtain best gene sets
rnkGS <- sort(abs(zS), decreasing = TRUE)
head(rnkGS)

#9. Scatter plot

plotGS <- function(se, gs, pheno, ...) {
  l <- levels(colData(se)[, pheno])
  idxSamples1 <- colData(se)[, pheno] == l[1]
  idxSamples2 <- colData(se)[, pheno] == l[2]
  exps1 <- rowMeans(assays(se)$logCPM[gs, idxSamples1])
  exps2 <- rowMeans(assays(se)$logCPM[gs, idxSamples2])
  rng <- range(c(exps1, exps2))
  plot(exps1, exps2, pch = 21, col = "black", bg = "black", xlim = rng, ylim = rng,
       xlab = l[1], ylab = l[2], ...)
  abline(a = 0, b = 1, lwd = 2, col = "red")
}

genesGS1 <- colnames(Im)[which(Im[names(rnkGS)[1], ] == 1)]
genesGS2 <- colnames(Im)[which(Im[names(rnkGS)[2], ] == 1)]
par(mfrow = c(1, 2), mar = c(4, 5, 3, 4))
plotGS(se.filt, genesGS1, "type", main = names(rnkGS)[1], cex.lab = 2, las = 1)
plotGS(se.filt, genesGS2, "type", main = names(rnkGS)[2], cex.lab = 2, las = 1)


#10. One-sample normal distribution
pv <- pmin(pnorm(zS), 1 - pnorm(zS))
sum(pv < 0.05)

pvadj <- p.adjust(pv, method = "fdr")
DEgs <- names(pvadj)[which(pvadj < 0.01)]
length(DEgs)

head(DEgs, n = 3)

# 11. Overlap of gene sets check. Build a table ranking pairs of more overlapped gene sets
library(GSVA)
gsov <- computeGeneSetsOverlap(gsc[DEgs], rownames(se.filt))

trimask <- upper.tri(gsov)
rnkOv <- data.frame(gs1 = row(gsov)[trimask], gs2 = col(gsov)[trimask], ov = gsov[trimask])
rnkOv <- rnkOv[order(rnkOv$ov, decreasing = TRUE), ]
rnkOv$gs1 <- rownames(gsov)[rnkOv$gs1]
rnkOv$gs2 <- rownames(gsov)[rnkOv$gs2]
sum(rnkOv$ov == 1)# How many gene sets are identical?

# 12. Change in scale: to avoid positive and negative values of one gene set to anulate each other
library(Category)
xS <- applyByCategory(tGSgenes,
                      Im,
                      function(x) (sum((x - mean(x))^2) - (length(x) - 1))/(2 * (length(x) - 1)))

rnkGS <- sort(abs(xS), decreasing = TRUE)

# 13. Obtain p-values and check results

pv <- pmin(pnorm(xS), 1 - pnorm(xS))
pvadj <- p.adjust(pv)
DEgsByScale <- names(pvadj)[which(pvadj < 0.01)]
length(DEgsByScale)

length(intersect(DEgs, DEgsByScale))

head(setdiff(DEgsByScale, DEgs))# Whoa! many cancer, so wow

#14. Top 3 plots
topgs1genes <- colnames(Im)[which(Im[names(rnkGS)[1], ] == 1)]
topgs2genes <- colnames(Im)[which(Im[names(rnkGS)[2], ] == 1)]
topgs3genes <- colnames(Im)[which(Im[names(rnkGS)[3], ] == 1)]
par(mfrow = c(1, 3))
plotGS(se.filt, topgs1genes, "type", main = names(rnkGS)[1], cex.lab = 2, las = 1)
plotGS(se.filt, topgs2genes, "type", main = names(rnkGS)[2], cex.lab = 2, las = 1)
plotGS(se.filt, topgs3genes, "type", main = names(rnkGS)[3], cex.lab = 2, las = 1)

############# Gene set Variation analysis

# 1. Create Gene set expresssion matrix (with gene sets instead of genes)
library(GSVA)
GSexpr <- gsva(assays(se.filt)$logCPM, gsc,
               min.sz=5, max.sz=300, verbose=FALSE)
class(GSexpr)

#2. Perform a regular SVA-DE pipeline
svaobj <- sva(GSexpr, mod, mod0)
modSVs <- cbind(mod, svaobj$sv)

fit <- lmFit(GSexpr, modSVs)
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, n = Inf)
DEgs <- rownames(tt[tt$adj.P.Val < 0.01, , drop = FALSE])
DEgs

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            # 3. Inspect results with boxplot. CHange MSY and XiE for the ones you desire to analyze
par(mfrow = c(1, 2))
#boxplot(GSexpr["TRIM29", ] ~ lclse$type, main = "MSY", las = 1, cex.axis = 2)
#boxplot(GSexpr["XiE", ] ~ lclse$type, main = "XiE", las = 1, cex.axis = 2)

# 4. Volcano plot time!!
par(mfrow = c(1,1))
plot(tt5sv$logFC, -log10(tt5sv$P.Value), xlab="Log2 fold-change", ylab="-log10 P-value",
     pch=".", cex=5, col=grey(0.75), cex.axis=1.2, cex.lab=1.5, las=1)
posx <- tt5sv[tt5sv$adj.P.Val < 0.01, "logFC"] ; posy <- -log10(tt5sv[tt5sv$adj.P.Val < 0.01, "P.Value"])
points(posx, posy, pch=".", cex=5, col="red")
text(posx, posy, tt5sv$symbol[tt5sv$adj.P.Val < 0.01], pos=1)


############# Extra!!!!!! How many factors contain information about normal samples
# Response: Only type, gender and patient_barcode. So we can study nothing about covariates

for (fac in names(colData(se))){
  if(length(table(table(se[[fac]],se$type)[,1])) > 1){
    print(fac)
    print(table(se[[fac]],se$type))
  }
}

