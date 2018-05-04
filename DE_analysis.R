library(SummarizedExperiment)
library(edgeR)
library(ggplot2)
library(plyr)
library(sva)
library(grid)
library(geneplotter)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Introduce data
se <- readRDS( "seLUAD.rds")
se

#Process dge list
dge_luad <- DGEList(counts = assays(se)$counts, genes = as.data.frame(mcols(se)), group = se$type)

#CPM scaling
CPM <- t(t(dge_luad$counts)/(dge_luad$samples$lib.size/1e+06))
assays(se)$logCPM <- cpm(dge_luad, log = TRUE, prior.count = 0.25)
assays(se)$logCPM[1:3, 1:7]

# BARPLOTS
ord <- order(dge_luad$sample$lib.size)
i <- 0
ordsliced <- list(0)
gender <- list(0)
race <- list(0)
typ <- list(0)
hstory_malignancy <- list(0)
for(num in 1:6){
  ordsliced[[num]] <- c(ord[(i*100 +1):(num*100)])
  i <- num
  g_nas <- length(colData(se)$gender[ordsliced[[num]]][colData(se)$gender[ordsliced[[num]]] == "NA"])
  gender[[num]] <- c(length(colData(se)$gender[ordsliced[[num]]][colData(se)$gender[ordsliced[[num]]] == "MALE"])-g_nas,
                     length(colData(se)$gender[ordsliced[[num]]][colData(se)$gender[ordsliced[[num]]] == "FEMALE"])-g_nas, g_nas)
  
  race_nas <- length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "NA"])
  race_tnas <- length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "[Not Available]"]) +
    length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "[Not Evaluated]"]) +
    length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "[Unknown]"]) -2*race_nas
  race[[num]] <- c(length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "AMERICAN INDIAN OR ALASKAN NATIVE"]) - race_nas,
                   length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "ASIAN"]) - race_nas,
                   length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "BLACK OR AFRICAN AMERICAN"]) - race_nas,
                   length(colData(se)$race[ordsliced[[num]]][colData(se)$race[ordsliced[[num]]] == "WHITE"]) - race_nas,
                   race_tnas)
  hist_nas <- length(colData(se)$history_other_malignancy[ordsliced[[num]]][colData(se)$history_other_malignancy[ordsliced[[num]]] == "NA"])
  hstory_malignancy[[num]] <- c(length(colData(se)$history_other_malignancy[ordsliced[[num]]][colData(se)$history_other_malignancy[ordsliced[[num]]] == "Yes"]) - hist_nas +
                                  length(colData(se)$history_other_malignancy[ordsliced[[num]]][colData(se)$history_other_malignancy[ordsliced[[num]]] == "Yes, History of Prior Malignancy"]) - hist_nas + 
                                  length(colData(se)$history_other_malignancy[ordsliced[[num]]][colData(se)$history_other_malignancy[ordsliced[[num]]] == "Yes, History of Synchronous and or Bilateral Malignancy"]) - hist_nas,
                                length(colData(se)$history_other_malignancy[ordsliced[[num]]][colData(se)$history_other_malignancy[ordsliced[[num]]] == "No"]) - hist_nas,
                                hist_nas)
  typ_na <- length(colData(se)$type[ordsliced[[num]]][colData(se)$type[ordsliced[[num]]] == "NA"])
  typ[[num]] <- c(length(colData(se)$type[ordsliced[[num]]][colData(se)$type[ordsliced[[num]]] == "normal"]) - typ_na,
                  length(colData(se)$type[ordsliced[[num]]][colData(se)$type[ordsliced[[num]]] == "tumor"]) - typ_na,
                  typ_na)
}

# Gender BARPLOT
g_df <- data.frame(supp=c("MALE", "FEMALE", "NA"),
                   depth=rep(c('1st', '2nd', '3rd', '4th', '5th', '6th'),each = 3),
                   freq = c(gender[[1]],gender[[2]],gender[[3]], 
                            gender[[4]], gender[[5]],gender[[6]]))
df_sorted <- arrange(g_df, depth, supp) 
df_cumsum <- ddply(df_sorted, "depth",transform, label_ypos=cumsum(freq))

p1 <- ggplot( data=df_cumsum, aes(x=depth, y=freq, fill=supp)) +
  geom_bar(stat="identity") +
  geom_text(aes(y=107.5-label_ypos,label=freq), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5) +
  theme_minimal() + ggtitle("Barplot Gender depth distribution") +
  theme(plot.title = element_text(hjust = 0.5))

# race BARPLOT
r_df <- data.frame(supp = c("A. NATIVE", "ASIAN","BLACK", "WHITE", "NA"),
                   depth=rep(c('1st', '2nd', '3rd', '4th', '5th', '6th'),each = 5),
                   freq = c(race[[1]],race[[2]],race[[3]], race[[4]], race[[5]],race[[6]]))
df_sorted <- arrange(r_df, depth, supp) 
df_cumsum <- ddply(df_sorted, "depth",transform, label_ypos=cumsum(freq))

p2 <-ggplot( data=df_cumsum, aes(x=depth, y=freq, fill=supp)) +
  geom_bar(stat="identity") +
  geom_text(aes(y=107.5-label_ypos,label=freq), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5) +
  theme_minimal() + ggtitle("Barplot Race depth distribution") +
  theme(plot.title = element_text(hjust = 0.5))

# history_other_malignanci BAPLOT
h_df <- data.frame(supp = c("Yes", "No", "NA"),
                   depth=rep(c('1st', '2nd', '3rd', '4th', '5th', '6th'),each = 3),
                   freq = c(hstory_malignancy[[1]],hstory_malignancy[[2]],hstory_malignancy[[3]], hstory_malignancy[[4]], hstory_malignancy[[5]],hstory_malignancy[[6]]))
df_sorted <- arrange(h_df, depth, supp) 
df_cumsum <- ddply(df_sorted, "depth",transform, label_ypos=cumsum(freq))

p3 <- ggplot( data=df_cumsum, aes(x=depth, y=freq, fill=supp)) +
  geom_bar(stat="identity") +
  geom_text(aes(y=107.5-label_ypos,label=freq), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5) +
  theme_minimal() + ggtitle("Barplot history of other malignancy depth distribution") +
  theme(plot.title = element_text(hjust = 0.5))

#Type
t_df <- data.frame(supp = c("Normal", "Tumor", "NA"),
                   depth=rep(c('1st', '2nd', '3rd', '4th', '5th', '6th'),each = 3),
                   freq = c(typ[[1]],typ[[2]],typ[[3]], typ[[4]], typ[[5]],typ[[6]]))
df_sorted <- arrange(t_df, depth, supp) 
df_cumsum <- ddply(df_sorted, "depth",transform, label_ypos=cumsum(freq))

p4 <-ggplot( data=df_cumsum, aes(x=depth, y=freq, fill=supp)) +
  geom_bar(stat="identity") +
  geom_text(aes(y=107.5-label_ypos,label=freq), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5) +
  theme_minimal() + ggtitle("Barplot of type depth distribution") +
  theme(plot.title = element_text(hjust = 0.5))

multiplot(p1, p2, p3, p4, cols=2)

#Expression level by sample
par(mfrow=c(1, 2))
multidensity(as.list(as.data.frame(assays(se[, se$type == "tumor"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Tumor samples", las=1)
multidensity(as.list(as.data.frame(assays(se[, se$type == "normal"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Normal samples", las=1)

# Distribution of expression across genes
avgexp <- rowMeans(assays(se)$logCPM)
par(mfrow=c(1,1))
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

png(filename = "./img/expression_by_gene_corrected.png",width = 600, height = 500)
par(mar = c(5.1,5.1,5.1,2.1))
h <- hist(avgexp, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "", 
          las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(se.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))
dev.off()

### MA plots

# Vuvuzela plot
par(mar = c(5.1,5.1,4.1,2.1),mfrow = c(1,1))
plotSmear(dge_luad, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

# Vuvuzela plot refiltered
mask <- rowMeans(assays(se)$logCPM) > 1
se.filt <- se[mask, ]
dge_luad.filt <- dge_luad[mask, ]
dim(se.filt)

plotSmear(dge_luad.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

#Vuvuzela plot normalized (TMM)

dge_luad.filt <- calcNormFactors(dge_luad.filt, normalize.method="quantile")

png("./img/vuvuzelas.png",height = 500, width = 1000)
par(mar = c(5.1,4.1,4.1,2.1),mfrow = c(1, 2))
plotSmear(dge_luad, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "Preprocesed data")
abline(h = 0, col = "blue", lwd = 2)
plotSmear(dge_luad.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "Filtered and normalized data")
abline(h = 0, col = "blue", lwd = 2)
dev.off()

# Expression by sample (tumor in blue, normal in green)
par(mfrow=c(5, 4), mar = c(2,2,2,2))
for (i in c(1,581)) {
  inici <- i
  final <- i+19
  for (i in inici:final) {
    A <- rowMeans(assays(se.filt)$logCPM) ; M <- assays(se.filt)$logCPM[, i] - A; C <- se.filt$type[[i]];
    print(C)
    smoothScatter(A, M,
                  main=colnames(se.filt)[i],
                  las=1,
                  cex.axis=1.2,
                  cex.lab=1.5,
                  cex.main=1,
                  colramp = colorRampPalette(c("white",ifelse(C == "tumor", "blue","green")),space = "Lab"))
    abline(h=0, col="blue", lwd=2) ; lo <- lowess(M ~ A) ; lines(lo$x, lo$y, col="red", lwd=2)
  }
}

## MDS
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))
plotMDS(dge_luad.filt, col = c("red", "blue")[as.integer(dge_luad.filt$samples$group)], cex = 0.7)
legend("topright", c("Normal", "Tumor"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)

## Remove samples TCGA.64.5775.01A.01R.1628.07 and perhaps TCGA.49.AAR9.01A.21R.A41B.07
maskbad <- colnames(dge_luad) %in% c("TCGA.64.5775.01A.01R.1628.07", "TCGA.49.AAR9.01A.21R.A41B.07")
se.filt <- se.filt[, !maskbad]
dge_luad.filt <- dge_luad.filt[, !maskbad]

# faltaria fer per els altres 3 casos comprobats amb el GGplot

##############
## Session 2
#############

## Identifying batch effect: We'll use TCGA
tss <- substr(colnames(se.filt), 6, 7)
table(tss)

center <- substr(colnames(se.filt), 27, 28)
table(center)
#All samples on same center: we can delete

plate <- substr(colnames(se.filt), 22, 25)
table(plate)

portionanalyte <- substr(colnames(se.filt), 18, 20)
table(portionanalyte)

samplevial <- substr(colnames(se.filt), 14, 16)
table(samplevial)

gender <- unname(se.filt$gender)

race <- unname(se.filt$race)

histo <- unname(se.filt$histologic_diagnosis.1)

## Cross tabulation data

table(data.frame(TYPE=se.filt$type, TSS=tss))
# Too much zeros. Only six TSS have noral samples

table(data.frame(TYPE=se.filt$type, PLATE=plate))
# Too much zeros: only 5 plates with normal. 
# One of the plates (1758) only has normals!!!

table(data.frame(TYPE=se.filt$type, PORTIONALYTE=portionanalyte))
# All normal samples have extremly low values of portion analyte
# Not necessary: tumor and normal have separated values

table(data.frame(TYPE=se.filt$type, GENDER=gender))
#Seems pretty equilibred

table(data.frame(TYPE=se.filt$type, RACE=race))
#Mising too much dat

## Batch effects plot MDS

#Titles and objects lists
titles_batchs <- c("MDS plot: Tisue Source Site",
                   "MDS plot: Plate",
                   "MDS plot: Sex",
                   "MDS plot: Ethnicity",
                   "MDS plot: Histology")
objects_batchs <- list(tss, plate, gender,race, histo)
outcome <- paste(substr(colnames(se.filt), 9, 12), as.character(se.filt$type), sep="-")
for (a in 1:5){
  png(filename = paste("projct/img/", titles_batchs[a], sep = ""), width = 900, height = 800)
  par(mfrow=c(1,1))
  logCPM <- cpm(dge_luad.filt, log=TRUE, prior.count=3)
  d <- as.dist(1-cor(logCPM, method="spearman"))
  sampleClustering <- hclust(d)
  batch <- as.integer(factor(objects_batchs[[a]]))
  plotMDS(dge_luad.filt, labels=outcome, col=batch, main =titles_batchs[a])
  legend("bottomright", paste("Batch", sort(unique(batch)), levels(factor(a))),
        fill=sort(unique(batch)), inset=0.05, cex = 0.7)
  dev.off()
}

png(filename = "projct/img/MDS_histology_tumoronly.png", width = 900, height = 800)
se.tumor <- se.filt[, se$type == "tumor"]
par(mfrow=c(1,1))
tumor_dge <- DGEList(counts = assays(se.tumor)$counts, genes = as.data.frame(mcols(se.tumor)))
outcome <- paste(substr(colnames(se.tumor), 9, 12), as.character(se.tumor$type), sep="-")
logCPM <- cpm(tumor_dge, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(unname(se.tumor$histologic_diagnosis.1)))
plotMDS(tumor_dge, labels=outcome, col=batch, main ="MDS tumor by histology")
legend("topleft", paste("Batch", sort(unique(batch)), levels(factor(histo))),
       fill=sort(unique(batch)), inset=0.05, cex = 0.70)
dev.off()



########
##3. SVA
########

mod <- model.matrix(~ se.filt$type, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))
pv <- f.pvalue(assays(se.filt)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.01)
# We have 8076 differentially expressed genes

# Estimate surrogate variables
sv <- sva(assays(se.filt)$logCPM, mod, mod0)
sv$n
# 50 surrogate variables

#Adjusting surrogate variables
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(se.filt)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="fdr") < 0.01)
#9396 genes differentially expressed now