library(SummarizedExperiment)
library(edgeR)
library(ggplot2)
library(plyr)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
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






SE <- readRDS("seLUAD.rds")
dge <- DGEList(counts = assays(SE)$counts, group = SE$gender, genes = as.data.frame(rowData(SE)))
logCPM <- cpm(dge, log = TRUE, prior.count = 0.25)
ord <- order(dge$sample$lib.size)

# BARPLOTS
i <- 0
ordsliced <- list(0)
gender <- list(0)
race <- list(0)
typ <- list(0)
hstory_malignancy <- list(0)
for(num in 1:6){
  ordsliced[[num]] <- c(ord[(i*100 +1):(num*100)])
  i <- num
  g_nas <- length(colData(SE)$gender[ordsliced[[num]]][colData(SE)$gender[ordsliced[[num]]] == "NA"])
  gender[[num]] <- c(length(colData(SE)$gender[ordsliced[[num]]][colData(SE)$gender[ordsliced[[num]]] == "MALE"])-g_nas,
                     length(colData(SE)$gender[ordsliced[[num]]][colData(SE)$gender[ordsliced[[num]]] == "FEMALE"])-g_nas, g_nas)
  
  race_nas <- length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "NA"])
  race_tnas <- length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "[Not Available]"]) +
              length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "[Not Evaluated]"]) +
              length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "[Unknown]"]) -2*race_nas
  race[[num]] <- c(length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "AMERICAN INDIAN OR ALASKAN NATIVE"]) - race_nas,
                   length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "ASIAN"]) - race_nas,
                   length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "BLACK OR AFRICAN AMERICAN"]) - race_nas,
                   length(colData(SE)$race[ordsliced[[num]]][colData(SE)$race[ordsliced[[num]]] == "WHITE"]) - race_nas,
                   race_tnas)
  hist_nas <- length(colData(SE)$history_other_malignancy[ordsliced[[num]]][colData(SE)$history_other_malignancy[ordsliced[[num]]] == "NA"])
  hstory_malignancy[[num]] <- c(length(colData(SE)$history_other_malignancy[ordsliced[[num]]][colData(SE)$history_other_malignancy[ordsliced[[num]]] == "Yes"]) - hist_nas +
                                length(colData(SE)$history_other_malignancy[ordsliced[[num]]][colData(SE)$history_other_malignancy[ordsliced[[num]]] == "Yes, History of Prior Malignancy"]) - hist_nas + 
                                length(colData(SE)$history_other_malignancy[ordsliced[[num]]][colData(SE)$history_other_malignancy[ordsliced[[num]]] == "Yes, History of Synchronous and or Bilateral Malignancy"]) - hist_nas,
                                length(colData(SE)$history_other_malignancy[ordsliced[[num]]][colData(SE)$history_other_malignancy[ordsliced[[num]]] == "No"]) - hist_nas,
                                hist_nas)
  typ_na <- length(colData(SE)$type[ordsliced[[num]]][colData(SE)$type[ordsliced[[num]]] == "NA"])
  typ[[num]] <- c(length(colData(SE)$type[ordsliced[[num]]][colData(SE)$type[ordsliced[[num]]] == "normal"]) - typ_na,
                  length(colData(SE)$type[ordsliced[[num]]][colData(SE)$type[ordsliced[[num]]] == "tumor"]) - typ_na,
                  typ_na)
  }

# Gender BARPLOT
g_df <- data.frame(supp=c("MALE", "FEMALE", "NA"),
                 depth=rep(c('1vlow', '2low', '3midlow', '4mid', '5hih', '6vhigh'),each = 3),
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
                   depth=rep(c('1vlow', '2low', '3midlow', '4mid', '5hih', '6vhigh'),each = 5),
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
                   depth=rep(c('1vlow', '2low', '3midlow', '4mid', '5hih', '6vhigh'),each = 3),
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
                   depth=rep(c('1vlow', '2low', '3midlow', '4mid', '5hih', '6vhigh'),each = 3),
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



meanUnloggedExp <- rowMeans(assays(lclse)$count)
sdUnloggedExp <- apply(assays(lclse)$counts, 1, sd)
meanLoggedExp <- rowMeans(logCPM)
sdLoggedExp <- apply(logCPM, 1, sd)


#variance correction // negative binomial distribution when the variance is much more of the mean
par(mfrow=c(1, 2))
plot(meanUnloggedExp, sdUnloggedExp, pch=".", cex=4, xlab="Unlogged mean expression", ylab="SD")
lines(lowess(meanUnloggedExp, sdUnloggedExp), lwd=2, col="red")
plot(meanLoggedExp, sdLoggedExp, pch=".", cex=4, xlab="Logged mean expression", ylab="SD")
lines(lowess(meanLoggedExp, sdLoggedExp, f=0.25), lwd=2, col="red")

#Fold-change

maleExp <- rowMeans(logCPM[, lclse$gender == "MALE"])
femaleExp <- rowMeans(logCPM[, lclse$gender == "FEMALE"])

par(mfrow = c(1, 2))
plot(maleExp, femaleExp, xlab = "Male", ylab = "Female", pch = ".", cex = 4, las = 1, xlim = c(-5, 10000), ylim = c(-5, 10000))
plot((femaleExp + maleExp)/2, femaleExp - maleExp, pch = ".", cex = 4, las = 1)

log2fc <- femaleExp - maleExp
ranking <- order(abs(log2fc), decreasing = TRUE)
head(data.frame(Log2FC = round(log2fc[ranking], digits = 3), FC = round(2^log2fc[ranking],
                digits = 3), `1/FC` = round(2^(-log2fc[ranking]), digits = 3), row.names = rowData(lclse)$symbol[ranking], check.names = FALSE), n = 10)
