library(dplyr)
library(vegan)
library(reshape2)
library(gplots)
library(RColorBrewer)


source("D:/Dropbox/R programming/My functions/namesplit.R")
source("D:/Dropbox/R programming/My functions/dist2df.R")

df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)

df.cdr3counts <- df %>% group_by(cdr3, mouseNum, population, tissue) %>% summarise(count = sum(numReads))
df.cdr3counts$length <- nchar(as.character(df.cdr3counts$cdr3))
df.cdr3length <- df.cdr3counts %>% group_by(length, mouseNum, population, tissue) %>% summarise(count = sum(count))
df.cdr3length <- df.cdr3length %>% group_by(mouseNum, population, tissue) %>% mutate(total = sum(count), frequency = count/total)

tis <- "Colon"

df.c <- subset(df.cdr3length, tissue == tis)
df.c <- dcast(df.c, mouseNum + population ~ length, value.var = "frequency", fun.aggregate = sum)
rownames(df.c) <- paste(df.c[,1],df.c[,2],sep = "_")

color.key <- c("CD44hi" = brewer.pal(4, "Set1")[1], 
               "Foxp3+" = brewer.pal(4, "Set1")[2],
               "RORgt+" = brewer.pal(4, "Set1")[3],
               "Foxp3+RORgt+" = brewer.pal(4, "Set1")[4])
col.color <- color.key[df.c$population]

hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(1000)

heatmap.2(as.matrix(vegdist(df.c[,c(-1,-2)], method = "horn")), 
          distfun = function(x) as.dist(x),
          col = hmcol,
          trace = "none",
          margins = c(15,15),
          ColSideColors = col.color)
legend("topright",
       legend = names(color.key),
       col = color.key,
       lty= 1, lwd = 5,
       border=FALSE, bty="n", y.intersp = 0.8, cex=0.8)
