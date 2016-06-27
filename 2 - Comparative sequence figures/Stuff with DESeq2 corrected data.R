library(DESeq2); library(reshape2); library(ggplot2)
source("D:\\Dropbox\\R programming\\My functions\\namesplit.R")

df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)
tissue <- "mLN"
df.1 <- subset(df, tissue == tissue)
# df.1 <- subset(df.1, population == "Foxp3+RORgt+" | population == "Foxp3+")
df.c <- dcast(df.1, trvcdr3 ~ mouseNum + population, fun.aggregate = sum, value.var = "numReads")
tcrs <- df.c$trvcdr3
samples <- colnames(df.c)[-1]
samples <- namesplit.f(samples, "_")[,2]
samples <- gsub("\\+", "", samples)
tcr.m <- as.matrix(df.c[,-1])
rownames(tcr.m) <- tcrs
tcr.m <- (tcr.m) + 1 #transforms to avoid 0 count error in DESeq()
# tcr.m[1,] <- tcr.m[1,]+1
colData <- data.frame("population" = samples)
rownames(colData) <- colnames(tcr.m)
dds <- DESeqDataSetFromMatrix(countData = tcr.m,
                              colData = colData,
                              design = ~ population)
plot(assay(dds)[,1:3], cex = 0.5)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
log.norm <- normTransform(dds)
# rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
plot(assay(dds)[,1:2], cex = 0.5)
plot(assay(log.norm)[,1:2], cex = 0.5)
# plot(assay(rld)[,1:2], cex = 0.5)
plot(assay(vsd)[,1:2], cex = 0.5)

plotSparsity(dds)

plot(hclust(dist(t(assay(vsd)))))

library(pvclust)
x <- pvclust(assay(vsd), nboot = 50)
plot(x)
pvrect(x)
vegdist(t(assay(vsd)), method = "horn")
y <- decostand(assay(vsd), "range", 2)
vegdist(t(y), method = "horn")
pv.veg <- function(x) {return(vegdist(t(x), method = "horn"))}

z <- pvclust(y, method.dist = pv.veg, nboot = 50)
head(y)



df.1 <- subset(df, tissue == "mLN")

ls <- split(df.1, df.1$mouseNum)
ls.c <- lapply(ls, dcast, trvcdr3 ~ population, fun.aggregate = sum, value.var = "numReads")

create.dds.list <- function(x){
  tcr.m <- as.matrix(x[,-1])
  tcr.m <- (tcr.m) + 1
  samples <- colnames(x)[-1]
  colData <- data.frame("population" = samples)
  rownames(colData) <- colData$population
    dds <- DESeqDataSetFromMatrix(countData = tcr.m,
                                colData = colData,
                                design = ~ population)
  return(dds)
}
ls.dds <- lapply(ls.c, create.dds.list)
ls.vsd <- lapply(ls.dds, varianceStabilizingTransformation)

norm.dds.ls <- function(x) return(t(decostand(assay(x), "total", 2)))

ls.norm <- lapply(ls.dds, norm.dds.ls)
ls.mh <- lapply(ls.norm, vegdist, method = "horn")
ls.mh <- lapply(ls.mh, dist2df.f)
df.mh <- do.call("rbind", ls.mh)
df.mh$comp <- paste(df.mh$Var1, df.mh$Var2, sep = "_")
df.mh$value <- 1-df.mh$value
levels(factor(df.mh$comp))


se <- function(x) {
  s <- sd(x)/sqrt(length(x))
  m <- mean(x)
  data.frame("ymin" = m-s, "ymax" = m+s)  
}

ggplot(df.mh, aes(x = comp, y = value)) +
  stat_summary(fun.y = mean, geom = "bar") +
  stat_summary(fun.data = se, geom = "errorbar", width = 0.5, size = 1) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

TukeyHSD(aov(value ~ comp, data = df.mh))$comp
