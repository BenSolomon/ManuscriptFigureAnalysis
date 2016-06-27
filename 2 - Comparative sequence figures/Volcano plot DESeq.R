library(DESeq2); library(reshape2); library(ggplot2)

source("D:\\Dropbox\\R programming\\My functions\\namesplit.R")

df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)
tis <- "Colon"
a <- 0.05

df.1 <- subset(df, tissue == tis)
df.1 <- subset(df.1, population == "Foxp3+RORgt+" | population == "Foxp3+")
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
levels(dds$population)
dds$population <- relevel(dds$population, "Foxp3")
levels(dds$population)
dds <- DESeq(dds)
res.1 <- results(dds, alpha = a)
# res
res1.df <- as.data.frame(res.1)
res1.df$tcr <- rownames(res1.df)



df.1 <- subset(df, tissue == tissue)
df.1 <- subset(df.1, population == "Foxp3+RORgt+" | population == "RORgt+")
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
levels(dds$population)
dds$population <- relevel(dds$population, "RORgt")
levels(dds$population)
dds <- DESeq(dds)
res.2 <- results(dds, alpha = a)
res2.df <- as.data.frame(res.2)
res2.df$tcr <- rownames(res2.df)

  

res.1.s <- subset(res.1, padj <= a & log2FoldChange >= 1)
res.2.s <- subset(res.2, padj <= a & log2FoldChange >= 1)

intsct <- intersect(rownames(res.1.s), rownames(res.2.s))

res1.df$hit <- res1.df$tcr %in% intsct 
res2.df$hit <- res2.df$tcr %in% intsct

ct2 <- "TRAV14-3*01_FCAASAIWNTGYQNFYFG"
res1.df$ct2 <- res1.df$tcr == ct2
res2.df$ct2 <- res2.df$tcr == ct2

res1.df$ct2hit <- apply(res1.df[,c("hit","ct2")], 1, all)
res2.df$ct2hit <- apply(res2.df[,c("hit","ct2")], 1, all)


ggplot(res1.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 1) +
  geom_point(data = res1.df[res1.df$ct2hit,], color = "white", size = 2, shape = 4, stroke = 4) +
  geom_point(data = res1.df[res1.df$ct2hit,], color = "blue", size = 2, shape = 4, stroke = 3) +
  geom_point(data = res1.df[res1.df$hit,], color = "red") +
  theme_bw() +
  ggtitle(label = paste0("DP vs. Foxp3 -",tis))
ggplot(res2.df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 1) +
  geom_point(data = res2.df[res2.df$ct2hit,], color = "white", size = 2, shape = 4, stroke = 4) +
  geom_point(data = res2.df[res2.df$ct2hit,], color = "blue", size = 2, shape = 4, stroke = 3) +
  geom_point(data = res2.df[res2.df$hit,], color = "red") +
  theme_bw() +
  ggtitle(label = paste0("DP vs. RORgt - ",tis))

FRvFtcrs <- res1.df[res1.df$hit,c("log2FoldChange", "padj", "tcr")]
FRvRtcrs <- res2.df[res2.df$hit,c("log2FoldChange", "padj", "tcr")]


TCRsummary <- merge(FRvFtcrs, FRvRtcrs, by = "tcr")
names(TCRsummary) <- c("tcr", "FRvF.log2FoldChange", "FRvF.padj", "FRvR.log2FoldChange", "FRvR.padj")
write.csv(x = TCRsummary, paste("DESeq2 TCR list - ", tissue, ".csv", sep = ""), row.names = F)
