library(reshape2); library(dplyr); library(vegan)
source("D:\\Dropbox\\R programming\\My functions\\namesplit.R")



df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)
tis <- "Colon"
df.s <- subset(df, tissue == tis)
df.c <- dcast(df.s, trvcdr3 ~ mouseNum + population, fun.aggregate = sum, value.var = "numReads")
cdr3 <- df.c[,1]
df.c <- decostand(df.c[,-1], "total", 2)
rownames(df.c) <- cdr3
df.t <- t(df.c)
df.t <- as.data.frame(df.t)
df.t$population <- namesplit.f(rownames(df.t), "_")[,2]
# summary <- df.t %>% group_by(population) %>% summarise_each(funs = c("mean"))
# f.t$group <- paste(namesplit.f(rownames(df.t), "_")[,2], namesplit.f(rownames(df.t), "_")[,3], sep = "_")
split.f <- function(x){
  y <- split(x, x[,length(x)])
  y <- sapply(y, function(x) apply(x[,-length(x)], 2, mean))
  return(y)
}
summary <- split.f(df.t)
summary <- as.data.frame(summary)
summary$tcr <- row.names(summary)


select.tcrs <- read.csv("D:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/2 - Comparative sequence figures/TCR tables/DESeq2 TCR list - Colon.csv", header = T)
# varimp.tcrs <- read.csv("D:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/2 - Comparative sequence figures/Random forest/Colon TCR imp.csv", header = T)


summary.s <- summary[summary$tcr %in% select.tcrs$tcr,] 
compiled.df <- merge(summary.s, select.tcrs, by = "tcr")

write.csv(compiled.df, "TCR table - Colon.csv", row.names = F)
# 
# varimp.s <- varimp.tcrs[varimp.tcrs$tcr %in% select.tcrs$tcr,] 
# compiled.df.2 <- merge(compiled.df, varimp.s, by = "tcr")

df <- read.csv("DESeq2 TCR list - Colon.csv", header = T)
df$cdr3 <- namesplit.f(df$tcr, "_")$V2

key <- read.csv("namedTCRkey.csv", header = T)
key$cdr3 <- paste("FC",key$cdr3,"FG",sep = "")

merge(df, key, by = "cdr3")
