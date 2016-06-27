library(ggplot2);library(vegan);library(reshape2);library(pvclust);library(dplyr)
source("D:/Dropbox/R programming/My functions/namesplit.R")

# df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)
# df.c <- dcast(df, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
# df.c <- decostand(df.c[,-1], "total", 2)
# samples <- names(df.c)
# df.t <- t(df.c)
# pv.veg <- function(x) {return(vegdist(t(x), method = "horn"))}
# y <- pvclust(df.c, method.dist = pv.veg, nboot = 100)
# y <- pvclust(df.c, nboot = 5)
# plot(y)
# pvrect(y)
# 
# # pv <- function(x){
# #   z <- df.c[,which(grepl("AR922", names(df.c)))]
# #   vegdist(t(z), method = "horn")
# # }
# 
# 
# 
# 
# df$split <- paste(df$mouseNum, df$tissue, df$tissue, sep = "_")
# ls <- split(df, df$split)
# ls.c <- lapply(ls, dcast, trvcdr3 ~ split + population, fun.aggregate = sum, value.var = "numReads")
# ls.c <- lapply(ls.c, function(x) decostand(x[,-1], "total", 2))
# ls.c <- lapply(ls.c, t)
# ls.mh <- lapply(ls.c, vegdist, method = "horn")
# ls.mh <- lapply(ls.mh, dist2df.f)
# df.mh <- do.call("rbind", ls.mh)
# df.mh$mouseNum <- namesplit.f(df.mh$Var1, "_")[,1]
# df.mh$tissue <- namesplit.f(df.mh$Var1, "_")[,2]
# df.mh$comp <- paste(namesplit.f(df.mh$Var1, "_")[,4], namesplit.f(df.mh$Var2, "_")[,4], sep = "_")

# levels(factor(df$mouseNum))
# x <- subset(df, mouseNum == "AR962/964")
# df.c <- dcast(x, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
# df.c <- decostand(df.c[,-1], "total", 2)
# plot(hclust(vegdist(t(df.c), method = "horn")))
# plot(hclust(dist(as.matrix(vegdist(t(df.c), method = "horn")))))


# for (i in 1:length(levels(factor(df$mouseNum)))){
#   x <- subset(df, mouseNum == levels(factor(df$mouseNum))[i])
#   df.c <- dcast(x, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
#   df.c <- decostand(df.c[,-1], "total", 2)
#   print(dim(as.matrix(vegdist(t(df.c), method = "horn"))))
# }

#Only mice with all samples
df.s <- subset(df, mouseNum == levels(factor(df$mouseNum))[c(3,5,6)])
for (i in 1:length(levels(factor(df.s$mouseNum)))){
  x <- subset(df.s, mouseNum == levels(factor(df.s$mouseNum))[i])
  df.c <- dcast(x, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
  df.c <- decostand(df.c[,-1], "total", 2)
  print(dim(as.matrix(vegdist(t(df.c), method = "horn"))))
}

df.s$split <- factor(df.s$mouseNum)
ls <- split(df.s, factor(df.s$mouseNum))
ls.c <- lapply(ls, dcast, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
ls.c <- lapply(ls.c, function(x) t(x[,-1]))
ls.mh <- lapply(ls.c, vegdist, method = "horn")
ls.mh <- lapply(ls.mh, as.matrix)
cnames <- paste(namesplit.f(colnames(ls.mh[[1]]), "_")[,2], namesplit.f(colnames(ls.mh[[1]]), "_")[,3], sep = "_")
ls.mh <- lapply(ls.mh, function(x) {colnames(x) <- cnames; return(x)})
df.mh <- do.call("rbind", ls.mh)

y <- pvclust(df.mh, n = 500)
plot(y)
pvrect(y)


pcoa(df.mh)
PC <- pcoa(dist(t(df.mh))) 
PC.df <- as.data.frame(PC$vectors) 
PC.df$sample <- rownames(PC.df) 
ggplot(data = PC.df, aes(x = Axis.1, y = Axis.2)) +
  geom_point() +
  geom_text(aes(label = sample))


#Uses all mice, w/ NAs for missing samples
df$split <- factor(df$mouseNum)
ls <- split(df, factor(df$mouseNum))
ls.c <- lapply(ls, dcast, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
ls.c <- lapply(ls.c, function(x) t(x[,-1]))
ls.mh <- lapply(ls.c, vegdist, method = "horn")
ls.mh <- lapply(ls.mh, as.data.frame(as.matrix))
list.rename <- function(x){
  colnames(x) <- paste(namesplit.f(colnames(x), "_")[,2], namesplit.f(colnames(x), "_")[,3], sep = "_")
  return(x)
}
all.names <- colnames(list.rename(ls.mh[[3]]))
ls.mh <- lapply(ls.mh, list.rename)

list.expand <- function(x){
  missing.names <- all.names[!all.names %in% names(x)]
  missing.mat <- matrix(NA, nrow = nrow(x), ncol = length(missing.names))
  colnames(missing.mat) <- missing.names
  full.df <- cbind(x, as.data.frame(missing.mat))
  ordered.df <- select_(full.df, .dots = all.names)
  return(ordered.df)
}

ls.mh <- lapply(ls.mh, list.expand)


# cnames <- paste(namesplit.f(colnames(x), "_")[,2], namesplit.f(colnames(x), "_")[,3], sep = "_")
# ls.mh <- lapply(ls.mh, function(x) {colnames(x) <- cnames; return(x)})
df.mh <- do.call("rbind", ls.mh)
y <- pvclust(df.mh, n = 1000)
# plot(y, hang = -1)
plot(y)
pvrect(y)

pcoa(df.mh)
PC <- pcoa(dist(t(df.mh))) #recall dm.trgv is a distance matrix, which can be used by pcoa
# #Class "pcoa"
# #Function in {ape}
PC.df <- as.data.frame(PC$vectors) #extract PC axis info
PC.df$sample <- rownames(PC.df) #create variable of trgv names
ggplot(data = PC.df, aes(x = Axis.1, y = Axis.2)) +
  geom_point() +
  geom_text(aes(label = sample))
