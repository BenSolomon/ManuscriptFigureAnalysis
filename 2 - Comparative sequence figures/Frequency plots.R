library(ggplot2);library(rehsape2);library(vegan)

#1 mouse
df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)
df.s <- subset(df, tissue == "Colon" & mouseNum == "AR5/6")
df.c <- dcast(df, trvcdr3 ~ population, fun.aggregate = sum, value.var = "numReads")
df.c <- df.c[,-1] + 1
df.c <- decostand(df.c, "total", 2)
plot(log10(df.c$`RORgt+`), log10(df.c$`Foxp3+RORgt+`))
cor.test(df.c$`Foxp3+`, df.c$`Foxp3+RORgt+`, method = "pearson")

df.c.2 <- log2(df.c)



#All mouse average
tissue <- "mLN"
df.s <- subset(df, tissue == tissue)
df.c <- dcast(df.s, trvcdr3 ~ mouseNum + population, fun.aggregate = sum, value.var = "numReads")
cdr3 <- df.c[,1]
df.c <- decostand(df.c[,-1], "total", 2)
df.t <- t(df.c)
df.t <- as.data.frame(df.t)
df.t$population <- namesplit.f(rownames(df.t), "_")[,2]
summary <- df.t %>% group_by(population) %>% summarise_each(funs = c("mean"))
# se.f <- function(x) sd(x)/sqrt(length(x))

# summary.se <- df.t %>% group_by(population) %>% summarise_each(funs("se.f"))
populations <- summary$population
summary <- t(summary[,-1])
colnames(summary) <- populations
summary <- as.data.frame(summary)
# qplot(x = log10(summary$`RORgt+`), y = log10(summary$`Foxp3+RORgt+`)) + theme_bw()
qplot(x = summary$`Foxp3+`, y = summary$`Foxp3+RORgt+`) + 
  scale_x_log10(breaks = c(1e-4,1e-3,1e-2)) + 
  scale_y_log10(breaks = c(1e-4,1e-3,1e-2)) + 
  theme_bw()
qplot(x = summary$`RORgt+`, y = summary$`Foxp3+RORgt+`) + 
  scale_x_log10(breaks = c(1e-4,1e-3,1e-2)) + 
  scale_y_log10(breaks = c(1e-4,1e-3,1e-2)) + 
  theme_bw()

write.csv(unlist(cor.test(summary$`Foxp3+`, summary$`Foxp3+RORgt+`, method = "pearson")), file = paste(tissue, "- FvsFR.csv"))
write.csv(unlist(cor.test(summary$`RORgt+`, summary$`Foxp3+RORgt+`, method = "pearson")), file = paste(tissue, "- RvsFR.csv"))




# #Simulation test

# tally <- vector()
# for (i in 1:10000){
#   x <- df.c$`Foxp3+RORgt+`
#   y.f <- apply(df.c, 1, function(x) x[sample(c(2,4), 1)])
#   y.r <- apply(df.c, 1, function(x) x[sample(c(2,4), 1)])
#   cor.f <- cor(x, y.f)
#   cor.r <- cor(x, y.r)
#   diff <- cor.r - cor.f
#   tally <- c(tally, diff)
# }
# tally
# hist(tally)
# 
# 
# actual.diff <- cor(df.c$`RORgt+`, df.c$`Foxp3+RORgt+`) - cor(df.c$`Foxp3+`, df.c$`Foxp3+RORgt+`)
# prop.table(table(abs(tally) >= actual.diff))
