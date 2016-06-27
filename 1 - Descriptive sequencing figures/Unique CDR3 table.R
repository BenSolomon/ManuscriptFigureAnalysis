library(reshape2);library(dplyr)
df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)
summary <- df %>% group_by(population, tissue) %>% summarise(count = length(levels(factor(trvcdr3))))
summary <- dcast(summary, population ~ tissue, fun.aggregate = sum, value.var = "count")
row.names(summary) <- summary$population
summary <- summary[,-1]
write.csv(summary, "UniqueCDR3table.csv")
# cbind(summary, apply(summary, 1, "sum"))
# rbind(summary, apply(summary, 2, "sum"))
