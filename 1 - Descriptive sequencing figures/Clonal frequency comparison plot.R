# library(plyr)
library(dplyr)
library(vegan)
library(reshape2)
library(ggplot2)

##Functions
source("D:\\Dropbox\\R programming\\My functions\\namesplit.R")

df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)


df.1 <- df %>%
  group_by(mouseNum, tissue, population) %>%
  mutate(total = sum(numReads), frequency = numReads/total) %>%
  group_by(trvcdr3, tissue, population) %>%
  summarise(frequency = mean(frequency)) %>%
  group_by(tissue, population) %>%
  arrange(desc(frequency)) %>%
  mutate(rank = 1:length(frequency)) %>%
  filter(population == "Foxp3+RORgt+" | population == "RORgt+") %>%
  filter(tissue == "Colon")

df.order <- dcast(df.1, trvcdr3 ~ population, value.var = "rank")
df.order <- arrange(df.order, `Foxp3+RORgt+`, desc(`RORgt+`))
# y <- y[apply(y, 1, function(x) !any(is.na(x))),]
df.order$order <- 1:nrow(df.order)

df.1 <- merge(dcast(df.1, trvcdr3~population, value.var = "frequency"), df.order[,c(1,4)], by = "trvcdr3")
df.1[is.na(df.1)] <- 0
df.1 <- df.1 %>% select(`Foxp3+RORgt+`, `RORgt+`, order)
df.1 <- df.1[apply(df.1[,1:2], 1, function(x) any(x > 0.01)),]
df.1$order[order(df.1$order, decreasing = F)]
df.1 <- df.1[order(df.1$order, decreasing = F),]
df.1$order <- 1:nrow(df.1)
df.1$`RORgt+` <- -df.1$`RORgt+`
df.1 <- melt(df.1, id.vars = c("order"), value.name = "frequency", variable.name = "population")


# ggplot(df.1, aes(x = order, y = frequency)) +
#   geom_bar(aes(fill = population), stat = "identity", position = "dodge") +
#   theme_bw() +
#   facet_grid(population ~ ., scales = "free")


ggplot() +
  geom_bar(data=subset(df.1, population == "Foxp3+RORgt+"), aes(x=order, y = frequency, fill=population), stat="identity") +
  geom_bar(data=subset(df.1, population == "RORgt+"), aes(x=order, y = frequency, fill=population), stat="identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw()

df.2 <- dcast(df.1, order ~ population, value.var = "frequency", fun.aggregate = sum)
df.2$`RORgt+` <- -df.2$`RORgt+`
df.2$`RORgt+`[df.2$`RORgt+` == 0] <- min(df.2$`RORgt+`[df.2$`RORgt+` > 0])
df.2$`Foxp3+RORgt+`[df.2$`Foxp3+RORgt+` == 0] <- min(df.2$`Foxp3+RORgt+`[df.2$`Foxp3+RORgt+` > 0])
df.2$logR <- log(df.2$`Foxp3+RORgt+`/df.2$`RORgt+`)

ggplot(df.2, aes(x = order, y = logR)) + geom_bar(stat="identity") + theme_bw()
