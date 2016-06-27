library("ggplot2")
library("reshape2")
library("Hmisc")
library("dplyr")

source("D:\\Dropbox\\R programming\\My functions\\namesplit.R")

df <- read.csv("160504 - Compiled CX3CR1 transfer - All - CXCR3all.csv", header = T)
names(df) <- c("sample", "foxp3.freq", "rorgt.freq", "experiment", "phenotype", "tissue")
df <- df[1:38,1:6]
df <- melt(df, id.vars = c("sample","experiment", "phenotype", "tissue"))
df$phenotype <- factor(df$phenotype, levels = c("WT", "CKO"))


se <- function(x) {
  s <- sd(x)/sqrt(length(x))
  m <- mean(x)
  data.frame("ymin" = m-s, "ymax" = m+s)  
}

# df[,2] <- factor(df[,2], c("WT", "KO"))

dot.size <- 1*max(df$value)/100
cb <- "Set1"

ggplot(df, aes(x = phenotype, y = value)) + 
  theme_classic() +
  stat_summary(fun.y = mean, geom = "bar") +
  stat_summary(fun.data = se, geom = "errorbar", width = 0.5, size = 1) +
  geom_point() +
  facet_grid(variable ~ tissue, scales = "free")
# 
# ggplot(df, aes(x = tissue, y = frequency)) + 
#   theme_classic() +
#   scale_fill_brewer(palette=cb) +
#   stat_summary(aes(fill = phenotype), fun.y = mean, geom = "bar", position = position_dodge(width = 0.9)) +
#   geom_dotplot(aes(fill = phenotype), position = "dodge", binaxis = "y", binwidth = dot.size, stackdir = "center") +
#   stat_summary(aes(fill = phenotype), fun.data = se, geom = "errorbar", width = 0.5, size = 1, position = position_dodge(width = 0.9)) +
#   facet_wrap(~ gate, scales = "free")

test <- df %>% group_by(tissue, variable) %>% do(w = wilcox.test(value ~ phenotype, data = .)) %>% summarise(tissue, variable, w$p.value)
write.csv(test, "wilcox.csv")
