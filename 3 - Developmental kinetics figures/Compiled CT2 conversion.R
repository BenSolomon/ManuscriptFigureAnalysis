library("ggplot2")
library("reshape2")
library("Hmisc")
library("dplyr")

se <- function(x) {
  s <- sd(x)/sqrt(length(x))
  m <- mean(x)
  data.frame("ymin" = m-s, "ymax" = m+s)  
}

df <- read.csv("D:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/3 - Developmental kinetics figures/CT2 transfer kinetics - All - CXCR3all.csv", header = T)
df <- df[,-12]
df <- df[1:70,]
df <- melt(df, id.vars = c("X", "RAG", "Tissue", "Day"), variable.name = "gate", value.name = "frequency")
names(df) <- c("sample", "rag", "tissue", "day", "gate", "frequency")
df$day <- gsub("d", "", df$day)
df$day <- as.numeric(df$day)
df$day[df$day == 23] <- 21 #Combines d21 and d23
df <- subset(df, day == 21)
df <- subset(df, gate == "RORinFOX" | gate == "FOXall")
df.tally <- df %>% select(sample, tissue, gate, frequency)
names(df.tally) <- c("sample", "tissue", "variable", "value")

df <- read.csv("D:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/4 - APC figures/Transfer experiments/CX3CR1 transfer/160504 - Compiled CX3CR1 transfer - All - CXCR3all.csv", header = T)
names(df) <- c("sample", "foxp3.freq", "rorgt.freq", "experiment", "phenotype", "tissue")
df <- df[1:38,1:6]
df <- melt(df, id.vars = c("sample","experiment", "phenotype", "tissue"))
df$phenotype <- factor(df$phenotype, levels = c("WT", "CKO"))
df <- subset(df, phenotype == "WT")
df.tally <- rbind(df.tally, df %>% select(sample, tissue, variable, value))

df <- read.csv("D:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/4 - APC figures/Transfer experiments/Notch2 transfer/Compiled Notch2 transfer - All - CXCR3all.csv", header = T)
names(df) <- c("sample", "foxp3.freq", "rorgt.freq", "experiment", "phenotype", "tissue")
df <- df[1:26,1:6]
df <- melt(df, id.vars = c("sample","experiment", "phenotype", "tissue"))
df$phenotype <- factor(df$phenotype, levels = c("WT", "KO"))
df <- subset(df, phenotype == "WT")
df.tally <- rbind(df.tally, df %>% select(sample, tissue, variable, value))

df <- read.csv("D:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/4 - APC figures/Transfer experiments/IRF4 transfer/Compiled IRF4 transfer - All - CXCR3all.csv", header = T)
names(df) <- c("sample", "foxp3.freq", "rorgt.freq", "phenotype", "experiment", "tissue")
df <- df[1:38,1:6]
df <- melt(df, id.vars = c("sample","experiment", "phenotype", "tissue"))
df$phenotype <- factor(df$phenotype, levels = c("WT", "KO"))
df <- subset(df, phenotype == "WT")
df.tally <- rbind(df.tally, df %>% select(sample, tissue, variable, value))

df.tally$variable[df.tally$variable == "FOXall"] <- "foxp3.freq"
df.tally$variable[df.tally$variable == "RORinFOX"] <- "rorgt.freq"

df.tally$variable <- droplevels(df.tally$variable)
table(df.tally$variable)

ggplot(df.tally, aes(x = tissue, y = value)) + 
  theme_classic() +
  stat_summary(fun.y = mean, geom = "bar") +
  geom_point() +
  stat_summary(fun.data = se, geom = "errorbar", width = 0.25, size = 1) +
  facet_wrap(~ variable, scales = "free")
