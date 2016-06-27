library("ggplot2")
library("reshape2")
library("Hmisc")
library("dplyr")

se <- function(x) {
  s <- sd(x)/sqrt(length(x))
  m <- mean(x)
  data.frame("ymin" = m-s, "ymax" = m+s)  
}

#CT2 data
df <- read.csv("D:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/3 - Developmental kinetics figures/CT2 transfer kinetics - All - CXCR3all.csv", header = T)
df <- df[,-12]
df <- df[1:70,]
df <- melt(df, id.vars = c("X", "RAG", "Tissue", "Day"), variable.name = "gate", value.name = "frequency")
names(df) <- c("sample", "rag", "tissue", "day", "gate", "frequency")
df$day <- gsub("d", "", df$day)
df$day <- as.numeric(df$day)
df$day[df$day == 23] <- 21 #Combines d21 and d23
df <- subset(df, day == 21 | day == 16)
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
df.tally$variable <- as.character(df.tally$variable)
df.tally$variable[df.tally$variable == "FOXall"] <- "foxp3.freq"
df.tally$variable[df.tally$variable == "RORinFOX"] <- "rorgt.freq"
# df.tally$variable <- droplevels(df.tally$variable)
df.ct2 <- df.tally
df.ct2 <- dcast(df.ct2, sample + tissue ~ variable, fun.aggregate = sum)
df.ct2$DP <- ((df.ct2$foxp3.freq/100)*(df.ct2$rorgt.freq/100))*100
df.ct2 <- melt(df.ct2, id.vars = c("sample", "tissue"))
df.ct2 <- df.ct2 %>% filter(tissue == "cdLN") %>% select(sample, variable, value)
names(df.ct2) <- c("sample", "population", "frequency")
df.ct2$tcr <- "CT2"

#Polyclonal data
df <- read.csv("d:/Dropbox/Hsieh Lab/Data Analysis/Foxp3 RORgt Figures/3 - Developmental kinetics figures/Polyclonal transfer kinetics - All - CXCR3all.csv", header = T)
df <- df[,-12]
df <- df[1:23,]
df <- melt(df, id.vars = c("X", "Experiment", "Tissue", "Day"), variable.name = "gate", value.name = "frequency")
names(df) <- c("sample", "experiment", "tissue", "day", "gate", "frequency")
df$day <- gsub("d", "", df$day)
df$day <- as.numeric(df$day)
df$day[df$day == 9] <- 6 #Combines d6 and d9
df <- subset(df, day !=16) #Removes d16, since it only has data from 1 out of 3 replicates
df <- subset(df, (day == 20 | day == 15 | day == 16) & (gate == "RORinFOX" | gate == "FOXall" | gate == "FOXandROR"))
# df <- subset(df, day == 20 & (gate == "RORinFOX" | gate == "FOXall" | gate == "FOXandROR"))

df.poly <- df
df.poly <- df.poly %>% filter(tissue == "cdLN") %>% select(sample, gate, frequency)
names(df.poly) <- c("sample", "population", "frequency")
df.poly$population <- as.character(df.poly$population)
df.poly$population[df.poly$population == "FOXall"] <- "foxp3.freq"
df.poly$population[df.poly$population == "RORinFOX"] <- "rorgt.freq"
df.poly$population[df.poly$population == "FOXandROR"] <- "DP"
df.poly$tcr <- "polyclonal"

#CT6 data
df <- read.csv("CT6 conversion.csv", header = T)
df <- df[1:5, -5]
df$rorgt.freq <- apply(df[,c(-1,-2)], 1, max, na.rm = T)
df$foxp3.freq <- df[,2]
df$foxp3.freq[is.na(df$foxp3.freq)] <- mean(df$foxp3.freq, na.rm = T)
df <- df %>% select(X, foxp3.freq, rorgt.freq)
df$DP <- ((df$foxp3.freq/100)*(df$rorgt.freq/100))*100
df <- melt(df, id.vars = "X")
names(df) <- c("sample", "population", "frequency")
df$tcr <- "CT6"

#Combine
df <- rbind(df, df.ct2, df.poly)
df$tcr <- factor(df$tcr, levels = c("CT2", "CT6", "polyclonal"))
df$population <- factor(df$population, levels = c("rorgt.freq", "DP", "foxp3.freq"))

se <- function(x) {
  s <- sd(x)/sqrt(length(x))
  m <- mean(x)
  data.frame("ymin" = m-s, "ymax" = m+s)  
}

pal <- brewer.pal(9, "Greys")[c(1,5,9)]

ggplot(df, aes(x = tcr, y = frequency)) + 
  theme_classic() +
  stat_summary(fun.y = mean, geom = "bar", color = "black", aes(fill = tcr)) +
  scale_fill_manual(values = pal) +
  # geom_point() +
  stat_summary(fun.data = se, geom = "errorbar", width = 0.25, size = 1) +
  facet_wrap(~population, scales = "free")
