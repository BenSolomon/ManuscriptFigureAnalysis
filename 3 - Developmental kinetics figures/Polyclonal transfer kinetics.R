library("ggplot2")
library("reshape2")
library("Hmisc")
library("plyr")
library("dplyr")

se <- function(x) {
  s <- sd(x)/sqrt(length(x))
  m <- mean(x)
  data.frame("ymin" = m-s, "ymax" = m+s)  
}

df <- read.csv("Polyclonal transfer kinetics - All - CXCR3all.csv", header = T)
df <- df[,-12]
df <- df[1:23,]
df <- melt(df, id.vars = c("X", "Experiment", "Tissue", "Day"), variable.name = "gate", value.name = "frequency")
names(df) <- c("sample", "experiment", "tissue", "day", "gate", "frequency")
df$day <- gsub("d", "", df$day)
df$day <- as.numeric(df$day)
# df$day[df$day == 8] <- 7 #Combines d7 and d8
# df$day[df$day == 23] <- 21 #Combines d21 and d23
# df <- subset(df, day !=16) #Removes d16, since it only has data from 1 out of 3 replicates

df$day[df$day == 9] <- 6 #Combines d6 and d9
# df$day[df$day == 23] <- 21 #Combines d21 and d23
df <- subset(df, day !=16) #Removes d16, since it only has data from 1 out of 3 replicates

# #Plotting all gates, via facet
# ggplot(df, aes(x = day, y = frequency)) + 
#   theme_classic() +
#   geom_point() +
#   stat_summary(fun.y = mean, geom = "line") +
#   # stat_summary(fun.data = se, geom = "errorbar", width = 0.25, size = 1)  +
#   facet_grid(.~gate)
# 
# summary <-df %>% group_by(gate, day) %>% summarise(mean(frequency))

#Plotting Foxp3 vs RORgt
df.1 <- subset(df, gate %in% c("FOXonly", "FOXandROR", "RORonly"))
ggplot(df.1, aes(x = day, y = frequency)) + 
  theme_classic() +
  # geom_point() +
  stat_summary(fun.data = se, geom = "ribbon", aes(group = gate)) +
    stat_summary(fun.y = mean, geom = "line", aes(color = gate)) +
  facet_wrap(~ tissue, scales = "free")
   
#   
#   
# #Plotting various possible Foxp3 subsets
# df.2 <- subset(df, gate %in% c("RORinFOX", "FOXall", "FOXonly"))
# ggplot(df.2, aes(x = day, y = frequency)) + 
#   theme_classic() +
#   # geom_point() +
#   stat_summary(fun.y = mean, geom = "line", aes(color = gate)) +
#   stat_summary(fun.data = se, geom = "ribbon", alpha = 0.1, aes(group = gate)) +
#   facet_wrap(~ tissue)
