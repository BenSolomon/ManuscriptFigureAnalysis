library(limma); library(reshape2); library(dplyr); library(vegan)
hsb2 <- read.csv("http://www.ats.ucla.edu/stat/data/hsb2.csv")
attach(hsb2)

hw <- (write >= 60)
hm <- (math >= 60)
hr <- (read >= 60)
c3 <- cbind(hw, hm, hr)
a <- vennCounts(c3)
vennDiagram(a)

df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)
mice <- levels(factor(df$mouseNum))
# df.s <- subset(df, tissue == "mLN" & mouseNum == mice[5])
df.s <- subset(df, tissue == "mLN")
# ls <- split(df.s, df.s$population)
# ls <- lapply()

df.s %>% group_by(population) %>% summarise(length(unique(trvcdr3)))


df.c <- dcast(df.s, trvcdr3 ~ population, fun.aggregate = sum, value.var = "numReads")
df.c <- df.c[,c(-1,-2)] > 0
df.c <- df.c[apply(df.c,1, any), ]
vc <- vennCounts(df.c)
# vennDiagram(vc)
vc.r <- vc
vc.r[,4] <- round(100*vc.r[,4]/sum(vc.r[,4]), 1)
vennDiagram(vc.r)


df.c <- dcast(df.s, trvcdr3 ~ population, fun.aggregate = sum, value.var = "numReads")
cdr3 <- df.c$trvcdr3
df.c <- decostand(df.c[-1], "total", 2)
df.c$trvcdr3 <- cdr3

##bootstrap function

bootstrap.venn <- function(data, n=500, b=100, replace = F, weighted = T, as.freq = F){
  se <- function(x) sd(x)/sqrt(length(x))
  tally <- list()
  for (i in 1:b){
    if (replace == F){
      if (weighted == T) {
        FR.rs <- sample(data$trvcdr3, size = n, replace = F, prob = data$`Foxp3+RORgt+`)
        F.rs <- sample(data$trvcdr3, size = n, replace = F, prob = data$`Foxp3+`)
        R.rs <- sample(data$trvcdr3, size = n, replace = F, prob = data$`RORgt+`)
      } else {
        FR.rs <- sample(data$trvcdr3, size = n, replace = F)
        F.rs <- sample(data$trvcdr3, size = n, replace = F)
        R.rs <- sample(data$trvcdr3, size = n, replace = F)
      }
    } else {
      if (weighted == T) {
        FR.rs <- sample(data$trvcdr3, size = n, replace = T, prob = data$`Foxp3+RORgt+`)
        F.rs <- sample(data$trvcdr3, size = n, replace = T, prob = data$`Foxp3+`)
        R.rs <- sample(data$trvcdr3, size = n, replace = T, prob = data$`RORgt+`)
      } else {
        FR.rs <- sample(data$trvcdr3, size = n, replace = T)
        F.rs <- sample(data$trvcdr3, size = n, replace = T)
        R.rs <- sample(data$trvcdr3, size = n, replace = T)
      }
    }
    cdr3.FR <- levels(factor(FR.rs))
    cdr3.F <- levels(factor(F.rs))
    cdr3.R <- levels(factor(R.rs))
    combined <- levels(factor(c(cdr3.FR, cdr3.F, cdr3.R)))
    log.FR <- combined %in% cdr3.FR
    log.F <- combined %in% cdr3.F
    log.R <- combined %in% cdr3.R
    tally <- c(tally, list(vennCounts(cbind(log.FR, log.F, log.R))))
  }
  venn.values <- lapply(tally, function(x) x[,4])
  venn.matrix <- do.call("rbind", venn.values)
  if (as.freq == T){
    venn.matrix <- decostand(venn.matrix, "total", 1) * 100
  }  
  venn.mean <- tally[[1]]
  venn.se <- tally[[1]]
  venn.mean[,4] <- round(apply(venn.matrix, 2, mean), digits = 1)

  # return(venn.matrix)
  venn.se[,4] <- round(apply(venn.matrix, 2, se), digits = 2)
  par(mfrow=c(1,2), oma = c(0,0,5,0))
  vennDiagram(venn.mean, mar = c(6,0,5,0), cex = 0.9, main = "MEAN")
  vennDiagram(venn.se, mar = c(6,0,5,0), cex = 0.9, main = "SE")
  title <- paste0("Replace = ", replace, ", Weighted = ", weighted, ", \nPercentage = ", as.freq)
  mtext(title, outer = T, cex = 1.5)
# }

bootstrap.venn(df.c, n=500, b=100, replace = F, as.freq = T)





#Randomly sample TCRs according to frequency, without replacement
#Every population will get 700 unique TCRs, with those of higher frequencies more likely to appear on the list
#Similar to just taking top 700 TCRs for each sample?
FR.rs <- sample(df.c$trvcdr3, size = 700, replace = F, prob = df.c$`Foxp3+RORgt+`)
F.rs <- sample(df.c$trvcdr3, size = 700, replace = F, prob = df.c$`Foxp3+`)
R.rs <- sample(df.c$trvcdr3, size = 700, replace = F, prob = df.c$`RORgt+`)

cdr3.FR <- levels(factor(FR.rs))
cdr3.F <- levels(factor(F.rs))
cdr3.R <- levels(factor(R.rs))
combined <- levels(factor(c(cdr3.FR, cdr3.F, cdr3.R)))

log.FR <- combined %in% cdr3.FR
log.F <- combined %in% cdr3.F
log.R <- combined %in% cdr3.R
vc.rs <- vennCounts(cbind(log.FR, log.F, log.R))
# vennDiagram(vc.rs)
vc.rs.r <- vc.rs
vc.rs.r[,4] <- round(100*vc.rs[,4]/sum(vc.rs[,4]), 1)
vennDiagram(vc.rs.r)



#Randomly sample TCRs, regardless of frequency, without replacement
#Every population will get 700 unique TCRs, regardless of TCR frequency
FR.rs <- sample(df.c$trvcdr3, size = 700, replace = F)
F.rs <- sample(df.c$trvcdr3, size = 700, replace = F)
R.rs <- sample(df.c$trvcdr3, size = 700, replace = F)

cdr3.FR <- levels(factor(FR.rs))
cdr3.F <- levels(factor(F.rs))
cdr3.R <- levels(factor(R.rs))
combined <- levels(factor(c(cdr3.FR, cdr3.F, cdr3.R)))

log.FR <- combined %in% cdr3.FR
log.F <- combined %in% cdr3.F
log.R <- combined %in% cdr3.R
vc.rs <- vennCounts(cbind(log.FR, log.F, log.R))
vennDiagram(vc.rs)
vc.rs.r <- vc.rs
vc.rs.r[,4] <- round(100*vc.rs[,4]/sum(vc.rs[,4]), 1)
# vennDiagram(vc.rs.r)



#Randomly sample TCRs according to frequency, with replacement
#Every population will get 700 reads, with variable numbers of unique TCRs relative to the frequencies of each TCR
FR.rs <- sample(df.c$trvcdr3, size = 700, replace = T, prob = df.c$`Foxp3+RORgt+`)
F.rs <- sample(df.c$trvcdr3, size = 700, replace = T, prob = df.c$`Foxp3+`)
R.rs <- sample(df.c$trvcdr3, size = 700, replace = T, prob = df.c$`RORgt+`)

cdr3.FR <- levels(factor(FR.rs))
cdr3.F <- levels(factor(F.rs))
cdr3.R <- levels(factor(R.rs))
combined <- levels(factor(c(cdr3.FR, cdr3.F, cdr3.R)))

log.FR <- combined %in% cdr3.FR
log.F <- combined %in% cdr3.F
log.R <- combined %in% cdr3.R
vc.rs <- vennCounts(cbind(log.FR, log.F, log.R))
vennDiagram(vc.rs)
vc.rs.r <- vc.rs
vc.rs.r[,4] <- round(100*vc.rs[,4]/sum(vc.rs[,4]), 1)
vennDiagram(vc.rs.r)



#Randomly sample TCRs, regardless of frequency, with replacement
#Every population will get 700 reads, but variable numbers of unique TCRs, regardless of TCR frequency
FR.rs <- sample(df.c$trvcdr3, size = 700, replace = T)
F.rs <- sample(df.c$trvcdr3, size = 700, replace = T)
R.rs <- sample(df.c$trvcdr3, size = 700, replace = T)

cdr3.FR <- levels(factor(FR.rs))
cdr3.F <- levels(factor(F.rs))
cdr3.R <- levels(factor(R.rs))
combined <- levels(factor(c(cdr3.FR, cdr3.F, cdr3.R)))

log.FR <- combined %in% cdr3.FR
log.F <- combined %in% cdr3.F
log.R <- combined %in% cdr3.R
vc.rs <- vennCounts(cbind(log.FR, log.F, log.R))
vennDiagram(vc.rs)
vc.rs.r <- vc.rs
vc.rs.r[,4] <- round(100*vc.rs[,4]/sum(vc.rs[,4]), 1)
# vennDiagram(vc.rs.r)
