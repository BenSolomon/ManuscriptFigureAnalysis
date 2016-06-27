library(caret); library(reshape2); library(vegan)

df <- read.csv("D:\\Dropbox\\Hsieh Lab\\Data Analysis\\Foxp3 RORgt Figures\\CleanedFoxp3RORgt.csv", header = T)


df.1 <- subset(df, tissue == "Colon")
df.c <- dcast(df.1, trvcdr3 ~ mouseNum + tissue + population, fun.aggregate = sum, value.var = "numReads")
row.names(df.c) <- df.c$trvcdr3
df.c <- as.data.frame(t(df.c[,-1]))
df.c <- decostand(df.c, "total", 1)
df.c$population <- factor(namesplit.f(rownames(df.c), "_")[,3])

# rf <- randomForest(population ~ ., data=df.c,
#                          importance=TRUE)
# varImpPlot(rf)
# 
# inTrain <- createDataPartition(y=df.c$population,
#                                p=0.75, list=FALSE)
# training <- spam[inTrain,]
# testing <- spam[-inTrain,]
# dim(training)
# set.seed(32343)
# 
modelFit <- train(df.c[,-8768], factor(df.c[,8768]), method="rf")
imp.df <- varImp(modelFit)$importance
imp.df$tcr <- rownames(imp.df)
write.csv(imp.df, "Colon TCR imp.csv", row.names = F)
# x <- predict(modelFit, df.c[,-8768])
# table(df.c$population, x)
# 
# rf <- randomForest(df.c[,-8769], factor(df.c[,8768]), importance = T)

# 
# modelFit
# modelFit$finalModel
# predictions <- predict(modelFit, newdata = testing)
# predictions
# confusionMatrix(predictions, testing$type)
# 
# 
# x <- data.frame("disease" = sample(c(TRUE, FALSE), 100, T), "test" = sample(c(TRUE, FALSE), 100, T))
# table(x)
# confusionMatrix(x$disease, x$test)
# confusionMatrix(x$test, x$disease, positive = "TRUE")


i <- 1
modelFit[i]
i <- i + 1