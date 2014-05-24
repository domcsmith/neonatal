###################################################################
# Classification of neonatal data via machine learning techniques
# Dominic Smith 2014/05/24
###################################################################

setwd("~/academic/projects/infant_blood_project/")
#source("./code/libraryInstall.R")
source("./code/neonatalPreprocessing.R")

library("Matrix")
library("caret")

data <- rbind(data.t, data.nt)
data.df <- data.frame(data)

# take data from neo.df and create a response variable of classes "neo.y"
# and a data matrix neo.x of explanatory variables.
neo.y <- as.factor(data.df[,1])

neo.x <- data.df[,2:32]
neo.x.matrix <- sparse.model.matrix(object=Class~.,data=data.df)
neo.y <- as.numeric(data.df[,1])
neo.y.factor <- as.factor(neo.y)

# for K-fold cross validation, K = 10% * n
cv.number <- as.integer(round(0.1*length(neo.y),0))

fitControl <- trainControl(
  method = "repeatedcv",
  number = cv.number,
  ## number of repeats set to 10
  repeats = 10
  )

source("./code/ibp_svm.R")

svmFit
