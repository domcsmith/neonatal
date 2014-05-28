###################################################################
# Classification of neonatal data via support vector machines
# Dominic Smith 2014/05/24
###################################################################

library("e1071")
library("caret")
library("kernlab")
library("pROC")

fitControl <- trainControl(
  method = "repeatedcv",
  number = cv.number,
  ## number of repeats set to 10
  repeats = 10,
  classProbs=TRUE
  )

svmFit <- train(neo.x, neo.y.factor,
                method="svmRadial",
                trControl=fitControl,
                preProc = c("center", "scale"),
                tuneLength = 8,
                metric = "ROC")

