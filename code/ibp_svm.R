###################################################################
# Classification of neonatal data via support vector machines
# Dominic Smith 2014/05/24
###################################################################

library("e1071")
library("caret")
library("kernlab")

svmFit <- train(neo.x, neo.y.factor,
                method="svmRadial",
                trControl=fitControl,
                preProc = c("center", "scale"),
                tuneLength = 8,
                metric = "ROC")

