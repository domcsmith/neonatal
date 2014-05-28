###################################################################
# Classification of neonatal data via logistic regression
# Dominic Smith 2014/05/28
###################################################################

library("MASS")
library("caret")
library("Matrix")
library("glmnet")
library("e1071")

fitControl <- trainControl(
  method = "repeatedcv",
  number = cv.number,
  ## number of repeats set to 10
  repeats = 10
  )

tuningGrid = expand.grid(alpha=(1:20)*0.05,
                         lambda=exp(seq(from=-10,to=1,length.out=10))
                         )

training.model <- train(neo.x, neo.y.factor,
                        method="glmnet",
                        trControl=fitControl,
                        metric="Kappa",
                        tuneGrid=tuningGrid,
                        maximize=TRUE)

neo.lognet_elasticnet_fit <- glmnet(neo.x.matrix,neo.y.factor,
                                    family="multinomial",
                                    alpha=as.data.frame(training.model[6])[1,1],
                                    lambda=as.data.frame(training.model[6])[1,2]
                                    )

# return coefficients for optimal lambda
y.response <- predict(neo.lognet_elasticnet_fit,
                      newx=neo.x.matrix,
                      type="response")

y.class <- predict(neo.lognet_elasticnet_fit,
                   newx=neo.x.matrix,
                   type="class")

y.actual <- as.matrix(neo.y)

# show misclassifications and response probabilities
misclassification.comparison <- cbind(y.class,
                                      y.actual,
                                      y.response)
print(misclassification.comparison)

misclassification.comparison.sorted <- misclassification.comparison[sort.list(as.numeric(misclassification.comparison[,3])),]

dummy.x <- seq(1,length(misclassification.comparison[,3]))

misclassified.y <- misclassification.comparison.sorted[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2],]

misclassified.x <- as.numeric(dummy.x[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2]])

png("neo_ds_dst_sigmoid.png",
    width=4,
    height=4,
    units="in",
    res=200,
    pointsize=10)

plot(dummy.x,
     misclassification.comparison.sorted[,3],
     pch=1,
     xlab="Rank of Sorted Response",
     ylab="Response")

points(misclassified.x,
       as.numeric(misclassified.y[,3]),
       pch=16,
       col="red")

lines(x=c(-100,length(neo.y)),
      y=c(0.5,0.5),
      lty=2,
      col="gray")

text(misclassified.x,
     as.numeric(misclassified.y[,3]),
     misclassified.y[,2],
     cex=0.6,
     pos=4)

text(x=c(20),y=c(0.2), "No TAM", pos="4", col="blue")
text(x=c(100),y=c(0.8), "TAM", pos="4", col="blue")

dev.off()
rm("dummy.x", "dummy.y", "misclassified.x", "misclassified.y")

