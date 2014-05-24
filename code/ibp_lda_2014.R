# Classification analysis of Irene Roberts' 2014/05/07 neonatal DS data set
# Linear discriminant analysis
# MPV was ignored in this analysis owing to too many missing values
# Data were cleaned on 12/05/13

setwd("~/academic/projects/infant_blood_project/")
#source("./code/libraryInstall.R")
source("./code/neonatalPreprocessing.R")

library("MASS")
library("klaR")
library("glmnet")
library("denstrip")
library("Matrix")
library("fields")
library("e1071")

# bind norm and ds data frames on shared columns

data <- rbind(data.t, data.nt)
data.df <- data.frame(data)

#####################################################################
# Perform LDA on DS-DST
neo.ldafit <- lda(formula=Class ~ ., data=data.df, na.action=na.omit)
a <- cbind(neo.ldafit$scaling, colnames(data)[-1])
a[order(abs(as.numeric(a[,1]))),]

par(mar = rep(2, 4))
partimat(formula=Class ~ ., data=ds_dst.df, na.action=na.omit)

# Perform LDA on DS-DST; include jackknifing (cross validation)

# plotting
plot(neo.ldafit, panel=panel.lda, cex = 0.7, dimen=2,
     abbrev = FALSE, xlab = "LD1", ylab = "LD2", type="eqscplot")
#####################################################################

#####################################################################
#Next steps
# Elastic net regularized logistic regression on TAM and no-TAM
# uses glmnet and caret for training elastic net models
#####################################################################

library("glmnet")
library("caret")

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
                                    family="binomial",
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


######################################################################
#OLD
######################################################################

neo.y <- neo.y[as.numeric(rownames(neo.x))] # remove NA rows from Y


for(i in 0:10) {

alpha <- i/10

# run logistic regression with lasso regularisation on the full ds vs dst data set
neo.lognet_lasso_fit <- glmnet(neo.x,neo.y,
family="binomial", alpha=i)
# run cross-validated logistic regression with lasso regularisation for the ds vs dst data
neo.lognet_lasso_fit.cv <- cv.glmnet(neo.x,neo.y,
family="binomial", alpha=i)




}


png(paste("../output/figures/neo_ds-dst_cv_fit_lasso.png"), width=4, height=4, units="in", pointsize=10, res=200)
plot(neo.lognet_lasso_fit.cv, xlab=expression(log(lambda)), main=paste("Selecting ", expression(lambda), " through cross validation")) # plot log(lambda) against deviance for 10-fold cross validation
dev.off()


print(paste("Lambda min = ", toString(neo.lognet_lasso_fit.cv$lambda.min)))

# predict DS or DST for an example taken from dataset
s <- sample(x=nrow(neo.x),size=nrow(neo.x))
neo.x.pred_sample <- as.matrix(neo.x[s,])

# return coefficients for optimal lambda
y.response <- predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.matrix, s=neo.lognet_lasso_fit.cv$lambda.min, type="response")
y.class <- predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.matrix, s=neo.lognet_lasso_fit.cv$lambda.min, type="class")
y.actual <- as.matrix(neo.y)

# show misclassifications and response probabilities
misclassification.comparison <- cbind(y.samp.class, y.samp.actual, y.samp.response)
print(misclassification.comparison)
misclassification.comparison.sorted <- misclassification.comparison[sort.list(as.numeric(misclassification.comparison[,3])),]

dummy.x <- seq(1,length(misclassification.comparison[,3]))

misclassified.y <- misclassification.comparison.sorted[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2],]

misclassified.x <- as.numeric(dummy.x[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2]])

png("../output/figures/neo_ds_dst_sigmoid.png", width=4, height=4, units="in", res=200, pointsize=10)
plot(dummy.x, misclassification.comparison.sorted[,3], pch=1, xlab="Rank of Sorted Response", ylab="Response")
points(misclassified.x, as.numeric(misclassified.y[,3]), pch=16, col="red")
lines(x=c(-100,200), y=c(0.5,0.5), lty=2, col="gray")
text(misclassified.x, as.numeric(misclassified.y[,3]), misclassified.y[,2], cex=0.6, pos=4)
text(x=c(20),y=c(0.2), "DSnT", pos="4", col="blue")
text(x=c(100),y=c(0.8), "DST", pos="4", col="blue")
dev.off()
rm("dummy.x", "dummy.y", "misclassified.x", "misclassified.y")

ds_dst.coeffs <- predict(neo.lognet_lasso_fit,
                         newx=neo.x, s=neo.lognet_lasso_fit.cv$lambda.min, type="coefficients")


ncoeffs<-row.names(ds_dst.coeffs)[3:28]

ds_dst.coeffs <- ds_dst.coeffs[3:28]

lda_coeffs<-neo.ldafit$scaling

############ PLOT COEFFICIENTS OF DDST ANALYSIS #############
png("../output/figures/ds_dst_coeffs.png", width=4, height=4, units="in", res=400, pointsize=10)
plot(ds_dst.coeffs, ylim=c(-4, 3), ylab="Coefficient Value", xlab="", axes=FALSE)
points(lda_coeffs, col="blue")
axis(2)
text(ncoeffs, x=1:26, y=rep(-2,26), pos=1, cex=0.6, srt=90)
lines(x=c(-1,28), y=c(0,0), col="black")
dev.off()

# store columns with non-zero coefficients from lasso logistic regression for each variable
subset.vars <- rownames(predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.pred_sample, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef"))[as.vector(predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.pred_sample, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef")!=0)]

# Run LDA on subset variables taken from logistic regression
ds_dst_subset.df <- ds_dst.df[,c("Class", "Preg_med", "MCV", "Blast", "AbnRBC", "AbnNeut")]
neo_subset.ldafit <- lda(formula=Class ~ ., data=ds_dst_subset.df, na.action=na.omit)
plot(neo_subset.ldafit)
plot(neo_subset.ldafit, type="density")

subset_vars.lognet_data <- as.matrix(predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.pred_sample, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef"))

names <- as.vector(rownames(neo_subset.ldafit$scaling))

print("The coefficients of the lasso logistic regression and the LDA function for DS vs DST")
print(rbind(subset_vars.lognet_data[names,], neo_subset.ldafit$scaling[names,]))
rm("names")

png(paste("../output/figures/partition_plot_subset_ds_dst.png"), width=4, height=4, units="in", res=1600, pointsize=4)
par(cex=1)
partimat(formula=Class ~ ., data=ds_dst_subset.df, na.action=na.omit, nplots.vert=5, nplots.hor=5, mar=rep(2,4), name=colnames(ds_dst_subset.df))
dev.off()

################### SUBSAMPLING CODE ####################
ds_dst.df <- read.csv("neo_data_cleaned20130514_ds-dst.csv", header=TRUE, sep=",", dec=".")

count <- c() # initialize cumulative counter
total <- c() # initialize counter denominator
count.false.ds <- c()
count.false.dst <- c()
total.resp.ds <- c()
total.resp.dst <- c()
min.sigmoid <- c()
max.sigmoid <- c()
sigmoid.false.ds <- c()
sigmoid.false.dst <- c()
count.min <- 0
count.max <- 0
vars <- list()
count.all.dst <- rep(0,50)
count.all.ds <- rep(0,50)
count.mis.dst <- rep(0,50)
count.mis.ds <- rep(0,50)


for(i in 1:2000) {
  s <- sample(1:nrow(ds_dst.df), size=nrow(ds_dst.df), replace=FALSE)
  n <- round(3*length(s)/4,0)
  ds_dst.df.train <- ds_dst.df[s[1:n],]
  ds_dst.df.pred <- ds_dst.df[s[n+1:length(s)],]
  
  neo.x.train <- sparse.model.matrix(object=Class~.,data=ds_dst.df.train, na.action=na.omit)
  neo.x.pred <- sparse.model.matrix(object=Class~.,data=ds_dst.df.pred, na.action=na.omit)
  neo.y <- as.factor(ds_dst.df[,1])
  
  neo.y.train <- neo.y[as.numeric(rownames(neo.x.train))] # remove NA rows from Y
  neo.y.pred <- neo.y[as.numeric(rownames(neo.x.pred))]
  
  # run logistic regression with lasso regularisation on the full ds vs dst data set
  neo.lognet_lasso_fit <- glmnet(neo.x.train,neo.y.train,
                                 family="binomial", alpha=alpha)
  # run cross-validated logistic regression with lasso regularisation for the ds vs dst data
  neo.lognet_lasso_fit.cv <- cv.glmnet(neo.x.train,neo.y.train,
                                       family="binomial", alpha=alpha)
  
  # print(paste("Lambda min = ", toString(neo.lognet_lasso_fit.cv$lambda.min)))
  
  # return coefficients for optimal lambda
  y.samp.response <- predict(neo.lognet_lasso_fit,
                             newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="response")
  y.samp.class <- predict(neo.lognet_lasso_fit,
                          newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="class")
  y.samp.actual <- as.matrix(neo.y.pred)
  
  # show misclassifications and response probabilities
  misclassification.comparison <- cbind(y.samp.class, y.samp.actual, y.samp.response)
  # print(misclassification.comparison)
  misclassification.comparison.sorted <- misclassification.comparison[sort.list(as.numeric(misclassification.comparison[,3])),]
  
  dummy.x <- seq(1,length(misclassification.comparison[,3]))
  
  misclassified.y <- misclassification.comparison.sorted[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2],]
  misclassified.x <- as.numeric(dummy.x[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2]])
  
  all.ds <- as.numeric(misclassification.comparison.sorted[misclassification.comparison.sorted[,1]=="DSnT",3])
  all.dst <- as.numeric(misclassification.comparison.sorted[misclassification.comparison.sorted[,1]=="DST",3])
  
  count.all.dst<-count.all.dst + stats.bin(all.dst, all.dst, breaks=seq(0.5,1,.01))$stats["N",]
  count.all.ds<-count.all.ds + stats.bin(all.ds, all.ds, breaks=seq(0,0.5,.01))$stats["N",]
  
  count <- c(count,length(misclassified.x))
  if(length(misclassified.x)==1)
  {
    count.false.ds <- c(count.false.ds,sum(as.numeric(misclassified.y[1]=="DSnT")))
    count.false.dst <- c(count.false.dst,sum(as.numeric(misclassified.y[1]=="DST")))
    count.mis.ds<-count.mis.ds + stats.bin(misclassified.y[3], misclassified.y[3], breaks=seq(0,0.5,.01))$stats["N",]
    count.mis.dst<-count.mis.dst + stats.bin(misclassified.y[3], misclassified.y[3], breaks=seq(0.5,1,.01))$stats["N",]
    if(misclassified.y[1]=="DSnT")
    {
      min.sigmoid <- c(min.sigmoid, as.numeric(misclassified.y[3]))
      count.min <- count.min+1
      
    }
    
    if(misclassified.y[1]=="DST")
    {
      max.sigmoid <- c(max.sigmoid, as.numeric(misclassified.y[3]))
      count.max <- count.max+1
    }
  }
  else if(length(misclassified.x)>1)
  {
    count.false.ds <- c(count.false.ds,sum(as.numeric(misclassified.y[,1]=="DSnT")))
    count.false.dst <- c(count.false.dst,sum(as.numeric(misclassified.y[,1]=="DST")))
    count.mis.ds<-count.mis.ds + stats.bin(misclassified.y[,3], misclassified.y[,3], breaks=seq(0,0.5,.01))$stats["N",]
    count.mis.dst<-count.mis.dst + stats.bin(misclassified.y[,3], misclassified.y[,3], breaks=seq(0.5,1,.01))$stats["N",]
    
    if(misclassified.y[1,1]=="DSnT")
    {
      min.sigmoid <- c(min.sigmoid, as.numeric(misclassified.y[1,3]))
      sigmoid.false.ds <- c(sigmoid.false.ds, as.numeric(misclassified.y[misclassified.y[,1]=="DSnT",3]))
      count.min <- count.min+1
    }
    
    if(misclassified.y[length(misclassified.x),1]=="DST")
    {
      max.sigmoid <- c(max.sigmoid, as.numeric(misclassified.y[length(misclassified.x),3]))
      sigmoid.false.norm <- c(sigmoid.false.dst, as.numeric(misclassified.y[misclassified.y[,1]=="DST",3]))
      count.max <- count.max+1
    }
  }
  else
  {
    count.false.ds <- c(count.false.ds,0)
    count.false.dst <- c(count.false.dst,0)
  }
  
  total <- c(total,nrow(neo.x.pred))
  total.resp.ds <- c(total.resp.ds,sum(as.numeric(y.samp.class=="DSnT")))
  total.resp.dst <- c(total.resp.dst,sum(as.numeric(y.samp.class=="DST")))
  
  # store columns with non-zero coefficients from lasso logistic regression for each variable
  subset.vars <- rownames(predict(neo.lognet_lasso_fit.cv,
                                  newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef"))[as.vector(predict(neo.lognet_lasso_fit.cv,
                                                                                                                         newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef")!=0)]
  subset.vars <- subset.vars[2:length(subset.vars)]
  vars[[i]] <- subset.vars # gather a list of selected variables
  # to be used in model stability study
  
  print(i)
}

error.pct <- count/total
error.pct.resp.ds <- count.false.ds/total.resp.ds
error.pct.resp.dst <- count.false.dst/total.resp.dst
max.s <- max.sigmoid
min.s <- min.sigmoid
c.max <- count.max
c.min <- count.min

tabDSnT <- data.frame(R=seq(0.005,0.495,0.01), freq=count.mis.ds/count.all.ds)
tabDST <- data.frame(R=seq(0.505,0.995,.01), freq=count.mis.dst/count.all.dst)

tabDSnT <- transform(tabDSnT, freq=freq/sum(freq))
tabDST <- transform(tabDST, freq=freq/sum(freq))


##############################################################################################
# norm-DS analysis
norm_ds.df <- read.csv("neo_data_cleaned20130516_norm-ds.csv", header=TRUE, sep=",", dec=".")
#norm_ds.df <- norm_ds.df[,index]
#rm("index")

computeC <- function(vector){
    r=sum(scale*vector)
    return(r)
}

# Perform LDA on Norm-DS
neo.norm_ds.ldafit <- lda(formula=Class ~ ., data=norm_ds.df, na.action=na.omit)
neo.x <- model.matrix(object=Class~.,data=norm_ds.df, na.action=na.omit)
neo.x <- neo.x[,2:ncol(neo.x)]
lda_coeffs<-neo.norm_ds.ldafit$scaling
#lda_out <- neo.x %*% coeff

dall1<-predict(neo.norm_ds.ldafit, norm_ds.df)$x

dnorm <- density(dall1[1:123], na.rm=T)
dds <- density(dall1[124:323], na.rm=T)

# Perform LDA on DS-DST
ds_dst.df <- read.csv("neo_data_cleaned20130514_ds-dst.csv", header=TRUE, sep=",", dec=".")
neo.ds_dst.ldafit <- lda(formula=Class ~ ., data=ds_dst.df, na.action=na.omit)

dall2<-predict(neo.ds_dst.ldafit, ds_dst.df)$x
ddsnt<-density(dall2[1:183],na.rm=T)
ddst<-density(dall2[184:200],na.rm=T)

png("../output/figures/lda_all_dens.png", width=4, height=2, units="in", pointsize=12)
par(mfrow=c(1,2))
plot(dds, col="blue", xlim=c(-4, 6), main="Normal and DS groups", xlab="Linear discriminant")
lines(dnorm, col="red")
text(x=c(-1.5),y=c(0.55), "Normal", pos=4)
text(x=c(5),y=c(0.25), "DS", pos=4)
plot(ddsnt, col="blue", xlim=c(-4, 15), main="DST and DSnT groups", xlab="Linear discriminant")
lines(ddst, col="red")
text(x=c(0),y=c(0.55), "DSnT", pos=4)
text(x=c(5.5),y=c(0.12), "DST", pos=4)
dev.off()

misclass.n_ds <- cbind(predict(neo.norm_ds.ldafit, norm_ds.df)$class, norm_ds.df[,1])
misclass.ds_dst <- cbind(predict(neo.ds_dst.ldafit, ds_dst.df)$class, ds_dst.df[,1])


############
# repeat analysis training on half the data set and predicting on the other half
s <- sample(1:nrow(norm_ds.df), size=nrow(norm_ds.df), replace=FALSE)
n <- round(length(s)/2,0)
norm_ds.df.train <- norm_ds.df[s[1:n],]
norm_ds.df.pred <- norm_ds.df[s[n+1:length(s)],]
############

#Next steps
# REGULARIZED LOGISTIC REGRESSION ON NORM-DS DATA CLASS
# uses glmnet
########################
# SET REGRESSION PARAMETERS
norm_ds.df <- read.csv("neo_data_cleaned20130516_norm-ds.csv", header=TRUE, sep=",", dec=".")
alpha <- 1 # Elastic net mixing parameter range [0,1): 1 = lasso, 0 = ridge
########################

neo.x <- sparse.model.matrix(object=Class~.,data=norm_ds.df, na.action=na.omit)
neo.x <- neo.x[,2:ncol(neo.x)]
neo.y <- as.factor(norm_ds.df[,1])

neo.y <- neo.y[as.numeric(rownames(neo.x))] # remove NA rows from Y

# run logistic regression with lasso regularisation on the full ds vs dst data set
neo.lognet_lasso_fit <- glmnet(neo.x,neo.y,
                               family="binomial", alpha=alpha)
# run cross-validated logistic regression with lasso regularisation for the ds vs dst data
neo.lognet_lasso_fit.cv <- cv.glmnet(neo.x,neo.y,
                                     family="binomial", alpha=alpha)

norm_ds.coeffs <- predict(neo.lognet_lasso_fit,
                           newx=neo.x, s=neo.lognet_lasso_fit.cv$lambda.min, type="coefficients")

ncoeffs<-row.names(norm_ds.coeffs)[2:15]

norm_ds.coeffs <- norm_ds.coeffs[2:15]

############ PLOT COEFFICIENTS OF NORM-DS ANALYSIS #############
png("../output/figures/norm_ds_coeffs.png", width=4, height=4, units="in", res=400, pointsize=10)
plot(norm_ds.coeffs, ylim=c(-8, 5), ylab="Coefficient Value", xlab="", axes=FALSE)
points(lda_coeffs, col="blue")
axis(2)
text(ncoeffs, x=1:14, y=norm_ds.coeffs, pos=1, cex=0.6)
lines(x=c(-1,15), y=c(0,0), col="black")
dev.off()

library("glmnet")

# take data from neo.df and create a response variable of classes "neo.y"
# and a data matrix neo.x of explanatory variables.

count <- c() # initialize cumulative counter
total <- c() # initialize counter denominator
count.false.ds <- c()
count.false.norm <- c()
total.resp.ds <- c()
total.resp.norm <- c()
min.sigmoid <- c()
max.sigmoid <- c()
sigmoid.false.ds <- c()
sigmoid.false.norm <- c()
count.min <- 0
count.max <- 0
vars <- list()
count.all.norm.cum <- rep(0,50)
count.all.ds.cum <- rep(0,50)
count.mis.norm.cum <- rep(0,50)
count.mis.ds.cum <- rep(0,50)

for(i in 1:500) {
s <- sample(1:nrow(norm_ds.df), size=nrow(norm_ds.df), replace=FALSE)
n <- round(length(s)/2,0)
norm_ds.df.train <- norm_ds.df[s[1:n],]
norm_ds.df.pred <- norm_ds.df[s[n+1:length(s)],]

neo.x.train <- sparse.model.matrix(object=Class~.,data=norm_ds.df.train, na.action=na.omit)
neo.x.pred <- sparse.model.matrix(object=Class~.,data=norm_ds.df.pred, na.action=na.omit)
neo.y <- as.factor(norm_ds.df[,1])

neo.y.train <- neo.y[as.numeric(rownames(neo.x.train))] # remove NA rows from Y
neo.y.pred <- neo.y[as.numeric(rownames(neo.x.pred))]

# run logistic regression with lasso regularisation on the full ds vs dst data set
neo.lognet_lasso_fit <- glmnet(neo.x.train,neo.y.train,
family="binomial", alpha=alpha)
# run cross-validated logistic regression with lasso regularisation for the ds vs dst data
neo.lognet_lasso_fit.cv <- cv.glmnet(neo.x.train,neo.y.train,
family="binomial", alpha=alpha)

# print(paste("Lambda min = ", toString(neo.lognet_lasso_fit.cv$lambda.min)))

# return coefficients for optimal lambda
y.samp.response <- predict(neo.lognet_lasso_fit,
               newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="response")
y.samp.class <- predict(neo.lognet_lasso_fit,
               newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="class")
y.samp.actual <- as.matrix(neo.y.pred)

# show misclassifications and response probabilities
misclassification.comparison <- cbind(y.samp.class, y.samp.actual, y.samp.response)
# print(misclassification.comparison)
misclassification.comparison.sorted <- misclassification.comparison[sort.list(as.numeric(misclassification.comparison[,3])),]

dummy.x <- seq(1,length(misclassification.comparison[,3]))

misclassified.y <- misclassification.comparison.sorted[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2],]
misclassified.x <- as.numeric(dummy.x[misclassification.comparison.sorted[,1]!=misclassification.comparison.sorted[,2]])

all.norm <- as.numeric(misclassification.comparison.sorted[misclassification.comparison.sorted[,1]=="Norm",3])
all.ds <- as.numeric(misclassification.comparison.sorted[misclassification.comparison.sorted[,1]=="DS",3])

count.all.norm<-count.all.norm + stats.bin(all.norm, all.norm, breaks=seq(0.5,1,.01))$stats["N",]
count.all.ds<-count.all.ds + stats.bin(all.ds, all.ds, breaks=seq(0,0.5,.01))$stats["N",]

count <- c(count,length(misclassified.x))
if(length(misclassified.x)==1)
{
  count.false.ds <- c(count.false.ds,sum(as.numeric(misclassified.y[1]=="DS")))
  count.false.norm <- c(count.false.norm,sum(as.numeric(misclassified.y[1]=="Norm")))
  count.mis.ds<-count.mis.ds + stats.bin(misclassified.y[3], misclassified.y[3], breaks=seq(0,0.5,.01))$stats["N",]
  count.mis.norm<-count.mis.norm + stats.bin(misclassified.y[3], misclassified.y[3], breaks=seq(0.5,1,.01))$stats["N",]
  if(misclassified.y[1]=="DS")
  {
    min.sigmoid <- c(min.sigmoid, as.numeric(misclassified.y[3]))
    count.min <- count.min+1
    
  }

  if(misclassified.y[1]=="Norm")
  {
    max.sigmoid <- c(max.sigmoid, as.numeric(misclassified.y[3]))
    count.max <- count.max+1
  }
}
else if(length(misclassified.x)>1)
{
  count.false.ds <- c(count.false.ds,sum(as.numeric(misclassified.y[,1]=="DS")))
  count.false.norm <- c(count.false.norm,sum(as.numeric(misclassified.y[,1]=="Norm")))
  count.mis.ds<-count.mis.ds + stats.bin(misclassified.y[,3], misclassified.y[,3], breaks=seq(0,0.5,.01))$stats["N",]
  count.mis.norm<-count.mis.norm + stats.bin(misclassified.y[,3], misclassified.y[,3], breaks=seq(0.5,1,.01))$stats["N",]
  
  if(misclassified.y[1,1]=="DS")
  {
    min.sigmoid <- c(min.sigmoid, as.numeric(misclassified.y[1,3]))
    sigmoid.false.ds <- c(sigmoid.false.ds, as.numeric(misclassified.y[misclassified.y[,1]=="DS",3]))
    count.min <- count.min+1
  }
  
  if(misclassified.y[length(misclassified.x),1]=="Norm")
  {
    max.sigmoid <- c(max.sigmoid, as.numeric(misclassified.y[length(misclassified.x),3]))
    sigmoid.false.norm <- c(sigmoid.false.norm, as.numeric(misclassified.y[misclassified.y[,1]=="Norm",3]))
    count.max <- count.max+1
  }
}
else
{
  count.false.ds <- c(count.false.ds,0)
  count.false.norm <- c(count.false.norm,0)
}

total <- c(total,nrow(neo.x.pred))
total.resp.ds <- c(total.resp.ds,sum(as.numeric(y.samp.class=="DS")))
total.resp.norm <- c(total.resp.norm,sum(as.numeric(y.samp.class=="Norm")))

# store columns with non-zero coefficients from lasso logistic regression for each variable
subset.vars <- rownames(predict(neo.lognet_lasso_fit.cv,
                                newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef"))[as.vector(predict(neo.lognet_lasso_fit.cv,
                                                                                                                       newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef")!=0)]
subset.vars <- subset.vars[2:length(subset.vars)]
vars[[i]] <- subset.vars # gather a list of selected variables
# to be used in model stability study

print(i)
}

tabDS <- data.frame(R=seq(0.005,0.495,0.01), freq=count.mis.ds/count.all.ds)
tabNorm <- data.frame(R=seq(0.505,0.995,.01), freq=count.mis.norm/count.all.norm)

tabDS <- transform(tabDS, freq=freq/sum(freq))
tabNorm <- transform(tabNorm, freq=freq/sum(freq))

# PLOT sigmoid with misclassification densities
# change DSnT to DS and DST to Norm to plot for the other comparison
# you will need to recalculate the densities at which the misclassification
# probability falls below 0.5% and enter those into the "lines" code
library("denstrip")
png("../output/figures/neo_ds_dst_sigmoid.png", width=4, height=4, units="in", res=400, pointsize=10)
plot(dummy.x, misclassification.comparison.sorted[,3], pch=1, xlab="Response Rank", ylab="Response")
denstrip(x=tabDSnT$R, dens=tabDSnT$freq, horiz=FALSE, at=length(dummy.x)/2, width=length(dummy.x), gamma=1, colmax="orange", kernel="gaussian", cut=0)
denstrip(x=tabDST$R, dens=tabDST$freq, horiz=FALSE, at=length(dummy.x)/2, width=length(dummy.x), gamma=1, colmax="lightblue", kernel="gaussian", cut=0)
points(dummy.x, misclassification.comparison.sorted[,3], pch=1)
points(misclassified.x, as.numeric(misclassified.y[,3]), pch=16, col="red")
lines(c(0,length(dummy.x)), c(0.995,0.995), lty=2, col="black")
lines(c(0,length(dummy.x)), c(0.075,0.075), lty=2, col="black")
text(misclassified.x, as.numeric(misclassified.y[,3]), misclassified.y[,2], cex=0.6, pos=4)
text(25,0.8,"DST", cex=1.5, pos=4)
text(10,0.3,"DSnT", cex=1.5, pos=4)
dev.off()


error.pct <- count/total
error.pct.resp.ds <- count.false.ds/total.resp.ds
error.pct.resp.norm <- count.false.norm/total.resp.norm
max.s <- max.sigmoid
min.s <- min.sigmoid
c.max <- count.max
c.min <- count.min

# response level for misclassification prob of 0.5% for DS
print(quantile(sigmoid.false.ds, probs=c(.005/mean(error.pct.resp.ds))))
print(quantile(sigmoid.false.norm, probs=c((1-(.005/mean(error.pct.resp.norm)))))

ecdf(sigmoid.false.norm)
quantile(sigmoid.false.norm, probs=(0.999))

png(paste("../output/figures/neo_norm-ds_cv_fit_alpha_lasso.png"), width=4, height=4, units="in", res=800, pointsize=10)
plot(neo.lognet_lasso_fit.cv, xlab=expression(log(lambda))) # plot log(lambda) against deviance for 10-fold cross validation
dev.off()

png(paste("../output/figures/neo_norm-ds_lambda_fit_alpha_", toString(i), ".png"), width=800, height=800, units="px")
plot(neo.lognet_lasso_fit, xvar="lambda") # plot lambda against coeffs for full data fit
dev.off()

library("denstrip")
png("../output/figures/neo_norm_ds_sigmoid.png", width=4, height=4, units="in", res=400, pointsize=10)
plot(dummy.x, misclassification.comparison.sorted[,3], pch=1, xlab="Response Rank", ylab="Response")
denstrip(x=tabDS$R, dens=tabDS$freq, horiz=FALSE, at=length(dummy.x)/2, width=length(dummy.x), gamma=1, colmax="orange", kernel="gaussian", cut=0)
denstrip(x=tabNorm$R, dens=tabNorm$freq, horiz=FALSE, at=length(dummy.x)/2, width=length(dummy.x), gamma=1, colmax="lightblue", kernel="gaussian", cut=0)
points(dummy.x, misclassification.comparison.sorted[,3], pch=1)
points(misclassified.x, as.numeric(misclassified.y[,3]), pch=16, col="red")
lines(c(0,length(dummy.x)), c(0.895,0.895), lty=2, col="black")
lines(c(0,length(dummy.x)), c(0.095,0.095), lty=2, col="black")
text(misclassified.x, as.numeric(misclassified.y[,3]), misclassified.y[,2], cex=0.6, pos=4)
dev.off()

p <- vector(mode="character", length=500)
for (i in 1:500)
{
  q <- sort(vars[[i]], decreasing=FALSE)
  p[i]<-paste(q, collapse="")
}

sort(table(p), decreasing=TRUE)

a <- as.vector(unlist(vars))
png("../output/figures/ds_dst_variable_selection_freq.png", width=4, height=4, units="in", res=400, pointsize=10)
plot(sort(table(a), decreasing=TRUE), ylim=c(0,2000), xlab="Selection Frequency Rank", ylab="Selection Count")
text(row.names(sort(table(a), decreasing=TRUE)), x=1:26, y=sort(table(a), decreasing=TRUE), pos=1, cex=0.3)
dev.off()

plot(ecdf(x=sigmoid.false.ds))

plot(density(sigmoid.false.norm, kernel="cosine", cut=0))

dev.off()
rm("dummy.x", "misclassified.x", "misclassified.y")

# store columns with non-zero coefficients from lasso logistic regression for each variable
subset.vars <- rownames(predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef"))[as.vector(predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.pred, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef")!=0)]
subset.vars <- subset.vars[2:length(subset.vars)]
vars <- c(vars, subset.vars)

subset.vars <- c("Class", "GA","Hb","MCV","Blast","Neut","Bas","Eos","GP")


# Run LDA on subset variables taken from logistic regression
norm_ds_subset.df <- norm_ds.df[,subset.vars[1:length(subset.vars)]]
neo_subset.ldafit <- lda(formula=Class ~ ., data=norm_ds_subset.df, na.action=na.omit)
plot(neo_subset.ldafit)
plot(neo_subset.ldafit, type="density")

subset_vars.lognet_data <- as.matrix(predict(neo.lognet_lasso_fit.cv,
               newx=neo.x.pred_sample, s=neo.lognet_lasso_fit.cv$lambda.min, type="coef"))


# plot a 1d graph demonstrating the relative size of response to each explanatory variable
names <- as.vector(rownames(neo_subset.ldafit$scaling))
y <- rep(0, length(names))
png("../output/figures/neo_norm_ds_1d_ldafit_subset.png", width=4, height=2, units="in", res=1600, pointsize=4)
plot(neo_subset.ldafit$scaling, y, pch=18, col="blue", ylim=c(-0.2,0.2), xlim=c(-1.5,1.5))
text(neo_subset.ldafit$scaling, y, names, cex=0.6, pos=3, col ="black") 
axis(side=2, pos=0, labels=FALSE)

print("The coefficients of the lasso logistic regression and the LDA function for DS vs DST")
print(rbind(subset_vars.lognet_data[names,], neo_subset.ldafit$scaling[names,]))
rm("names")



##############################################################################################

# Hct-only analysis
# NB: separated from other variables because records were from different patients
hct.df <- read.csv("neo_data_cleaned20130516_HCTonly.csv", header=TRUE, sep=",", dec=".")

hct.norm <- as.vector(hct.df[1:118,1])
hct.ds <- as.vector(hct.df[,2])
hct.dst <- as.vector(hct.df[1:17,3])


d.n <- density(hct.norm)
d.ds <- density(hct.ds)
d.dst <- density(hct.dst)
d.dsdst <- density(c(hct.ds,hct.dst))

png("../output/figures/neo_hct_histograms.png", width=4, height=4, units="in", res=400, pointsize=10)
plot(density(hct.norm), main="Hct", xlab="Hct", xlim=c(0,1), ylim=c(0, max(c(d.n$y, d.ds$y, d.dst$y, d.dsdst$y))), col="black")
lines(d.ds, col="blue")
lines(d.dst, col="red")
lines(d.dsdst, col="green")
dev.off()

# test the stochastic equivalence of Hct for each dataset using the Wilcoxon rank sum test
print(wilcox.test(hct.norm, hct.ds, alternative="two.sided"))
print(wilcox.test(hct.norm, hct.dst, alternative="two.sided"))
print(wilcox.test(hct.ds, hct.dst, alternative="two.sided"))
print(wilcox.test(hct.norm, rbind(hct.ds, hct.dst), alternative="two.sided"))
