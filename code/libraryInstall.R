##########################################################################
# Package installer script for neonatal data code
##########################################################################

r <- "http://cran.ma.imperial.ac.uk"
install.packages("gdata", repos=r)
install.packages("MASS", repos=r)
install.packages("klaR", repos=r)
install.packages("glmnet", repos=r)
install.packages("denstrip", repos=r)
install.packages("Matrix", repos=r)
install.packages("fields", repos=r)
install.packages("e1071", repos=r)
install.packages("caret", repos=r)
install.packages("kernlab", repos=r)

rm(r)
