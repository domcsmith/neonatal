##########################################################################
# Preprocessing script for NeonatalData.xlsx
# Converts data from NeonatalData files into a form that can be used
# in classification analysis in R
##########################################################################

library("gdata")
xls <- read.xls(xls = "../data2014/source/NeonatalData.xlsx",
          sheet = "No TAM", method = "tab", header = TRUE, dec = ".")

# Drop columns at the end
xls <- xls[,1:36]
# Drop unwanted columns
unwanted.cols <- c("no", "Gender","New.Number", "Gata1", "Preg.med", "Med.class")
xls <- xls[,!(names(xls) %in% unwanted.cols)]

data.nt <- na.omit(apply(xls, 2, as.numeric))

print(paste("Imported ", dim(data.nt)[1], " from ", dim(xls)[1], " records"))

############
# TAM data
############

xls <- read.xls(xls = "../data2014/source/NeonatalData.xlsx",
                sheet = "TAM", method = "tab", header = TRUE, dec = ".")

# Drop columns at the end
xls <- xls[,1:36]
# Drop unwanted columns
unwanted.cols <- c("no", "Gender","New.Number", "Gata1", "Preg.med", "Med.class")

xls <- xls[,!(names(xls) %in% unwanted.cols)]

data.t <- na.omit(apply(xls, 2, as.numeric))

print(paste("Imported ", dim(data.t)[1], " from ", dim(xls)[1], " records"))

data.nt <- cbind(rep(0, dim(data.nt)[1]), data.nt)
colnames(data.nt)[1] <- "Class"
data.t <- cbind(rep(1, dim(data.t)[1]), data.t)
colnames(data.t)[1] <- "Class"
