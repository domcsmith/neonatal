##########################################################################
# Preprocessing script for NeonatalData.xlsx
# Converts data from NeonatalData files into a form that can be used
# in classification analysis in R
##########################################################################

library("gdata")
xls <- read.xls(xls = "~/Documents/academic/infant_blood_project/data/data2014/source/NeonatalData.xlsx",
          sheet = "No TAM", method = "tab", header = TRUE, dec = ".")

# Drop unwanted columns
unwanted.cols <- c("no", "Gender","New.Number", "Gata1", "Preg.med", "Med.class",
           "X", "Aut.Neut", "Aut.mon", "Autbas", "Aut.eo", "MPV", "Day",
           "Preg.med.not.considered", "X.1")
xls <- xls[,!(names(xls) %in% unwanted.cols)]

data.nt <- na.omit(apply(xls, 2, as.numeric))

print(paste("Imported ", dim(data.nt)[1], " from ", dim(xls)[1], " records"))

############
# TAM data
############

xls <- read.xls(xls = "~/Documents/academic/infant_blood_project/data/data2014/source/NeonatalData.xlsx",
                sheet = "TAM", method = "tab", header = TRUE, dec = ".")

# Drop unwanted columns
unwanted.cols <- c("no", "Gender","New.Number", "Gata1", "Preg.med", "Med.class",
                   "X", "Aut.Neut", "Aut.mon", "Autbas", "MPV", "Day",
                   "Aut.eo","Med","Preg.med.not.considered", "X.1", "Clone.size.old")
xls <- xls[,!(names(xls) %in% unwanted.cols)]

data.t <- na.omit(apply(xls, 2, as.numeric))

print(paste("Imported ", dim(data.t)[1], " from ", dim(xls)[1], " records"))

data.nt <- cbind(rep(0, dim(data.nt)[1]), data.nt)
colnames(data.nt)[1] <- "Class"
data.t <- cbind(rep(1, dim(data.t)[1]), data.t)
colnames(data.t)[1] <- "Class"