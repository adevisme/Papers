library(caret)
library(openxlsx)
library(randomForest)
library(mlbench)
library(ROCR)


set.seed(147896325)

source("~/Random_Forest_Function.R")

Data.Ann <- read.xlsx("~/Annotations_Files.xlsx"
Data <- read.xlsx("~/Data_Random_Forest.xlsx"

### pre processing
Data.rf <- t(Data)
Data.Ann <- Data.Ann[grep("T", Data.Ann[,"Type"]),] # T is the annotation for Stroma samples  
Stroma <-  rep_len(2, length(Data.Ann[,"Stroma"]))
Stroma[grep("^mature", Data.Ann[,"Stroma"])] <- 1



### remove samples with a variance close to zero

nzv <- nearZeroVar(Data.rf)
if(length(nzv >0)){
	Data.rf <- Data.rf[, -nzv]
}

### centering and scaling

preProcValues <- preProcess(Data.rf, method = c("center", "scale"))      
Data.rf <- predict(preProcValues, Data.rf)

#############

Data.rf <- cbind(Data.rf, Stroma)

RandomForest(Data=Data.rf, Comparaison="Stroma", Split=FALSE, Nb=3, Repeats=10, PathSave="Path/to/save", SetSeed=147896325)

