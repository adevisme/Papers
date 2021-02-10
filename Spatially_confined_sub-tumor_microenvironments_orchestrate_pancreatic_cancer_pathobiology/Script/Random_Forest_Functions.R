ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

entrez2name <- function(entrez)
{
  gname <- mget(as.character(entrez), org.Hs.egGENENAME, ifnotfound=NA)
  gname <- unlist(lapply(gname, function(i) return(i[1])))
  return(gname)
}

RandomForest <- function(Data, Comparaison, Nb=3, Repeats=10, PathSave, SetSeed=147896325){
require(ggplot2)
require(ROCR)
require(caret)


	set.seed(SetSeed)
	
	Data.train <- Data
	Data.ctrl <- Data
	
	Stroma <- Data.train[,Comparaison]
	
	##########################################
	### First run with all the genes
	########################################## 
	
	PathSave.tmp <- paste0(PathSave, "/Before_threshold")
	          
	### fit_control
	
	fit_control <- trainControl(method = "LOOCV", 
                        	savePredictions = "all",
                            returnResamp = "all",
                            search = 'grid')

	
	### rf_fit

	rf_fit <- train(as.factor(Stroma) ~ ., 
              data =  Data.train, 
              method = "rf",
              metric="Accuracy",
              trControl = fit_control, 
              na.action=na.omit)


	### Important variable 

	rf_Imp <- varImp(rf_fit, scale = FALSE)         
         
         
	Data_rf_pred <- predict(rf_fit, Data.ctrl)
	Con.Matrix <- confusionMatrix(Data_rf_pred, as.factor(Data.ctrl[,Comparaison]))

	setwd(PathSave.tmp)
	save(Data.train, file="Data_train.RData")
	save(Data.ctrl, file="Data_ctrl.RData")
	save(rf_fit, file="rf_fit.RData")
	save(Con.Matrix, file="Con.Matrix.RData")
	save(rf_Imp, file="rf_Imp.RData")


	####################################################
	### Second run with only the significant features
	####################################################
	
	PathSave.tmp <- paste0(PathSave, "/After_threshold")
		
	Data.tmp <- t(Data.train)
	Data.tmp <- Data.tmp[-dim(Data.tmp)[1],] # Remove Stroma row
		
	rf_Imp.matrix <- as.matrix(rf_Imp$importance)
				
	# Filter rf_imp
	rf_Imp.matrix <- rf_Imp.matrix[order(rf_Imp.matrix, decreasing = TRUE),]
	rf_Imp.matrix.filter <- rf_Imp.matrix[rf_Imp.matrix>0.030]
		
	#Save Features
	require(org.Hs.eg.db)
	Matrix.Save <- as.data.frame(rf_Imp.matrix.filter)
	colnames(Matrix.Save) <- "Genes_Impact"
	EntrezID <- ensembl2entrez(rownames(Matrix.Save))
	Symbol <- entrez2symbol(EntrezID)
	Name <- entrez2name(EntrezID)
	Matrix.Save <- cbind(Matrix.Save, EntrezID)
	Matrix.Save <- cbind(Matrix.Save, Symbol)
	Matrix.Save <- cbind(Matrix.Save, Name)
	write.xlsx(Matrix.Save, file=paste0(PathSave.tmp, "/Genes_used.xlsx"), row.names = TRUE, firstRow=TRUE)	
		
	#Filter Features
	Data.tmp <- Data.tmp[match(names(rf_Imp.matrix.filter), rownames(Data.tmp)),]
	Data.train <- t(Data.tmp)
	Data.train <- cbind(Data.train, Stroma)
			
			
	### fit_control
	

	fit_control <- trainControl(method = "LOOCV", 
    	                    savePredictions = "all",
    	                    returnResamp = "all",
    	                    search = 'grid')

		
	### rf_fit
		
	rf_fit <- train(as.factor(Stroma) ~ ., 
   	          data =  Data.train, 
   	          method = "rf",
   	          metric="Accuracy",
   	          trControl = fit_control, 
   	          na.action=na.omit)
	
		
	### Important variable 

	rf_Imp <- varImp(rf_fit, scale = FALSE)         
         
         
	Data_rf_pred <- predict(rf_fit, Data.ctrl)
	Con.Matrix <- confusionMatrix(Data_rf_pred, as.factor(Data.ctrl[,Comparaison]))
	
	setwd(PathSave.tmp)

	save(Data.train, file="Data_train.RData")
	save(Data.ctrl, file="Data_ctrl.RData")
	save(rf_fit, file="rf_fit.RData")
	save(Con.Matrix, file="Con.Matrix.RData")
	save(rf_Imp, file="rf_Imp.RData")
	
}

