# TME PHENOtyper

TME PHENOtyper is a random forest model built for discriminate immature and mature stroma subtypes from PDAC cancer samples. It has been built with pure stroma samples but it has been shown that the model prediction has an accuracy of 74,8% when samples have less than 50% of stroma content. 

## How to use TME PHENOtyper

In the first time, the packages, the data and the functions necessery have to be loaded
<pre><code>load("/Path/to/TME PHENOtyper/rf_fit.RData")
load("/Path/to/data/Data.RData")

PreProcess <- function(Data){
	  require(caret)
	  nzv <- nearZeroVar(Data)
	  if(length(nzv >0)){
	  	Data <- Data[, -nzv]
	  }

	  preProcValues <- preProcess(Data, method = c("center", "scale"))      
	  Data <- predict(preProcValues, Data)

	  return(Data)
  }</code></pre>

Data is a dataframe with samples as column, Ensembl ID of genes as rows and CPM as value. 
In the second time, data is preprocessed and TME PHENOtyper is used to determine subTMEs score for each samples

<pre><code>Data  <- PreProcess(t(Data))
subTMEs_Score <- predict(rf_fit, Data, type = "prob")
</code></pre>

The first column of subTMEs_Score, which is labelled 1, is the deserted subTME score and the second column is the reactive subTME score. The closer the score is 1, the higher the probability is. 


