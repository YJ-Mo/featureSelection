library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(broom)
library(caret)
library(limma)
library(glmnet)
library(NMF)
library(doParallel)
library(pROC)

source("your_dir/Classifiers.R")
Combined=read.table("your_dir/pico.txt",sep='\t',header=TRUE,row.names=1)
Input=read.table("your_dir//pico_label.txt",sep='\t',header=TRUE)
Features.CVparam<- trainControl(method = "repeatedcv",number = 8, repeats =2,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=FALSE)

Cluster <-makeCluster(5)
Cluster <- registerDoParallel(Cluster)
AllIterations.onevEach <- list()
Splits<-SplitFunction(Combined,Input$label)
for(i in 1:100) {
  AllIterations.onevEach[[i]] <- OnevsEach(Combined, classes.df = Splits$df, Indices = Splits$samples[[i]], nDMR = 300) 
  message(i)
}

PredFunction <- function(ModelList, TestData, Indices, classes.df) { 
  TrainPheno <- classes.df[Indices,]
  TestData <- TestData[,!colnames(TestData) %in% TrainPheno$ID]
  TestPheno <- classes.df%>%filter(!ID %in% TrainPheno$ID)
  
  Predictions.list <- list()
  OutputNames <- names(ModelList)
  
  for(i in 1:length(ModelList)) {
    Features <- ModelList[[i]]$Model$finalModel$xNames
    TestDataNew <- TestData[match(Features,rownames(TestData)),]
    Model <- ModelList[[i]]$Model
    Prediction.classProbs <- predict(Model, newdata = t(TestDataNew), type = "prob")%>%
      data.frame
    Prediction.classProbs$ActualClass <- TestPheno$Classes
    Prediction.classProbs$PredictedClass <- predict(Model, newdata = t(TestDataNew), type = "raw")
    Predictions.list[[i]] <- Prediction.classProbs
    message(i)
  }
  
  names(Predictions.list) <- OutputNames
  return(Predictions.list)
} 


Classes.df <- Splits$df
TestPerformance.list <- list()
for(i in 1:100) {
  TestPerformance.list[[i]] <- PredFunction(ModelList = AllIterations.onevEach[[i]],
                                            TestData = Combined, Indices = Splits$samples[[i]], classes.df = Classes.df)
}
AUCs.DiscoveryCohort <- GetAUC.ClassWise2(TestPerformance.list)
Counts.Samples <- count(Input, SampGroups)%>%
  mutate(Frac = n/189)

Plot1 <- qplot(data = AUCs.DiscoveryCohort, y = AUC, x = ID, geom = "jitter", size = I(2), colour = I("orange"))+
  geom_boxplot(colour = I("black"), alpha = I(0))+
  theme_bw()+
  ylab("Discovery Cohort AUCs")+
  xlab("Class")
Plot1
