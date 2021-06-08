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
#install.packages('e1071', dependencies=TRUE)

#source the functions and load the data objects here
source("~/Desktop/LuLab/cancer/MLFinal/OnevAllClassifiers.R")
load("~/Desktop/LuLab/cancer/ProcessedDataArchive_ML/StartingPoints.RData")
Combined <- log2(Combined * 0.3 + 1e-6)
load("~/Desktop/LuLab/cancer/ProcessedDataArchive_ML/OnevAll_Splits.RData")

Features.CVparam<- trainControl(method = "repeatedcv",number = 10, repeats =3,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=FALSE)

library(doParallel)
Cluster <-makeCluster(5)
Cluster <- registerDoParallel(Cluster)

AllIterations.onevEach <- list()
Splits<-SplitFunction(Combined,Input$SampGroups) # added
for(i in 1:100) {
  AllIterations.onevEach[[i]] <- OnevsEach(Combined, classes.df = Splits$df, Indices = Splits$samples[[i]], nDMR = 300) 
  message(i)
}

save(AllIterations.onevEach, file = "AllIterations_300Features_OnevEach.RData")
load("~/Desktop/LuLab/cancer/ProcessedDataArchive_ML/AllIterations_300Features_OnevEach.RData") #added
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

#I can plot the random expected estimate too based on the prevalence. 

Counts.Samples <- count(Input, SampGroups)%>%
  mutate(Frac = n/189)


Plot1 <- qplot(data = AUCs.DiscoveryCohort, y = AUC, x = ID, geom = "jitter", size = I(2), colour = I("orange"))+
  geom_boxplot(colour = I("black"), alpha = I(0))+
  theme_bw()+
  ylab("Discovery Cohort AUCs")+
  xlab("Class")
Plot1

# on validation
load("~/Desktop/LuLab/cancer/ProcessedDataArchive_ML/Parsed_ValidationData.RData")

#AML Validation

AMLs <- c(lapply(AllIterations.onevEach, function(x) x$AML))

AMLs <- lapply(AMLs, function(x) x$Model)

AMLModels.Validation <- list() 

for(i in 1:length(AMLs)) {
  
  Features <- AMLs[[i]]$finalModel$xNames
  ValData <- Validation.mat[match(Features, rownames(Validation.mat)),]
  Predictions <- predict(AMLs[[i]], newdata = t(ValData), type ="prob")
  AMLModels.Validation[[i]] <- Predictions
  message(i) 
  
}

AMLModels.Validation <- lapply(AMLModels.Validation, function(x) x%>%as.data.frame%>%mutate(Class = ifelse(Validation.Pheno$Group == "AML",1,0)))

AMLs.AUCs <- lapply(AMLModels.Validation, function(x) with(x,auc(Class ~ One)))

rm(AMLs)


AMLs <- c(lapply(AllIterations.onevEach, function(x) x$AML))

AMLs <- lapply(AMLs, function(x) x$Model)

#PDAC

PDACs <- c(lapply(AllIterations.onevEach, function(x) x$PDAC))

PDACs <- lapply(PDACs, function(x) x$Model)

PDACModels.Validation <- list() 

for(i in 1:length(PDACs)) {
  
  Features <- PDACs[[i]]$finalModel$xNames
  ValData <- Validation.mat[match(Features, rownames(Validation.mat)),]
  Predictions <- predict(PDACs[[i]], newdata = t(ValData), type ="prob")
  PDACModels.Validation[[i]] <- Predictions
  message(i) 
  
}

PDACModels.Validation <- lapply(PDACModels.Validation, function(x) x%>%as.data.frame%>%mutate(Class = ifelse(Validation.Pheno$Group == "PDAC",1,0)))

PDACs.AUCs <- lapply(PDACModels.Validation, function(x) with(x,auc(Class ~ One)))


rm(PDACs)


#Normal


Normals <- c(lapply(AllIterations.onevEach, function(x) x$Normal))

Normals <- lapply(Normals, function(x) x$Model)

NormalModels.Validation <- list() 

for(i in 1:length(Normals)) {
  
  Features <- Normals[[i]]$finalModel$xNames
  ValData <- Validation.mat[match(Features, rownames(Validation.mat)),]
  Predictions <- predict(Normals[[i]], newdata = t(ValData), type ="prob")
  NormalModels.Validation[[i]] <- Predictions
  message(i) 
  
}

NormalModels.Validation <- lapply(NormalModels.Validation, function(x) x%>%as.data.frame%>%mutate(Class = ifelse(Validation.Pheno$Group == "Normal",1,0)))

Normals.AUCs <- lapply(NormalModels.Validation, function(x) with(x,auc(Class ~ One)))
rm(Normals)

Normals <- c(lapply(AllIterations.onevEach, function(x) x$Normal))

Normals <- lapply(Normals, function(x) x$Model)

#Lung

LUCs <- c(lapply(AllIterations.onevEach, function(x) x$LUC))

LUCs <- lapply(LUCs, function(x) x$Model)

LUCModels.Validation <- list() 

for(i in 1:length(LUCs)) {
  
  Features <- LUCs[[i]]$finalModel$xNames
  ValData <- Validation.mat[match(Features, rownames(Validation.mat)),]
  Predictions <- predict(LUCs[[i]], newdata = t(ValData), type ="prob")
  LUCModels.Validation[[i]] <- Predictions
  message(i) 
  
}

LUCModels.Validation <- lapply(LUCModels.Validation, function(x) x%>%as.data.frame%>%mutate(Class = ifelse(Validation.Pheno$Group == "LUC",1,0)))

LUCs.AUCs <- lapply(LUCModels.Validation, function(x) with(x,auc(Class ~ One)))
rm(LUCs)

#Bind and plot

AMLs.AUCs <- data.frame(AUC = unlist(AMLs.AUCs), stringsAsFactors = F)%>%
  mutate(Class = "AML")

PDACs.AUCs <- data.frame(AUC = unlist(PDACs.AUCs), stringsAsFactors = F)%>%
  mutate(Class = "PDAC")

Normals.AUCs <- data.frame(AUC = unlist(Normals.AUCs), stringsAsFactors = F)%>%
  mutate(Class = "Normal")

LUCs.AUCs <- data.frame(AUC = unlist(LUCs.AUCs), stringsAsFactors = F)%>%
  mutate(Class = "LUC")

AUCTab <- rbind(AMLs.AUCs, PDACs.AUCs, Normals.AUCs, LUCs.AUCs)
Plot2 <- qplot(data = AUCTab, y = AUC, x = Class, geom = "jitter", size = I(2), colour = I("orange"))+
  geom_boxplot(colour = I("black"), alpha = I(0))+
  theme_bw()+
  ylab("Validation Cohort AUCs")+
  ylim(c(0.4,1.0))
Plot2

library(gridExtra)

#pdf(width = 10, height = 5, file = "OnevsEachFeatures.pdf")
grid.arrange(Plot1, Plot2, nrow = 1, widths = c(6,3))

# roc curve
Class.LUC <- LUCModels.Validation[[1]]$Class
ClassProbs.LUC <- lapply(LUCModels.Validation, function(x) x$One)
ClassProbs.LUC <- do.call(cbind, ClassProbs.LUC)
ClassProbs.LUC <- rowSums(ClassProbs.LUC)/100
ClassProbs.LUC <- data.frame(Probability = ClassProbs.LUC, Classes =Class.LUC, stringsAsFactors = F)


Class.AML <- AMLModels.Validation[[1]]$Class
ClassProbs.AML <- lapply(AMLModels.Validation, function(x) x$One)
ClassProbs.AML <- do.call(cbind, ClassProbs.AML)
ClassProbs.AML <- rowSums(ClassProbs.AML)/100
ClassProbs.AML <- data.frame(Probability = ClassProbs.AML, Classes =Class.AML, stringsAsFactors = F)


Class.Normal <- NormalModels.Validation[[1]]$Class
ClassProbs.Normal <- lapply(NormalModels.Validation, function(x) x$One)
ClassProbs.Normal <- do.call(cbind, ClassProbs.Normal)
ClassProbs.Normal <- rowSums(ClassProbs.Normal)/100
ClassProbs.Normal <- data.frame(Probability = ClassProbs.Normal, Classes =Class.Normal, stringsAsFactors = F)



Class.PDAC <- PDACModels.Validation[[1]]$Class
ClassProbs.PDAC <- lapply(PDACModels.Validation, function(x) x$One)
ClassProbs.PDAC <- do.call(cbind, ClassProbs.PDAC)
ClassProbs.PDAC <- rowSums(ClassProbs.PDAC)/100
ClassProbs.PDAC <- data.frame(Probability = ClassProbs.PDAC, Classes =Class.PDAC, stringsAsFactors = F)


#pdf(width = 8, height = 8, file = "ROCCurves_ModelAverages.pdf")

par(mfrow = c(2,2))
with(ClassProbs.LUC, plot(roc(Classes ~ Probability), main = "LUC", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))

with(ClassProbs.AML, plot(roc(Classes ~ Probability), main = "AML", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))

with(ClassProbs.PDAC, plot(roc(Classes ~ Probability), main = "PDAC", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))

with(ClassProbs.Normal, plot(roc(Classes ~ Probability), main = "Normal", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))

### 画stage相关的图
ClassProbs.PDAC$Sample <- Validation.Pheno$Sample
ClassProbs.AML$Sample <- Validation.Pheno$Sample
ClassProbs.Normal$Sample <- Validation.Pheno$Sample
ClassProbs.LUC$Sample <- Validation.Pheno$Sample

#Demarcate early and late stage ; I will need this
ES.lung <- filter(Validation.Pheno, Stage2 == "ES" & Group == "LUC")
LS.lung <- filter(Validation.Pheno, Stage2 == "LS" & Group == "LUC")
ES.PDAC <- filter(Validation.Pheno, Stage2 == "ES" & Group == "PDAC")
LS.PDAC <- filter(Validation.Pheno, Stage2 == "LS" & Group == "PDAC")


#Calculate LUC ES vs others, and LUC LS vs other AUROC and AUPCR

Lung.EarlyStage <- filter(ClassProbs.LUC, !Sample %in% LS.lung$Sample)
Lung.LateStage <-  filter(ClassProbs.LUC, ! Sample %in% ES.lung$Sample)
PDAC.EarlyStage <- filter(ClassProbs.PDAC, !Sample %in% LS.PDAC$Sample)
PDAC.LateStage <-  filter(ClassProbs.PDAC, ! Sample %in% ES.PDAC$Sample)


#pdf(width = 10, height = 2.5, file = "ROCs_StageStratified.pdf")
par(mfrow = c(1,4))
with(Lung.EarlyStage, plot(roc(Classes ~ Probability), main = "LUC - Early Stage", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))

with(Lung.LateStage, plot(roc(Classes ~ Probability), main = "LUC - Late Stage", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))

with(PDAC.EarlyStage, plot(roc(Classes ~ Probability), main = "PDAC - Early Stage", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))

with(PDAC.LateStage, plot(roc(Classes ~ Probability), main = "PDAC - Late Stage", col = "dodgerblue3",smooth = TRUE, legacy = T, print.auc = T, ylim = c(0,1)))


#### 获取非0特征
#Load Models again 
AllIterations.onevEach <- do.call(c, AllIterations.onevEach)
Nonzeros <- lapply(AllIterations.onevEach, function(x) coef(x$Model$finalModel, s = x$Model$bestTune$lambda))

Nonzeros <- lapply(Nonzeros, function(x) as.matrix(x))
Nonzeros <- lapply(Nonzeros, function(x) data.frame(ID = rownames(x), Coef = as.numeric(x[,1])))
Nonzeros <- lapply(Nonzeros, function(x) filter(x, !Coef == 0))
Nonzeros <- do.call(rbind, Nonzeros)

#Get the union of windows here 

Union <- unique(Nonzeros$ID)
Union2 <- count(Nonzeros, ID)
save(Union, Nonzeros, Union2, file = "UnionDMRs.RData")
