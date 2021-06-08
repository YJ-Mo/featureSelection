##Algorithmic design 

##For Each training set, create class-specific models
#This involves - DMR selection - one-vs-all
#GLMnet fitting by CV - binomial 
#Estimating performance on test sets 

SplitFunction <- function(Mat, Classes) {
  
  require(dplyr)
  require(caret)
  require(glmnet)
  
  ##Split into training and test sets
  
  df <- data.frame(ID = colnames(Mat), Classes = Classes)
  samples <- createDataPartition(df$Classes, p = 0.8,times = 100)
  
  return(list(df = df, samples = samples))
  
  
}


ModFunction <- function(Mat, classes.df, Indices) {
  
    
  #Then, we apply an iterative process of DMR selection and model training for each 
  #training set for each one vs all comparison
  
  
    
    TrainData <- Mat[,Indices]
    TrainPheno <- classes.df[Indices,]
    
    TestData <- Mat[,!colnames(Mat) %in% TrainPheno$ID]
    TestPheno <- classes.df%>%filter(!ID %in% TrainPheno$ID)
    
    
  AllClasses.v <- unique(classes.df$Classes)  
  ModList <- list()
    
  #Limma-trend preselection on the training partition 
  
  for(i in 1: length(AllClasses.v)) {
    
    
    NewAnn <- ifelse(TrainPheno$Classes == AllClasses.v[[i]], "One","Others") 
    Des <- model.matrix(~0 + NewAnn)
    colnames(Des) <- levels(factor(NewAnn))
    
    LimmaFit <- lmFit(TrainData,Des)%>%
      contrasts.fit(., makeContrasts(One-Others, levels = Des))%>%
      eBayes(., trend = TRUE)%>%
      topTable(., number = nrow(TrainData))
    
    LimmaFit <- LimmaFit%>%.[order(.$t),]
    TotalRows <- nrow(LimmaFit) - 49
    Features <- rbind(LimmaFit[1:50,] ,
                      LimmaFit[TotalRows:nrow(LimmaFit),])
    
    Features <- rownames(Features)
    message("DMR preselection complete")
    #GLMnet training et cetera; inherits CV parameters from a Features.CVparam object in the script.
    
    Model <- train(x = t(TrainData[rownames(TrainData) %in% Features,]), y = factor(NewAnn), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0,0.2,0.5,0.8,1),.lambda = seq(0,0.05,by=0.01)))
       message("Model Selection Complete")
   Prediction.classProbs <- predict(Model, newdata = t(TestData), type = "prob")%>%
     data.frame
   
   Prediction.classProbs$ActualClass <- TestPheno$Classes
   Prediction.classProbs$PredictedClass <- predict(Model, newdata = t(TestData), type = "raw")
   
   
   CombinedOutput <- list(Model = Model, TestPred = Prediction.classProbs)
   ModList[[i]] <- CombinedOutput
   
  }
    
    
    names(ModList) <- AllClasses.v 
    return(ModList)
    
    
  }

      
  
ModFunction.varyFeatureN <- function(Mat, classes.df, Indices, nDMR) {
  
  
  #Then, we apply an iterative process of DMR selection and model training for each 
  #training set for each one vs all comparison
  
  
  
  TrainData <- Mat[,Indices]
  TrainPheno <- classes.df[Indices,]
  
  TestData <- Mat[,!colnames(Mat) %in% TrainPheno$ID]
  TestPheno <- classes.df%>%filter(!ID %in% TrainPheno$ID)
  
  
  AllClasses.v <- unique(classes.df$Classes)  
  ModList <- list()
  
  #Limma-trend preselection on the training partition 
  
  for(i in 1: length(AllClasses.v)) {
    
    
    NewAnn <- ifelse(TrainPheno$Classes == AllClasses.v[[i]], "One","Others") 
    Des <- model.matrix(~0 + NewAnn)
    colnames(Des) <- levels(factor(NewAnn))
    
    LimmaFit <- lmFit(TrainData,Des)%>%
      contrasts.fit(., makeContrasts(One-Others, levels = Des))%>%
      eBayes(., trend = TRUE)%>%
      topTable(., number = nrow(TrainData))
    
    LimmaFit <- LimmaFit%>%.[order(.$t),]
    
    nDMR.b <- nDMR/2
    
    TotalRows <- nrow(LimmaFit) - (nDMR.b + 1)
    Features <- rbind(LimmaFit[1:nDMR.b,] ,
                      LimmaFit[TotalRows:nrow(LimmaFit),])
    
    Features <- rownames(Features)
    message(length(Features))
    
    message("DMR preselection complete")
    #GLMnet training et cetera; inherits CV parameters from a Features.CVparam object in the script.
    
    Model <- train(x = t(TrainData[rownames(TrainData) %in% Features,]), y = factor(NewAnn), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0,0.2,0.5,0.8,1),.lambda = seq(0,0.05,by=0.01)))
    message("Model Selection Complete")
    Prediction.classProbs <- predict(Model, newdata = t(TestData), type = "prob")%>%
      data.frame
    
    Prediction.classProbs$ActualClass <- TestPheno$Classes
    Prediction.classProbs$PredictedClass <- predict(Model, newdata = t(TestData), type = "raw")
    
    
    CombinedOutput <- list(Model = Model, TestPred = Prediction.classProbs)
    ModList[[i]] <- CombinedOutput
    
  }
  
  
  names(ModList) <- AllClasses.v 
  return(ModList)
  
  
}


GetAUC.ClassWise <- function(Runs) {
  
  
  Normals <- lapply(Runs, function(x) x$Normal)
  Normals.predictions <- lapply(Normals, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "Normal","One","Others")))
  Normals.auc <- lapply(Normals.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  PDACs <- lapply(Runs, function(x) x$PDAC)
  PDACs.predictions <- lapply(PDACs, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "PDAC","One","Others")))
  PDACs.auc <- lapply(PDACs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  BRCAs <- lapply(Runs, function(x) x$BRCA)
  BRCAs.predictions <- lapply(BRCAs, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "BRCA","One","Others")))
  BRCAs.auc <- lapply(BRCAs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  BLCAs <- lapply(Runs, function(x) x$BLCA)
  BLCAs.predictions <- lapply(BLCAs, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "BLCA","One","Others")))
  BLCAs.auc <- lapply(BLCAs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  LUCs <- lapply(Runs, function(x) x$LUC)
  LUCs.predictions <- lapply(LUCs, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "LUC","One","Others")))
  LUCs.auc <- lapply(LUCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  AMLs <- lapply(Runs, function(x) x$AML)
  AMLs.predictions <- lapply(AMLs, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "AML","One","Others")))
  AMLs.auc <- lapply(AMLs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  CRCs <- lapply(Runs, function(x) x$CRC)
  CRCs.predictions <- lapply(CRCs, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "CRC","One","Others")))
  CRCs.auc <- lapply(CRCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  RCCs <- lapply(Runs, function(x) x$RCC)
  RCCs.predictions <- lapply(RCCs, function(x) x$TestPred%>%mutate(Class2 = ifelse(ActualClass == "RCC","One","Others")))
  RCCs.auc <- lapply(RCCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  #Aggregate into boxplots
  
  AMLs.auc <- data.frame(AUC = unlist(AMLs.auc))%>%
    mutate(ID = "AML")
  
  BLCAs.auc <- data.frame(AUC = unlist(BLCAs.auc))%>%
    mutate(ID = "BLCA")
  
  BRCAs.auc <- data.frame(AUC = unlist(BRCAs.auc))%>%
    mutate(ID = "BRCA")
  
  CRCs.auc <- data.frame(AUC = unlist(CRCs.auc))%>%
    mutate(ID = "CRC")
  
  LUCs.auc <- data.frame(AUC = unlist(LUCs.auc))%>%
    mutate(ID = "LUC")
  
  Normals.auc <- data.frame(AUC = unlist(Normals.auc))%>%
    mutate(ID = "Normal")
  
  PDACs.auc <- data.frame(AUC = unlist(PDACs.auc))%>%
    mutate(ID = "PDAC")
  
  RCCs.auc <- data.frame(AUC = unlist(RCCs.auc))%>%
    mutate(ID = "RCC")
  
  Bound <- rbind(AMLs.auc,BLCAs.auc,BRCAs.auc,CRCs.auc,LUCs.auc,Normals.auc,PDACs.auc,RCCs.auc)
  
  
  return(Bound)
  
  
}


#Select DMRs in one vs each class procedure and fit a GLMnet 

OnevsEach <- function(Mat, classes.df, Indices, nDMR) {
  
  
  #Training/test-set split 
  
  TrainData <- Mat[,Indices]
  TrainPheno <- classes.df[Indices,]
  
  TestData <- Mat[,!colnames(Mat) %in% TrainPheno$ID]
  TestPheno <- classes.df%>%filter(!ID %in% TrainPheno$ID)
  
  #Set up models list
  AllClasses.v <- unique(classes.df$Classes)  
  ModList <- list()
  
  for(i in 1:length(AllClasses.v)) {
    
    FixedClass <- which(TrainPheno$Classes == AllClasses.v[[i]])
    OtherClasses <- which(!TrainPheno$Classes == AllClasses.v[[i]])
    
    
      DMRList <- list()
      OtherClasses.vector <- unique(TrainPheno$Classes[OtherClasses])
      print(OtherClasses.vector)
    ##This loop does DMR preselection using a one vs each criterion
    for(j in 1:length(OtherClasses.vector)) {
      
      CurrentOtherClass <- which(TrainPheno$Classes == OtherClasses.vector[[j]] )
      
      FixedClass.matrix <- TrainData[,FixedClass]
      OtherMatrix <- TrainData[,CurrentOtherClass]
      DMR.classes <- c(rep("One",ncol(FixedClass.matrix)), rep("Others",ncol(OtherMatrix)))
      
      DMR.Data <- cbind(FixedClass.matrix, OtherMatrix)
      Des <- model.matrix(~0 + DMR.classes)
      colnames(Des) <- levels(factor(DMR.classes))
      
      LimmaFit <- lmFit(DMR.Data, Des)%>%
        contrasts.fit(., makeContrasts(One-Others, levels = Des))%>%
        eBayes(., trend = TRUE)%>%
        topTable(., number = nrow(FixedClass.matrix))

      LimmaFit <- LimmaFit%>%.[order(.$t),]
      
      nDMR.b <- nDMR/2
      
      TotalRows <- nrow(LimmaFit) - (nDMR.b + 1)
      Features <- rbind(LimmaFit[1:nDMR.b,] ,
                        LimmaFit[TotalRows:nrow(LimmaFit),])
      
      Features <- rownames(Features)
      DMRList[[j]] <- Features      
      message(paste0(j,"of each vs other classes DMR selection done"))
    }
    
    #This creates feature set 
    Features <- unlist(DMRList)
    
    #Here we fit the model and chuck it into the modlist 
    
    NewAnn <- ifelse(TrainPheno$Classes == AllClasses.v[[i]],"One","Others")
    
    Model <- train(x = t(TrainData[rownames(TrainData) %in% Features,]), y = factor(NewAnn), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0,0.2,0.5,0.8,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")
    message("Model Selection Complete")
    Prediction.classProbs <- predict(Model, newdata = t(TestData), type = "prob")%>%
      data.frame
    
    Prediction.classProbs$ActualClass <- TestPheno$Classes
    Prediction.classProbs$PredictedClass <- predict(Model, newdata = t(TestData), type = "raw")
    
    
    CombinedOutput <- list(Model = Model, TestPred = Prediction.classProbs)
    ModList[[i]] <- CombinedOutput
    
  }
  
  names(ModList) <- AllClasses.v  
  return(ModList)  
    
}
  
  
#Get AUC.classwise2

GetAUC.ClassWise2 <- function(Runs) {
  
  
  Normals <- lapply(Runs, function(x) x$Normal)
  Normals.predictions <- lapply(Normals, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "Normal","One","Others")))
  Normals.auc <- lapply(Normals.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  PDACs <- lapply(Runs, function(x) x$PDAC)
  PDACs.predictions <- lapply(PDACs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "PDAC","One","Others")))
  PDACs.auc <- lapply(PDACs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  BRCAs <- lapply(Runs, function(x) x$BRCA)
  BRCAs.predictions <- lapply(BRCAs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "BRCA","One","Others")))
  BRCAs.auc <- lapply(BRCAs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  BLCAs <- lapply(Runs, function(x) x$BLCA)
  BLCAs.predictions <- lapply(BLCAs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "BLCA","One","Others")))
  BLCAs.auc <- lapply(BLCAs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  LUCs <- lapply(Runs, function(x) x$LUC)
  LUCs.predictions <- lapply(LUCs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "LUC","One","Others")))
  LUCs.auc <- lapply(LUCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  AMLs <- lapply(Runs, function(x) x$AML)
  AMLs.predictions <- lapply(AMLs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "AML","One","Others")))
  AMLs.auc <- lapply(AMLs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  CRCs <- lapply(Runs, function(x) x$CRC)
  CRCs.predictions <- lapply(CRCs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "CRC","One","Others")))
  CRCs.auc <- lapply(CRCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  RCCs <- lapply(Runs, function(x) x$RCC)
  RCCs.predictions <- lapply(RCCs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "RCC","One","Others")))
  RCCs.auc <- lapply(RCCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  #Aggregate into boxplots
  
  AMLs.auc <- data.frame(AUC = unlist(AMLs.auc))%>%
    mutate(ID = "AML")
  
  BLCAs.auc <- data.frame(AUC = unlist(BLCAs.auc))%>%
    mutate(ID = "BLCA")
  
  BRCAs.auc <- data.frame(AUC = unlist(BRCAs.auc))%>%
    mutate(ID = "BRCA")
  
  CRCs.auc <- data.frame(AUC = unlist(CRCs.auc))%>%
    mutate(ID = "CRC")
  
  LUCs.auc <- data.frame(AUC = unlist(LUCs.auc))%>%
    mutate(ID = "LUC")
  
  Normals.auc <- data.frame(AUC = unlist(Normals.auc))%>%
    mutate(ID = "Normal")
  
  PDACs.auc <- data.frame(AUC = unlist(PDACs.auc))%>%
    mutate(ID = "PDAC")
  
  RCCs.auc <- data.frame(AUC = unlist(RCCs.auc))%>%
    mutate(ID = "RCC")
  
  Bound <- rbind(AMLs.auc,BLCAs.auc,BRCAs.auc,CRCs.auc,LUCs.auc,Normals.auc,PDACs.auc,RCCs.auc)
  
  
  return(Bound)
  
  
}


# Get AUPR classwise 


GetAUPR.ClassWise <- function(Runs) {
  
  require(PRROC)
  
  Normals <- lapply(Runs, function(x) x$Normal)
  Normals.predictions <-lapply(Normals, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "Normal",1,0))%>%mutate(Class2 = as.numeric(Class2)))
 Normals.auc <- lapply(Normals.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral)
  
  PDACs <- lapply(Runs, function(x) x$PDAC)
  PDACs.predictions <- lapply(PDACs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "PDAC",1,0))%>%mutate(Class2 = as.numeric(Class2)))
  
  PDACs.auc <- lapply(PDACs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral)
 
  
   BRCAs <- lapply(Runs, function(x) x$BRCA)
  BRCAs.predictions <-   lapply(BRCAs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "BRCA",1,0))%>%mutate(Class2 = as.numeric(Class2)))
  BRCAs.auc <- lapply(BRCAs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral)
  
   
  BLCAs <- lapply(Runs, function(x) x$BLCA)
  BLCAs.predictions <- lapply(BLCAs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "BLCA",1,0))%>%mutate(Class2 = as.numeric(Class2)))
  BLCAs.auc <- lapply(BLCAs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral)

  
    LUCs <- lapply(Runs, function(x) x$LUC)
  LUCs.predictions <- lapply(LUCs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "LUC",1,0))%>%mutate(Class2 = as.numeric(Class2)))
  LUCs.auc <- lapply(LUCs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral) 

  
  AMLs <- lapply(Runs, function(x) x$AML)
  AMLs.predictions <- lapply(AMLs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "AML",1,0))%>%mutate(Class2 = as.numeric(Class2)))
  AMLs.auc <- lapply(AMLs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral)
  
  CRCs <- lapply(Runs, function(x) x$CRC)
  CRCs.predictions <- lapply(CRCs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "CRC",1,0))%>%mutate(Class2 = as.numeric(Class2)))
  CRCs.auc <- lapply(CRCs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral)
  
  RCCs <- lapply(Runs, function(x) x$RCC)
  RCCs.predictions <- lapply(RCCs, function(x) x%>%mutate(Class2 = ifelse(ActualClass == "RCC",1,0))%>%mutate(Class2 = as.numeric(Class2)))
  RCCs.auc <- lapply(RCCs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.integral)
  
  
#Compute AUCPR by integration here. 
  
  AMLs.auc <- data.frame(AUC = unlist(AMLs.auc))%>%
    mutate(ID = "AML")
  
  BLCAs.auc <- data.frame(AUC = unlist(BLCAs.auc))%>%
    mutate(ID = "BLCA")
  
  BRCAs.auc <- data.frame(AUC = unlist(BRCAs.auc))%>%
    mutate(ID = "BRCA")
  
  CRCs.auc <- data.frame(AUC = unlist(CRCs.auc))%>%
    mutate(ID = "CRC")
  
  LUCs.auc <- data.frame(AUC = unlist(LUCs.auc))%>%
    mutate(ID = "LUC")
  
  Normals.auc <- data.frame(AUC = unlist(Normals.auc))%>%
    mutate(ID = "Normal")
  
  PDACs.auc <- data.frame(AUC = unlist(PDACs.auc))%>%
    mutate(ID = "PDAC")
  
  RCCs.auc <- data.frame(AUC = unlist(RCCs.auc))%>%
    mutate(ID = "RCC")
  
  
  Bound.integral <- rbind(AMLs.auc,BLCAs.auc,BRCAs.auc,CRCs.auc,LUCs.auc,Normals.auc,PDACs.auc,RCCs.auc)
  
  #Then calculate using the goadrich davis method 
  
  
  Normals.auc <- lapply(Normals.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich)
  PDACs.auc <- lapply(PDACs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich)
  BRCAs.auc <- lapply(BRCAs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich)
  BLCAs.auc <- lapply(BLCAs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich)
  LUCs.auc <- lapply(LUCs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich) 
  AMLs.auc <- lapply(AMLs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich)
  CRCs.auc <- lapply(CRCs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich)
  RCCs.auc <- lapply(RCCs.predictions, function(x) pr.curve(scores.class0 = x$One, weights.class0 = x$Class2)$auc.davis.goadrich)
  
  Bound.davis.goadrich <- rbind(AMLs.auc,BLCAs.auc,BRCAs.auc,CRCs.auc,LUCs.auc,Normals.auc,PDACs.auc,RCCs.auc)
  
  
  return(list(Integral = Bound.integral, DavisGoadrich = Bound.davis.goadrich))
  
}
