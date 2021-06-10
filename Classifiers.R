mPredict=function(ModelList,TestData,Indices,classes.df) { 
  TrainPheno=classes.df[Indices,]
  TestData=TestData[,!colnames(TestData) %in% TrainPheno$ID]
  TestPheno=classes.df%>%filter(!ID %in% TrainPheno$ID)
  Predictions.list=list()
  OutputNames=names(ModelList)
  for(i in 1:length(ModelList)){
    Features=ModelList[[i]]$Model$finalModel$xNames
    TestDataNew=TestData[match(Features,rownames(TestData)),]
    Model=ModelList[[i]]$Model
    Prediction.classProbs=predict(Model, newdata=t(TestDataNew), type="prob")%>%
      data.frame
    Prediction.classProbs$ActualClass=TestPheno$Classes
    Prediction.classProbs$PredictedClass=predict(Model, newdata=t(TestDataNew), type="raw")
    Predictions.list[[i]]=Prediction.classProbs}
  names(Predictions.list)=OutputNames
  return(Predictions.list)
} 

mValidate=function(type) {
Type=c(lapply(DiscoveryIteration, function(x) x[[type]]))
Type=lapply(Type, function(x) x$Model)
Models.Validation=list() 
for(i in 1:length(Type)) {
  Features <- Type[[i]]$finalModel$xNames
  ValData <- ValidationCount[match(Features, rownames(ValidationCount)),]
  Predictions <- predict(Type[[i]], newdata = t(ValData), type ="prob")
  Models.Validation[[i]] <- Predictions
  message(i)
}
Models.Validation <- lapply(Models.Validation, function(x) x%>%as.data.frame%>%mutate(Class = ifelse(ValidationLabel$label == [[type]],1,0)))
Type.AUCs <- lapply(Models.Validation, function(x) with(x,auc(Class ~ One)))
rm(Type)
return(Type.AUCs)
}

mValidate2=function(type) {
Type=c(lapply(DiscoveryIteration, function(x) x$type))
Type=lapply(Type, function(x) x$Model)
Models.Validation=list() 
for(i in 1:length(Type)) {
  Features <- Type[[i]]$finalModel$xNames
  ValData <- Validation.mat[match(Features, rownames(Validation.mat)),]
  Predictions <- predict(Type[[i]], newdata = t(ValData), type ="prob")
  Models.Validation[[i]] <- Predictions
  message(i)
}
Models.Validation <- lapply(Models.Validation, function(x) x%>%as.data.frame%>%mutate(Class = ifelse(Validation.Pheno$Group == type,1,0)))
Type.AUCs <- lapply(Models.Validation, function(x) with(x,auc(Class ~ One)))
rm(Type)
return(Type.AUCs)
}
mSplit=function(Count,Label,Num){
  require(dplyr)
  require(caret)
  require(glmnet)
  df=data.frame(ID=colnames(Count),Classes=Label)
  splited=createDataPartition(df$Classes,p=0.8,times=Num)
  return(list(df=df, samples=splited))
}


mMod=function(Mat,classes.df,Indices){
    TrainData=Mat[,Indices]
    TrainPheno=classes.df[Indices,]
    TestData=Mat[,!colnames(Mat) %in% TrainPheno$ID]
    TestPheno=classes.df%>%filter(!ID %in% TrainPheno$ID)
    AllClasses.v=unique(classes.df$Classes)  
    ModList=list()
    for(i in 1: length(AllClasses.v)) {
    NewAnn=ifelse(TrainPheno$Classes == AllClasses.v[[i]], "One","Others") 
    Des=model.matrix(~0 + NewAnn)
    colnames(Des)=levels(factor(NewAnn))
    LimmaFit=lmFit(TrainData,Des)%>%
      contrasts.fit(., makeContrasts(One-Others, levels=Des))%>%
      eBayes(., trend=TRUE)%>%
      topTable(., number=nrow(TrainData))
    LimmaFit=LimmaFit%>%.[order(.$t),]
    TotalRows=nrow(LimmaFit) - 49
    Features=rbind(LimmaFit[1:50,] ,
                      LimmaFit[TotalRows:nrow(LimmaFit),])
    Features=rownames(Features)
    Model=train(x=t(TrainData[rownames(TrainData) %in% Features,]), y=factor(NewAnn), trControl=Features.CVparam, method="glmnet" , tuneGrid=expand.grid(.alpha=c(0,0.2,0.5,0.8,1),.lambda=seq(0,0.05,by=0.01)))
   Prediction.classProbs=predict(Model, newdata=t(TestData), type="prob") %>% data.frame
   Prediction.classProbs$ActualClass=TestPheno$Classes
   Prediction.classProbs$PredictedClass=predict(Model, newdata=t(TestData), type="raw")
   CombinedOutput=list(Model=Model, TestPred=Prediction.classProbs)
   ModList[[i]]=CombinedOutput
 }
    names(ModList)=AllClasses.v 
    return(ModList)
}

mMod.varyFeatureN=function(Mat, classes.df, Indices, nDE) {
  TrainData=Mat[,Indices]
  TrainPheno=classes.df[Indices,]
  TestData=Mat[,!colnames(Mat) %in% TrainPheno$ID]
  TestPheno=classes.df%>%filter(!ID %in% TrainPheno$ID)
  AllClasses.v=unique(classes.df$Classes)  
  ModList=list()
  for(i in 1: length(AllClasses.v)) {
    NewAnn=ifelse(TrainPheno$Classes == AllClasses.v[[i]], "One","Others") 
    Des=model.matrix(~0 + NewAnn)
    colnames(Des)=levels(factor(NewAnn))
    LimmaFit=lmFit(TrainData,Des)%>%
      contrasts.fit(., makeContrasts(One-Others, levels=Des))%>%
      eBayes(., trend=TRUE)%>%
      topTable(., number=nrow(TrainData))
    LimmaFit=LimmaFit%>%.[order(.$t),]
    nDE.b=nDE/2
    TotalRows=nrow(LimmaFit) - (nDE.b + 1)
    Features=rbind(LimmaFit[1:nDE.b,] ,
                      LimmaFit[TotalRows:nrow(LimmaFit),])
    Features=rownames(Features)
    Model=train(x=t(TrainData[rownames(TrainData) %in% Features,]), y=factor(NewAnn), trControl=Features.CVparam, method="glmnet" , tuneGrid=expand.grid(.alpha=c(0,0.2,0.5,0.8,1),.lambda=seq(0,0.05,by=0.01)))
    message("Model Selection Complete")
    Prediction.classProbs=predict(Model, newdata=t(TestData), type="prob") %>% data.frame
    Prediction.classProbs$ActualClass=TestPheno$Classes
    Prediction.classProbs$PredictedClass=predict(Model, newdata=t(TestData), type="raw")
    CombinedOutput=list(Model=Model, TestPred=Prediction.classProbs)
    ModList[[i]]=CombinedOutput
 }
  names(ModList)=AllClasses.v 
  return(ModList) 
}

GetAUC.ClassWise=function(Runs) {
  NCs=lapply(Runs, function(x) x$NC)
  NCs.predictions=lapply(NCs, function(x) x$TestPred%>%mutate(Class2=ifelse(ActualClass == "NC","One","Others")))
  NCs.auc=lapply(NCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  ESCAs=lapply(Runs, function(x) x$ESCA)
  ESCAs.predictions=lapply(ESCAs, function(x) x$TestPred%>%mutate(Class2=ifelse(ActualClass == "ESCA","One","Others")))
  ESCAs.auc=lapply(ESCAs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  HCCs=lapply(Runs, function(x) x$HCC)
  HCCs.predictions=lapply(HCCs, function(x) x$TestPred%>%mutate(Class2=ifelse(ActualClass == "HCC","One","Others")))
  HCCs.auc=lapply(HCCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  LUADs=lapply(Runs, function(x) x$LUAD)
  LUADs.predictions=lapply(LUADs, function(x) x$TestPred%>%mutate(Class2=ifelse(ActualClass == "LUAD","One","Others")))
  LUADs.auc=lapply(LUADs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  STADs=lapply(Runs, function(x) x$STAD)
  STADs.predictions=lapply(STADs, function(x) x$TestPred%>%mutate(Class2=ifelse(ActualClass == "STAD","One","Others")))
  STADs.auc=lapply(STADs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  CRCs=lapply(Runs, function(x) x$CRC)
  CRCs.predictions=lapply(CRCs, function(x) x$TestPred%>%mutate(Class2=ifelse(ActualClass == "CRC","One","Others")))
  CRCs.auc=lapply(CRCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  LUADs.auc=data.frame(AUC=unlist(LUADs.auc))%>%
    mutate(ID="LUAD")
  
  HCCs.auc=data.frame(AUC=unlist(HCCs.auc))%>%
    mutate(ID="HCC")
  
  CRCs.auc=data.frame(AUC=unlist(CRCs.auc))%>%
    mutate(ID="CRC")
  
  STADs.auc=data.frame(AUC=unlist(STADs.auc))%>%
    mutate(ID="STAD")
  
  NCs.auc=data.frame(AUC=unlist(NCs.auc))%>%
    mutate(ID="NC")
  
  ESCAs.auc=data.frame(AUC=unlist(ESCAs.auc))%>%
    mutate(ID="ESCA")
                  
  Bound=rbind(LUADs.auc,HCCs.auc,CRCs.auc,STADs.auc,NCs.auc,ESCAs.auc)
  return(Bound)
}

OnevsEach=function(Mat, classes.df, Indices, nDE){
  TrainData=Mat[,Indices]
  TrainPheno=classes.df[Indices,]
  TestData=Mat[,!colnames(Mat) %in% TrainPheno$ID]
  TestPheno=classes.df%>%filter(!ID %in% TrainPheno$ID)
  AllClasses.v=unique(classes.df$Classes)  
  ModList=list()
  for(i in 1:length(AllClasses.v)) {
    FixedClass=which(TrainPheno$Classes == AllClasses.v[[i]])
    OtherClasses=which(!TrainPheno$Classes == AllClasses.v[[i]])
      DEList=list()
      OtherClasses.vector=unique(TrainPheno$Classes[OtherClasses])
      print(OtherClasses.vector)
    for(j in 1:length(OtherClasses.vector)) {
      CurrentOtherClass=which(TrainPheno$Classes == OtherClasses.vector[[j]] )
      FixedClass.matrix=TrainData[,FixedClass]
      OtherMatrix=TrainData[,CurrentOtherClass]
      DE.classes=c(rep("One",ncol(FixedClass.matrix)), rep("Others",ncol(OtherMatrix)))
      DE.Data=cbind(FixedClass.matrix, OtherMatrix)
      Des=model.matrix(~0 + DE.classes)
      colnames(Des)=levels(factor(DE.classes))
      LimmaFit=lmFit(DE.Data, Des)%>%
        contrasts.fit(., makeContrasts(One-Others, levels=Des))%>%
        eBayes(., trend=TRUE)%>%
        topTable(., number=nrow(FixedClass.matrix))
      LimmaFit=LimmaFit%>%.[order(.$t),]
      nDE.b=nDE/2
      TotalRows=nrow(LimmaFit) - (nDE.b + 1)
      Features=rbind(LimmaFit[1:nDE.b,] ,
                        LimmaFit[TotalRows:nrow(LimmaFit),])
      Features=rownames(Features)
      DEList[[j]]=Features      
      message(paste0(j,"of each vs other classes DE selection done"))
    }
    Features=unlist(DEList)
    NewAnn=ifelse(TrainPheno$Classes == AllClasses.v[[i]],"One","Others")
    Model=train(x=t(TrainData[rownames(TrainData) %in% Features,]), y=factor(NewAnn), trControl=Features.CVparam, method="glmnet" , tuneGrid=expand.grid(.alpha=c(0,0.2,0.5,0.8,1),.lambda=seq(0,0.05,by=0.01)), metric="Kappa")
    message("Model Selection Complete")
    Prediction.classProbs=predict(Model, newdata=t(TestData), type="prob") %>% data.frame
    Prediction.classProbs$ActualClass=TestPheno$Classes
    Prediction.classProbs$PredictedClass=predict(Model, newdata=t(TestData), type="raw")
    
    
    CombinedOutput=list(Model=Model, TestPred=Prediction.classProbs)
    ModList[[i]]=CombinedOutput
  
  }
  names(ModList)=AllClasses.v  
  return(ModList)  
}
  


GetAUC.ClassWise2=function(Runs) {
  
  NCs=lapply(Runs, function(x) x$NC)
  NCs.predictions=lapply(NCs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "NC","One","Others")))
  NCs.auc=lapply(NCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  ESCAs=lapply(Runs, function(x) x$ESCA)
  ESCAs.predictions=lapply(ESCAs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "ESCA","One","Others")))
  ESCAs.auc=lapply(ESCAs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  HCCs=lapply(Runs, function(x) x$HCC)
  HCCs.predictions=lapply(HCCs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "HCC","One","Others")))
  HCCs.auc=lapply(HCCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  LUADs=lapply(Runs, function(x) x$LUAD)
  LUADs.predictions=lapply(LUADs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "LUAD","One","Others")))
  LUADs.auc=lapply(LUADs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  STADs=lapply(Runs, function(x) x$STAD)
  STADs.predictions=lapply(STADs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "STAD","One","Others")))
  STADs.auc=lapply(STADs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  CRCs=lapply(Runs, function(x) x$CRC)
  CRCs.predictions=lapply(CRCs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "CRC","One","Others")))
  CRCs.auc=lapply(CRCs.predictions, function(x) with(x,roc(Class2 ~ One)$auc))
  
  LUADs.auc=data.frame(AUC=unlist(LUADs.auc))%>% mutate(ID="LUAD")
  
  HCCs.auc=data.frame(AUC=unlist(HCCs.auc))%>% mutate(ID="HCC")
  
  CRCs.auc=data.frame(AUC=unlist(CRCs.auc)) %>% mutate(ID="CRC")
  
  STADs.auc=data.frame(AUC=unlist(STADs.auc)) %>% mutate(ID="STAD")
  
  NCs.auc=data.frame(AUC=unlist(NCs.auc)) %>% mutate(ID="NC")
  
  ESCAs.auc=data.frame(AUC=unlist(ESCAs.auc)) %>% mutate(ID="ESCA")
  
  Bound=rbind(LUADs.auc,HCCs.auc,CRCs.auc,STADs.auc,NCs.auc,ESCAs.auc)
  return(Bound)
}

GetAUPR.ClassWise=function(Runs) {
  require(PRROC)
    
  NCs=lapply(Runs, function(x) x$NC)
  NCs.predictions <-lapply(NCs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "NC",1,0))%>%mutate(Class2=as.numeric(Class2)))
  NCs.auc=lapply(NCs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.integral)
  
  ESCAs=lapply(Runs, function(x) x$ESCA)
  ESCAs.predictions=lapply(ESCAs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "ESCA",1,0))%>%mutate(Class2=as.numeric(Class2)))
  ESCAs.auc=lapply(ESCAs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.integral)
 
  HCCs=lapply(Runs, function(x) x$HCC)
  HCCs.predictions=  lapply(HCCs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "HCC",1,0))%>%mutate(Class2=as.numeric(Class2)))
  HCCs.auc=lapply(HCCs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.integral)
   
  LUADs=lapply(Runs, function(x) x$LUAD)
  LUADs.predictions=lapply(LUADs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "LUAD",1,0))%>%mutate(Class2=as.numeric(Class2)))
  LUADs.auc=lapply(LUADs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.integral)
  
  STADs=lapply(Runs, function(x) x$STAD)
  STADs.predictions=lapply(STADs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "STAD",1,0))%>%mutate(Class2=as.numeric(Class2)))
  STADs.auc=lapply(STADs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.integral) 

  CRCs=lapply(Runs, function(x) x$CRC)
  CRCs.predictions=lapply(CRCs, function(x) x%>%mutate(Class2=ifelse(ActualClass == "CRC",1,0))%>%mutate(Class2=as.numeric(Class2)))
  CRCs.auc=lapply(CRCs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.integral)
  
  LUADs.auc=data.frame(AUC=unlist(LUADs.auc)) %>% mutate(ID="LUAD")
  HCCs.auc=data.frame(AUC=unlist(HCCs.auc)) %>% mutate(ID="HCC")
  CRCs.auc=data.frame(AUC=unlist(CRCs.auc)) %>% mutate(ID="CRC")
  STADs.auc=data.frame(AUC=unlist(STADs.auc)) %>% mutate(ID="STAD")
  NCs.auc=data.frame(AUC=unlist(NCs.auc)) %>% mutate(ID="NC")
  ESCAs.auc=data.frame(AUC=unlist(ESCAs.auc)) %>% mutate(ID="ESCA")
  
  Bound.integral=rbind(LUADs.auc,HCCs.auc,CRCs.auc,STADs.auc,NCs.auc,ESCAs.auc)
  NCs.auc=lapply(NCs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.davis.goadrich)
  ESCAs.auc=lapply(ESCAs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.davis.goadrich)
  HCCs.auc=lapply(HCCs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.davis.goadrich)
  LUADs.auc=lapply(LUADs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.davis.goadrich)
  STADs.auc=lapply(STADs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.davis.goadrich) 
  CRCs.auc=lapply(CRCs.predictions, function(x) pr.curve(scores.class0=x$One, weights.class0=x$Class2)$auc.davis.goadrich)

  Bound.davis.goadrich=rbind(LUADs.auc,HCCs.auc,CRCs.auc,STADs.auc,NCs.auc,ESCAs.auc)
  return(list(Integral=Bound.integral, DavisGoadrich=Bound.davis.goadrich))
}
