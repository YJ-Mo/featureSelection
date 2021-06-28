library(caret)
library(dplyr)
library(glmnet)
library(edgeR)
library(pROC)
library(ggplot2)
counts=read.table("/Users/yajin/Desktop/LuLab/multiomics/counts.txt",header = T,row.names = 1)
labels=read.table("/Users/yajin/Desktop/LuLab/multiomics/labels.txt",header = T)
df=data.frame(ID=colnames(counts),Classes=labels)

mPredict=funciton(counts,labels){
  prediction=data.frame()
for (i in 1:ncol(counts)){
train=counts[,-i]
train_label=labels[-i,]
test=as.data.frame(counts[,i])
rownames(test)=rownames(counts)
  
dgelist <- DGEList(counts = train, group = train_label$label)
group = train_label$label

keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
design <- model.matrix(~group)

dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
lrt=lrt[order(lrt$table$logFC),]
Features=lrt$table[1:200,] #取FC的top200
Features=rownames(Features)

Features.CVparam=trainControl(method="repeatedcv",number=1, repeats=1,verboseIter=TRUE,returnData=FALSE,classProbs=TRUE,savePredictions=FALSE)
Model=train(x=t(train[rownames(train) %in% Features,]), y=factor(train_label$label),  method="glmnet" , tuneGrid=expand.grid(.alpha=c(0,0.2,0.5,0.8,1),.lambda=seq(0,0.05,by=0.01)), metric="Kappa")

Prediction.classProbs=predict(Model, newdata=t(test), type="prob") %>% data.frame
Prediction.classProbs$ActualClass=labels$label[i]
Prediction.classProbs$PredictedClass=predict(Model, newdata=t(test), type="raw")
prediction=rbind(prediction,Prediction.classProbs)
message("finish",i)
}
  return(prediction)
}
