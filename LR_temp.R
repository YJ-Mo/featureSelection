library("broom")
library("caret")
library("doParallel")
library("dplyr")
library("edgeR")
library("extrafont")
library("ggplot2")
library("ggpubr")
library("ggsci")
library("ggsignif")
library("glmnet")
library("limma")
library("magrittr")
library("NMF")
library("optparse")
library("pROC")
library("reshape2")
library("tidyr")
option_list <- list(
    make_option(c("-c","--count"), help="path of raw count matrix",type="character"),
    make_option(c("-l","--label"), help="path of label file",type="character"),
    make_option(c("--validationCount"), help="path of raw count matrix for independent validation",type="character",default="stupid"),
    make_option(c("--validationLabel"), help="path of label file for validation count",type="character",default="stupid"),
    make_option(c("--classifier"), help="path of classifier file"),
    make_option(c("-f","--fold"), help="number of cross validation fold",type="numeric",default="10"),
    make_option(c("-i","--iteration"), help="number of iterations",type="numeric",default="3"),
    make_option(c("-p","--partition"), help="number of partition times of discovery count matrix",type="numeric",default="100"),
    make_option(c("-o","--output"), help="path to output AUC plots and data frames",type="character",default="."),
    make_option(c("-d","--degene"), help="number of DE genes",type="numeric",default="200")
)
opt <- parse_args(OptionParser(option_list=option_list))
source(opt$classifier)

### Step1: Input discovery count matrix file and label file
if (opt$validationCount != "stupid"){
cat("Bimode")
DiscoveryCount=read.table(opt$count,sep='\t',header=TRUE,row.names=1)
Labels=read.table(opt$label,sep='\t',header=TRUE)
ValidationCount=read.table(opt$validationCount,sep='\t',header=TRUE,row.names=1)
ValidationLabels=read.table(opt$validationLabel,sep='\t',header=TRUE)
} else {
cat("Unimode")
AllCount=read.table(opt$count,sep='\t',header=TRUE,row.names=1)
AllLabels=read.table(opt$label,sep='\t',header=TRUE)
DVSplit=mSplit(AllCount,AllLabels$label,1)
DiscoveryCount=AllCount[,DVSplit$samples$Resample]
Labels=AllLabels[DVSplit$samples$Resample,]
ValidationCount=AllCount[,-DVSplit$samples$Resample]
ValidationLabels=AllLabels[-DVSplit$samples$Resample,]
}

### Step2: Split discovery count matrix
Features.CVparam=trainControl(method="repeatedcv",number=opt$fold, repeats=opt$iteration,verboseIter=TRUE,returnData=FALSE,classProbs=TRUE,savePredictions=FALSE)
Cluster=makeCluster(5)
Cluster=registerDoParallel(Cluster)
DiscoveryIteration=list()
Splits=mSplit(DiscoveryCount,Labels$label,opt$partition)
for(i in 1:opt$partition) {
  DiscoveryIteration[[i]]=OnevsEach(DiscoveryCount,nDE=opt$degene,classes.df=Splits$df,Indices=Splits$samples[[i]])
  cat("Finish",i,"at",date())
}

### Step3: Cross validation on splited matrix
Classes.df <- Splits$df
TestPerformance.list <- list()
for(i in 1:opt$partition) {
    TestPerformance.list[[i]]=mPredict(ModelList=DiscoveryIteration[[i]],TestData=DiscoveryCount,Indices=Splits$samples[[i]], classes.df = Classes.df)}
AUCs.Discovery <- GetAUC.ClassWise2(TestPerformance.list)
write.table(AUCs.Discovery,paste(opt$output,"/AUC_Discovery.txt",sep=""),sep="\t",quote=FALSE)
Counts.Samples <- count(Labels,label) %>% mutate(Frac = n/373)

### Step4: Ploting inner AUC
png(file=paste(opt$output,"/inner_AUC.png",sep=""))
plot_innerAUC=ggboxplot(data=AUCs.Discovery,x='ID',y='AUC',fill='ID',bxp.errorbar=T,bxp.errorbar.width=0.2,palette="npg")+
  labs(titles='Inner AUC',
       subtitle='AUCs within DIscovery Cross Validation',
       caption='Data on gene set1',
       x='Cancer Type',
       y='AUC'
       )+
  theme(
    plot.title = element_text(colour = 'black', size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size=16,hjust = 0.5),
    plot.caption = element_text(color = 'black', size = 16, face = 'italic', hjust = 1),
    legend.title = element_text(color = 'black', size = 16),
    legend.text = element_text(color = 'black', size = 16),
    axis.line.x = element_line(colour = 'black', linetype = 'solid'),
    axis.line.y = element_line(colour = 'black', linetype = 'solid'),
    panel.border = element_rect(linetype = 'solid', size = 1.2, fill = NA)
)
plot_innerAUC
while (!is.null(dev.list()))  dev.off()

### Step5: Predict on validation matrix
NC_AUCs=mValidate("NC")
NC_AUCs=data.frame(AUC=unlist(NC_AUCs), stringsAsFactors=FALSE) %>% mutate(Class="NC")
HCC_AUCs=mValidate("HCC")
HCC_AUCs=data.frame(AUC=unlist(HCC_AUCs), stringsAsFactors=FALSE) %>% mutate(Class="HCC")
CRC_AUCs=mValidate("CRC")
CRC_AUCs=data.frame(AUC=unlist(CRC_AUCs), stringsAsFactors=FALSE) %>% mutate(Class="CRC")
STAD_AUCs=mValidate("STAD")
STAD_AUCs=data.frame(AUC=unlist(STAD_AUCs), stringsAsFactors=FALSE) %>% mutate(Class="STAD")
ESCA_AUCs=mValidate("ESCA")
ESCA_AUCs=data.frame(AUC=unlist(ESCA_AUCs), stringsAsFactors=FALSE) %>% mutate(Class="ESCA")
LUAD_AUCs=mValidate("LUAD")
LUAD_AUCs=data.frame(AUC=unlist(LUAD_AUCs), stringsAsFactors=FALSE) %>% mutate(Class="LUAD")

### Step6: Plotting external AUC
All_AUCs=rbind(NC_AUCs, HCC_AUCs, CRC_AUCs) #STAD_AUCs, ESCA_AUCs, LUAD_AUCs)
write.table(All_AUCs,paste(opt$output,"/AUC_Validation.txt",sep=""),sep="\t",quote=FALSE)
png(file=paste(opt$output,"/external_AUC.png",sep=""))
#Plot_externalAUC=mPlot(All_AUCs)
plot_externalAUC=ggboxplot(data=All_AUCs,x='Class',y='AUC',fill='Class',bxp.errorbar=T,bxp.errorbar.width=0.2,palette="npg")+
  labs(titles='Inner AUC',
       subtitle='AUCs within DIscovery Cross Validation',
       caption='Data on gene set1',
       x='Cancer Type',
       y='AUC'
       )+
  theme(
    plot.title = element_text(colour = 'black', size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size=16,hjust = 0.5),
    plot.caption = element_text(color = 'black', size = 16, face = 'italic', hjust = 1),
    legend.title = element_text(color = 'black', size = 16),
    legend.text = element_text(color = 'black', size = 16),
    axis.line.x = element_line(colour = 'black', linetype = 'solid'),
    axis.line.y = element_line(colour = 'black', linetype = 'solid'),
    panel.border = element_rect(linetype = 'solid', size = 1.2, fill = NA)
)
plot_externalAUC
while (!is.null(dev.list()))  dev.off()
