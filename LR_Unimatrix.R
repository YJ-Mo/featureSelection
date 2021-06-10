library("broom")
library("caret")
library("doParallel")
library("dplyr")
library("edgR")
library("extrafont")
library("ggplot2")
library("ggpubr")
library("ggsci")
library("glmnet")
library("limma")
library("NMF")
library("optparse")
library("pROC")
library("reshape2")
library("tidyr")
option_list <- list( 
    make_option(c("-c","--count"), help="path of raw count matrix"),
    make_option(c("-l","--label"), help="path of label file"),
    make_option(c("-cl","--classifier"), help="path of classifier file"),
    make_option(c("-f","--fold"), help="number of cross validation fold",type="numeric",default="10"),
    make_option(c("-i","--iteration"), help="number of iterations",type="numeric",default="3"),
    make_option(c("-p","--partition"), help="number of partition times of discovery count matrix",type="numeric",default="100"),
    make_option(c("-p1","--plot_innerAUC"), help="path of saved AUC plot of inner CV",type="character",default="./innerAUC.png"),
    make_option(c("-p2","--plot_externalAUC"), help="path of saved AUC plot of external CV",type="character",default="./externalAUC.png")
)
opt <- parse_args(OptionParser(option_list=option_list))
source(opt$classifier)

### Step1: Input discovery count matrix file and label file
AllCount=read.table(opt$count,sep='\t',header=TRUE,row.names=1)
AllLabels=read.table(opt$label,sep='\t',header=TRUE)

DVSplit=mSplit(AllCount,AllLabels$label,1)
DiscoveryCount=AllCount[,DVSplit$samples$Resample]
Labels=AllLabels[DVSplit$samples$Resample,]
ValidationCount=AllCount[,-DVSplit$samples$Resample]
ValidationLabels=AllLabels[-DVSplit$samples$Resample,]
#ValidationCount=read.table(opt$validation,sep='\t',header=TRUE,row.names=1)

### Step2: Split discovery count matrix
Features.CVparam=trainControl(method="repeatedcv",number=opt$fold, repeats=opt$iteration,verboseIter=TRUE,returnData=FALSE,classProbs=TRUE,savePredictions=FALSE)
Cluster=makeCluster(5)
Cluster=registerDoParallel(Cluster)
DiscoveryIteration=list()
Splits=mSplit(DiscoveryCount,Labels$label,opt$partition)
for(i in 1:opt$partition) {
  DiscoveryIteration[[i]]=OnevsEach(DiscoveryCount,nDE=100,classes.df=Splits$df,Indices=Splits$samples[[i]]) 
  cat("Finish",i,"at",date())
}

### Step3: Cross validation on splited matrix
Classes.df <- Splits$df
TestPerformance.list <- list()
for(i in 1:100) {
    TestPerformance.list[[i]]=mPredict(ModelList=DiscoveryIteration[[i]],TestData=DiscoveryCount,Indices=Splits$samples[[i]], classes.df = Classes.df)}
AUCs.DiscoveryCohort <- GetAUC.ClassWise2(TestPerformance.list)
Counts.Samples <- count(Labels,label) %>% mutate(Frac = n/373)

### Step4: Ploting inner AUC
png(file=opt$plot_innerAUC)
Plot_innerAUC=qplot(data=AUCs.DiscoveryCohort,y=AUC,x=ID,geom="jitter",size=I(2), colour=I("#FF6666"))+
  geom_boxplot(colour=I("black"), alpha=I(0))+
  theme_bw()+
  ylab("AUCs within Discovery")+
  xlab("Cancer types")
Plot_innerAUC
while (!is.null(dev.list()))  dev.off()

### Step5: Predict on validation matrix
NC_AUCs=mValidate("NC")
NC_AUCs=data.frame(AUC=unlist(NC_AUCs), stringsAsFactors=False) %>% mutate(Class="NC")
HCC_AUCs=mValidate("HCC")
HCC_AUCs=data.frame(AUC=unlist(HCC_AUCs), stringsAsFactors=False) %>% mutate(Class="HCC")
CRC_AUCs=mValidate("CRC")
CRC_AUCs=data.frame(AUC=unlist(CRC_AUCs), stringsAsFactors=False) %>% mutate(Class="CRC")
STAD_AUCs=mValidate("STAD")
STAD_AUCs=data.frame(AUC=unlist(STAD_AUCs), stringsAsFactors=False) %>% mutate(Class="STAD")
ESCA_AUCs=mValidate("ESCA")
ESCA_AUCs=data.frame(AUC=unlist(ESCA_AUCs), stringsAsFactors=False) %>% mutate(Class="ESCA")
LUAD_AUCs=mValidate("LUAD")
LUAD_AUCs=data.frame(AUC=unlist(LUAD_AUCs), stringsAsFactors=False) %>% mutate(Class="LUAD")

### Step6: Plotting external AUC
All_AUCs=rbind(NC_AUCs, HCC_AUCs, CRC_AUCs, STAD_AUCs, ESCA_AUCs, LUAD_AUCs)
png(file=opt$plot_externalAUC)
Plot_externalAUC=qplot(data = All_AUCs, y = AUC, x = Class, geom = "jitter", size = I(2), colour = I("#46AFFF"))+
  geom_boxplot(colour = I("black"), alpha = I(0))+
  theme_bw()+
  ylab("AUCs within Validation")+
  ylim(c(0.4,1.0))
Plot_externalAUC
dev.off
