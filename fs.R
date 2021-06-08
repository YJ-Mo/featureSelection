library("broom")
library("caret")
library("doParallel")
library("dplyr")
library("ggplot2")
library("glmnet")
library("limma")
library("NMF")
library("optparse")
library("pROC")
library("reshape2")
library("tidyr")
library("optparse")
option_list <- list( 
    make_option("--count", help="path of raw count matrix"),
    make_option("--label", help = "path of label file"),
    make_option("--classifier", help = "path of classifier file"),
    make_option("--fold", help = "number of cross validation fold"),
    make_option("--iteration", help = "number of iterations"),
    make_option("--partition", help = "number of partition times of discovery count matrix"),
    make_option("--plot_innerAUC", help = "path of saved AUC plot of inner CV")
)
opt <- parse_args(OptionParser(option_list=option_list))
source(opt$classifier)

### Step1: Input discovery count matrix file and label file
DiscoveryCount=read.table(opt$count,sep='\t',header=TRUE,row.names=1)
Labels=read.table(opt$label,sep='\t',header=TRUE)

### Step2: Split discovery count matrix
Features.CVparam=trainControl(method="repeatedcv",number=opt$fold, repeats=opt$iteration,verboseIter=TRUE,returnData=FALSE,classProbs=TRUE,savePredictions=FALSE)
Cluster=makeCluster(5)
Cluster=registerDoParallel(Cluster)
DiscoveryIteration=list()
Splits=mSplit(DiscoveryCount,Labels$label,opt$partition)
for(i in 1:opt$partition) {
  DiscoveryIteration[[i]]=OnevsEach(DiscoveryCount,nDE=200,classes.df=Splits$df,Indices=Splits$samples[[i]]) 
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
dev.off
