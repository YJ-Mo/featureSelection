library('RUVSeq')
library('stringr')
library("optparse")
option_list <- list( 
		    make_option("--count", help="path of raw count matrix"),
		    make_option("--label", help = "path of label file"),
        make_option("--ref", help = "path of ensembl references file"),
		    make_option("--tmm", help="ouput file of TMM_normalized matrix"),
		    make_option("--ruv",help="ouput file of TMM + RUVg_normalized matrix"),
		    make_option("--anova",help="output of anova table")
	           )
opt <- parse_args(OptionParser(option_list=option_list))

################################################
####  Step_1: Input count matrix and label  ####
################################################
count_matrix=read.csv(opt$count,sep="\t",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE,row.names=1)
labels=read.csv(opt$label,sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
ref=read.csv(opt$ref,sep = '\t',header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
sampleIDs=colnames(count_matrix)
featureIDs=rownames(count_matrix)
cat("Finish input at",date())

################################################
####    Step_2: Cross check sample names    ####
################################################
labels=subset(labels,labels$sample_id %in% sampleIDs)
sampleIDs=labels$sample_id
uni_label=unique(labels$label)
class=factor(labels$label)
group=as.factor(labels$label)
count_matrix=count_matrix[,sampleIDs]
cat("Finish cross check names at",date())

################################################
####     Step_3: Cross check gene names     ####
################################################
#fea_vec=c()
#for (fea in featureIDs){
#  if (TRUE %in% str_detect(fea, ref$`Gene stable ID`)){ fea_vec=append(fea_vec,fea) }}
# count_matrix=subset(count_matrix, rownames(count_matrix) %in% fea_vec)
cat("Finish cross check genes at",date())

################################################
####         Step_4: Filter sparsity        ####
################################################
for (i in 1: length(uni_label)){
  filter <- apply(count_matrix[,grep(uni_label[i], colnames(count_matrix))], 1, function(x) length(x[x>0])/length(x)>=0.2)
  count_matrix <- count_matrix[filter,]}
cat("Finish filtering sparsity at",date())

################################################
####     Step_5: Differential Analysis      ####
################################################
design=model.matrix(~group)
y=DGEList(counts=count_matrix, group=group)
y=calcNormFactors(y, method="TMM")
y=estimateDisp(y, design)
count_matrix.tmm=cpm(y)
fit=glmFit(y, design)
lrt=glmLRT(fit, coef=2:6)
top=topTags(lrt, n=nrow(count_matrix))$table
n.mask=as.integer(0.8*nrow(as.matrix(count_matrix)))
empirical=rownames(count_matrix)[which(!(rownames(count_matrix) %in% rownames(top)[1:n.mask]))]
cat("Finish DE at",date())

################################################
####          Step_6: Remove Batch          ####
################################################
res=RUVg(x=log(count_matrix.tmm+0.1), cIdx=empirical, k=2,isLog = T)
count_matrix.tmm.ruv=exp(res$normalizedCounts)
cat("Finish removing batch at",date())

################################################
####             Step_7: Output             ####
################################################
write.table(count_matrix.tmm,opt$tmm,quote=F,sep="\t")
write.table(count_matrix.tmm.ruv,opt$ruv,quote=F,sep="\t")
write.table(top,opt&anova,quote=F,sep="\t")
cat("Finish output at",date())
cat("Finish all at",date())
