library('RUVSeq')
library('stringr')
input_matrix='your_dir/pico.txt'
input_label='your_dir/pico_label.txt'
input_ref='your_dir/human_ref.txt'
tmm='your_dir/pico_tmm.txt'
ruv='your_dir/pico_ruv.txt'
anova='your_dir/pico_anova.txt'

################################################
####  Step_1: Input count matrix and label  ####
################################################
count_matrix=read.csv(input_matrix,sep="\t",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE,row.names=1)
labels=read.csv(input_label,sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
ref=read.csv(input_ref,sep = '\t',header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
sampleIDs=colnames(count_matrix)
featureIDs=rownames(count_matrix)

################################################
####    Step_2: Cross check sample names    ####
################################################
labels=subset(labels,labels$sample_id %in% sampleIDs)
sampleIDs=labels$sample_id
uni_label=unique(labels$label)
class=factor(labels$label)
group=as.factor(labels$label)
count_matrix=count_matrix[,sampleIDs]

################################################
####     Step_3: Cross check gene names     ####
################################################
fea_vec=c()
for (fea in featureIDs){
  if (TRUE %in% str_detect(fea, ref$`Gene stable ID`)){ fea_vec=append(fea_vec,fea) }}

count_matrix=subset(count_matrix, rownames(count_matrix) %in% fea_vec)

################################################
####         Step_4: Filter sparsity        ####
################################################
#for (i in 1: length(uni_label)){
#  filter <- apply(count_matrix[,grep(uni_label[i], colnames(count_matrix))], 1, function(x) length(x[x>0])/length(x)>=0.2)
#  count_matrix <- count_matrix[filter,]}

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

################################################
####          Step_6: Remove Batch          ####
################################################
res=RUVg(x=log(count_matrix.tmm+0.1), cIdx=empirical, k=2,isLog = T)
count_matrix.tmm.ruv=exp(res$normalizedCounts)

################################################
####             Step_7: Output             ####
################################################
write.table(count_matrix.tmm,tmm,quote=F,sep="\t")
write.table(count_matrix.tmm.ruv,ruv,quote=F,sep="\t")
write.table(top,anova,quote=F,sep="\t")
