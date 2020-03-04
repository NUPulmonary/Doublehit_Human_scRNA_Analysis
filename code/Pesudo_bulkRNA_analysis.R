#AM integrated pesudo bulk RNA-seq analysis
#6 dataset HPVMMR
library(edgeR)
setwd("/Users/zrt781/Desktop/doublehit_paper/scRNA_seq_analysis")
x<-read.csv("AM_HPVMMR_pesudo_bulkRNA.csv",header = T)
x<-x[!duplicated(x$X),]
row.names(x)<-x$X
x<-x[,2:53]
targets<-read.csv("AM_HPVMMR_pesudo_bulkRNA_targets.csv",header=T)
head(x)
targets$SampleName<-colnames(x)

targets2<-as.data.frame(cbind(targets$SampleName,as.character(targets$AgeGroup)))
names(targets2)<-c("SampleName","condition")
x<-round(x)
color<-c("Age55+"="red","Age<55"="blue")
getpca_rlog_matrix(x,targets2,color,"PCA_results_pesudo_AMbulk_020420")
library(ggplot2)


targets3<-as.data.frame(cbind(targets$SampleName,as.character(targets$Study)))
names(targets3)<-c("SampleName","condition")
color<-c("Reyfman"="orange","Morse"="blue","Valenzi"="red","Madissoon"="black","Raredon"="purple","Habermann"="pink")
getpca_rlog_matrix(x,targets3,color,"PCA_results_pesudo_AMbulk_study_020420")
library(ggplot2)

#######
#model V3
#######
x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_counts.csv",header=T,row.names=1)
#x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_SCT_counts.csv",header=T,row.names=1)
#x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_RNA_counts.csv",header=T,row.names=1)
#x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_counts.csv",header=T,row.names=1)


colnames(x3)
cellnames<-read.csv("colnames_6datasets.csv",header=T)
colnames(x3)<-colnames(cellnames)

targets<-read.csv("AM_targets_quest_heatmap.csv",header=T)
targets$Sex[targets$Sex=="Male"]<-"male"


targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(x3)==targets$SampleName)


remove.names<-c(grep("_24h",colnames(x3)),grep("_12h",colnames(x3)),grep("_72h",colnames(x3)))
x3<-x3[,-remove.names]
targets <- targets[-remove.names,-1]
sum(colnames(x3)==targets$SampleName)


row.names(targets2)<-NULL
targets3<-targets3[match(colnames(x3),targets3$SampleName),]
sum(colnames(x3)==targets3$SampleName)
row.names(targets3)<-NULL
#remove negative value;
x3[x3<0]<-0
#x3<-round(x3)
sum(x3<0)
#x3<-round(log2(x3+1))

color<-c("Age55+"="red","Age<55"="blue")
getpca_rlog_matrix(x3,targets2,color,"PCA_results_pesudo_AMbulk_6dataRNA_020520")


color2<-c("Reyfman"="orange","Morse"="blue","Valenzi"="red","Madissoon"="black","Raredon"="purple","Habermann"="pink")
getpca_rlog_matrix(x3,targets3,color2,"PCA_results_pesudo_AMbulk_study_6dataRNA_020520")
library(ggplot2)
sum(x3)

#run DE analysis on the integrated 
#targets$AgeGroup2[targets$Age<=30]<-"Age30less"
#targets$AgeGroup2[targets$Age>60]<-"Age60more"
#targets4<-targets[!is.na(targets$AgeGroup2),]
#x3<-x3[,colnames(x3) %in% targets4$SampleName]
#targets4<-targets4[,c(1,6)]
#colnames(targets4)<-c("SampleName","condition")
#targets4<-targets4[match(colnames(x3),targets4$SampleName),]
#sum(colnames(x3)==targets4$SampleName)

targets4<-targets[!is.na(targets$AgeGroup2),]
targets4<-targets4[,c(1,6)]
colnames(targets4)<-c("SampleName","condition")
x3<-x3[,colnames(x3) %in% targets4$SampleName]
sum(colnames(x3)==targets4$SampleName)

#group <- factor(paste0(targets4$condition)) #define groups based on targets
group <- factor(paste0(targets4$condition)) #define groups based on targets

samplenames <- colnames(x3)
y <- DGEList(counts=x3,group=group) #Define DGEList and design matrix
y <- calcNormFactors(y) #calculate normalization factors

#Filtering
#keep <- rowSums(cpm(y)>1) >= 6
#y <- y[keep, , keep.lib.sizes=FALSE]
#Recompute library size
y$samples$lib.size <- colSums(y$counts)
#Vizualize library size
barplot(y$samples$lib.size*1e-6, names=1:27, ylab="Library size (millions)", xlab="sample")
abline(h=4, col="red", lty=2, lwd=2)
#TMM normalization 
y <- calcNormFactors(y)
y$samples
#Vizualize effect of TMM normalization (inspect all samples)
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
#EstimatingDispersion, GLM
#To estimate common dispersion and tagwise dispersions in one run (recommended):
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmFit(y, design, robust=TRUE) #Fit a negative binomial generalized log-linear model to counts

#AMD vs. AMR
SvsM <- makeContrasts(Age60more - Age30less, levels=design)
lrt <- glmLRT(fit, contrast=SvsM)

is.de <- decideTestsDGE(lrt)
summary(is.de) #check how many sig genes


dataH<-as.data.frame(topTags(lrt,n=14000))
write.csv(dataH,file="Pesudo_bulk_AM_DEgene60vs30_list_filtered.csv")



#HC heatmap
logCPM <- cpm(y, prior.count = 2, log=TRUE)
#write.csv(logCPM,file="AM_logCPM_integrated_aging3000.csv")
rownames(logCPM) <- rownames(y$counts)
colnames(logCPM) <- targets4$SampleName
#colnames(logCPM) <- targets3$SampleName

o <- order(lrt$table$PValue)
logCPM <- logCPM[o[1:784],]


logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")

library(pheatmap)
sum(colnames(logCPM)==targets4$SampleName)
targets5<-targets[!is.na(targets$AgeGroup2),]
sum(targets5$SampleName==colnames(logCPM))

targets5<-targets5[,c(4:6)]
row.names(targets5)<-targets4$SampleName

pheatmap(logCPM,annotation_col = targets5,show_rownames=F,main="All DEgeen (784 genes) between 60+vs30- AM Human")

heatmap.2(logCPM, col=col.pan, Rowv=T, scale="none",
          trace="none", dendrogram="both", cexRow=0.8, cexCol=1, density.info="none",
          ylab="gene", labRow = F,margin=c(13,2), main=" \n1358 sig genes")


#calculate the DEgenes between all samples.
#output data the same DE gene list for sample with all ages.
gene.names<-row.names(logCPM)
x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_counts.csv",header=T,row.names=1)
cellnames<-read.csv("colnames_6datasets.csv",header=T)
colnames(x3)<-colnames(cellnames)
x3<-x3[row.names(x3) %in% gene.names,]
targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(logCPM)==targets$SampleName)
write.csv(logCPM,file="AM_logCPMscale_integrated_aging3000_quest_heatmap.csv")
sum(colnames(logCPM)==targets$SampleName)
write.csv(targets,file="AM_targets_quest_heatmap.csv")

####################################################
#use top genes from DE analysis to plot for all samples
####################################################
x4<-read.csv("AM_6dataset_integrated_pesudo_bulk_counts.csv",header=T,row.names=1)


colnames(x4)
cellnames<-read.csv("colnames_6datasets.csv",header=T)
colnames(x4)<-colnames(cellnames)

targets<-read.csv("AM_targets_quest_heatmap.csv",header=T)
targets$Sex[targets$Sex=="Male"]<-"male"


targets<-targets[match(colnames(x4),targets$SampleName),]
sum(colnames(x4)==targets$SampleName)


remove.names<-c(grep("_24h",colnames(x4)),grep("_12h",colnames(x4)),grep("_72h",colnames(x4)))
x4<-x4[,-remove.names]
targets <- targets[-remove.names,-1]
sum(colnames(x4)==targets$SampleName)
targets6<-targets[,c(1,3)]
colnames(targets6)<-c("SampleName","condition")
sum(colnames(x4)==targets6$SampleName)
targets6$condition[targets6$condition=="Age55+"]<-"Age55more"
targets6$condition[targets6$condition=="Age<55"]<-"Age55less"

#group <- factor(paste0(targets4$condition)) #define groups based on targets
group <- factor(paste0(targets6$condition)) #define groups based on targets
x4[x4<0]<-0
samplenames <- colnames(x4)
y <- DGEList(counts=x4,group=group) #Define DGEList and design matrix
y <- calcNormFactors(y) #calculate normalization factors

#Filtering
#keep <- rowSums(cpm(y)>1) >= 6
#y <- y[keep, , keep.lib.sizes=FALSE]
#Recompute library size
y$samples$lib.size <- colSums(y$counts)
#Vizualize library size
barplot(y$samples$lib.size*1e-6, names=1:27, ylab="Library size (millions)", xlab="sample")
abline(h=4, col="red", lty=2, lwd=2)
#TMM normalization 
y <- calcNormFactors(y)
y$samples
#Vizualize effect of TMM normalization (inspect all samples)
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
#EstimatingDispersion, GLM
#To estimate common dispersion and tagwise dispersions in one run (recommended):
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmFit(y, design, robust=TRUE) #Fit a negative binomial generalized log-linear model to counts

#AMD vs. AMR
SvsM <- makeContrasts( Age55more - Age55less, levels=design)
lrt <- glmLRT(fit, contrast=SvsM)

is.de <- decideTestsDGE(lrt)
summary(is.de) #check how many sig genes

dataH<-as.data.frame(topTags(lrt,n=14000))
write.csv(dataH,file="Pesudo_bulk_AM_DEgene55vs55less_list_filtered.csv")



#HC heatmap
logCPM <- cpm(y, prior.count = 2, log=TRUE)
#write.csv(logCPM,file="AM_logCPM_integrated_aging3000.csv")
rownames(logCPM) <- rownames(y$counts)
colnames(logCPM) <- targets6$SampleName
#colnames(logCPM) <- targets3$SampleName

o <- order(lrt$table$PValue)
#logCPM <- logCPM[o[1:843],]


logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")

library(pheatmap)
sum(colnames(logCPM)==targets$SampleName)
#targets5<-targets[!is.na(targets$AgeGroup2),]
#sum(targets5$SampleName==colnames(logCPM))

targets7<-targets[,c(3:6)]
targets7$AgeGroup2[is.na(targets$AgeGroup2)]<-"Age30to60"
row.names(targets7)<-targets$SampleName

pheatmap(logCPM,annotation_col = targets7,show_rownames=F,main="")

logCPM2<-logCPM[row.names(logCPM) %in% gene.names,]

p1<-pheatmap(logCPM2,annotation_col = targets7,show_rownames=F,main="DEgeen (784 genes) between 60+vs30-\n All AM Human(38samples)")

logCPM3<-logCPM2[,colnames(logCPM2) %in% targets4$SampleName]
sum(colnames(logCPM3)==row.names(targets5))
pheatmap(logCPM3,annotation_col = targets5,show_rownames=F,main="Core DEgeen (343 genes) AM Human Extreme (19samples)\n old vs young")


ziyou<-cutree(p1$tree_row,k=2)
plot(p1$tree_row)
logCPM2_up<-logCPM2[ziyou==1,]
logCPM2_down<-logCPM2[ziyou==2,]
up_mean<-colMeans(logCPM2_up)
down_mean<-colMeans(logCPM2_down)
sum(row.names(targets7)==targets$SampleName)

plot(as.numeric(as.character(targets$Age)),up_mean,xlab="Age",ylab="Avg Down Gene Cluster Z-score",pch=16,cex=3)
abline(lm(up_mean~targets$Age),lwd=3,col="red")
summary(lm(up_mean~targets$Age))
title("Exp=-0.012*age+0.57 p=0.005")

summary(lm(up_mean~targets$Age))$r.squared 


plot(as.numeric(as.character(targets$Age)),down_mean,xlab="Age",ylab="Avg Up Gene Cluster Z-score",pch=16,cex=3)
abline(lm(down_mean~targets$Age),lwd=3,col="red")
summary(lm(down_mean~targets$Age))
title("Exp=0.009*age-0.41 p=0.03")

write.csv(logCPM2_down,"Up_genelist_with_aging.csv")
write.csv(logCPM2_up,"Down_genelist_with_aging.csv")

ziyou2<-p1$tree_col
plot(ziyou2)
row.names(targets7_r)
targets7_r<-targets[ziyou2$order,]
write.csv(targets7_r,file="TMP_age.csv")
barplot(targets7_r$Age,names.arg = targets7_r$SampleName,las=2,ylab="Age",main="Sample Age")


#################################################
#analysis the result for each individual dataset and check for overlapping signatures 
#################################################
x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_RNA_counts.csv",header=T,row.names=1)
#x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_counts.csv",header=T,row.names=1)

#x3<-round(x3)
#sum(x3<0)
x3[x3<0]<-0
colnames(x3)
cellnames<-read.csv("colnames_6datasets.csv",header=T)
colnames(x3)<-colnames(cellnames)

targets<-read.csv("AM_targets_quest_heatmap.csv",header=T)
targets$Sex[targets$Sex=="Male"]<-"male"


targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(x3)==targets$SampleName)


remove.names<-c(grep("_24h",colnames(x3)),grep("_12h",colnames(x3)),grep("_72h",colnames(x3)))
x3<-x3[,-remove.names]
targets <- targets[-remove.names,-1]
sum(colnames(x3)==targets$SampleName)

Paul_targets<-targets[targets$Study %in% c("Raredon"),]
x_Paul<-x3[,colnames(x3) %in% Paul_targets$SampleName]
sum(x_Paul)
x_Paul<-x_Paul*1000
sum(x_Paul)

mean(targets$Age) #55
Paul_targets$AgeGroup3[Paul_targets$Age<55]<-"Age55less"
Paul_targets$AgeGroup3[Paul_targets$Age>=55]<-"Age55more"
targets6<-Paul_targets[,c(1,7)]
colnames(targets6)<-c("SampleName","condition")
library(edgeR)
group <- factor(paste0(targets6$condition)) #define groups based on targets

samplenames <- colnames(x_Paul)
y <- DGEList(counts=x_Paul,group=group) #Define DGEList and design matrix
y <- calcNormFactors(y) #calculate normalization factors

#Filtering
keep <- rowSums(cpm(y)>1) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
#Recompute library size
y$samples$lib.size <- colSums(y$counts)
#Vizualize library size
barplot(y$samples$lib.size*1e-6, names=1:6, ylab="Library size (millions)", xlab="sample")
abline(h=4, col="red", lty=2, lwd=2)
#TMM normalization 
y <- calcNormFactors(y)
y$samples
#Vizualize effect of TMM normalization (inspect all samples)
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
#EstimatingDispersion, GLM
#To estimate common dispersion and tagwise dispersions in one run (recommended):
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmFit(y, design, robust=TRUE) #Fit a negative binomial generalized log-linear model to counts

#AMD vs. AMR
SvsM <- makeContrasts(Age55more - Age55less, levels=design)
lrt <- glmLRT(fit, contrast=SvsM)

is.de <- decideTestsDGE(lrt)
summary(is.de) #check how many sig genes


dataH<-as.data.frame(topTags(lrt,n=20000))
write.csv(dataH,file="Pesudo_bulk_AM_Raredon_DEgene55cutoff_list_filtered.csv")

#test the overlap
list1<-read.csv("Pesudo_bulk_AM_Reyfman_DEgene55cutoff_list_filtered.csv",row.names=1)
list2<-read.csv("Pesudo_bulk_AM_Raredon_DEgene55cutoff_list_filtered.csv",row.names=1)
list3<-read.csv("Pesudo_bulk_AM_Morse_DEgene55cutoff_list_filtered.csv",row.names=1)
list1<-list1[list1$FDR<0.1,]
list2<-list2[list2$FDR<0.1,]
list3<-list3[list3$FDR<0.1,]

library(VennDiagram)

df1<-list("Reyfman"=as.matrix(row.names(list1)),"Raredon"=as.matrix(row.names(list2)),"Morse&Valenzi"=as.matrix(row.names(list3)))

venn.plot.df1 <- venn.diagram(x = df1, filename = NULL,
                              cat.cex=1,
                              cex=3,
                              cat.col = "black",
                              fill=c( "red","green3","blue"),main="Top DE genes 55+ vs 55- between three datasets FDR<0.1 cutoff")
overlap.df1 <- calculate.overlap(df1)
grid.draw(venn.plot.df1)


#################################################
#analysis using clean dataset with removing cluster 0 and 3 
#################################################
library(edgeR)
setwd("/Users/zrt781/Desktop/doublehit_paper/scRNA_seq_analysis")


x3<-read.csv("AM_6dataset_integrated_pesudo_anchor_counts_clean03.csv",header=T,row.names=1)


colnames(x3)
cellnames<-read.csv("colnames_6datasets.csv",header=T)
colnames(x3)<-colnames(cellnames)

targets<-read.csv("AM_targets_quest_heatmap.csv",header=T)
targets$Sex[targets$Sex=="Male"]<-"male"


targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(x3)==targets$SampleName)


remove.names<-c(grep("_24h",colnames(x3)),grep("_12h",colnames(x3)),grep("_72h",colnames(x3)))
x3<-x3[,-remove.names]
targets <- targets[-remove.names,-1]
sum(colnames(x3)==targets$SampleName)


row.names(targets2)<-NULL
targets3<-targets3[match(colnames(x3),targets3$SampleName),]
sum(colnames(x3)==targets3$SampleName)
row.names(targets3)<-NULL
#remove negative value;
x3[x3<0]<-0
#x3<-round(x3)
sum(x3<0)
#x3<-round(log2(x3+1))

library(ggplot2)
sum(x3)

#run DE analysis on the integrated 
#targets$AgeGroup2[targets$Age<=30]<-"Age30less"
#targets$AgeGroup2[targets$Age>60]<-"Age60more"
#targets4<-targets[!is.na(targets$AgeGroup2),]
#x3<-x3[,colnames(x3) %in% targets4$SampleName]
#targets4<-targets4[,c(1,6)]
#colnames(targets4)<-c("SampleName","condition")
#targets4<-targets4[match(colnames(x3),targets4$SampleName),]
#sum(colnames(x3)==targets4$SampleName)

targets4<-targets[!is.na(targets$AgeGroup2),]
targets4<-targets4[,c(1,6)]
colnames(targets4)<-c("SampleName","condition")
x3<-x3[,colnames(x3) %in% targets4$SampleName]
sum(colnames(x3)==targets4$SampleName)

#group <- factor(paste0(targets4$condition)) #define groups based on targets
group <- factor(paste0(targets4$condition)) #define groups based on targets

samplenames <- colnames(x3)
y <- DGEList(counts=x3,group=group) #Define DGEList and design matrix
y <- calcNormFactors(y) #calculate normalization factors

#Filtering
#keep <- rowSums(cpm(y)>1) >= 6
#y <- y[keep, , keep.lib.sizes=FALSE]
#Recompute library size
y$samples$lib.size <- colSums(y$counts)
#Vizualize library size
barplot(y$samples$lib.size*1e-6, names=1:19, ylab="Library size (millions)", xlab="sample")
abline(h=4, col="red", lty=2, lwd=2)
#TMM normalization 
y <- calcNormFactors(y)
y$samples
#Vizualize effect of TMM normalization (inspect all samples)
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
#EstimatingDispersion, GLM
#To estimate common dispersion and tagwise dispersions in one run (recommended):
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmFit(y, design, robust=TRUE) #Fit a negative binomial generalized log-linear model to counts

#AMD vs. AMR
SvsM <- makeContrasts(Age60more - Age30less, levels=design)
lrt <- glmLRT(fit, contrast=SvsM)

is.de <- decideTestsDGE(lrt)
summary(is.de) #check how many sig genes

#-1*Age30less 1*Age60more
#Down                        783
#NotSig                     2002
#Up                          215

dataH<-as.data.frame(topTags(lrt,n=14000))
write.csv(dataH,file="Pesudo_bulk_AM_DEgene60vs30_list_filtered_clean_021920.csv")



#HC heatmap
logCPM <- cpm(y, prior.count = 2, log=TRUE)
#write.csv(logCPM,file="AM_logCPM_integrated_aging3000.csv")
rownames(logCPM) <- rownames(y$counts)
colnames(logCPM) <- targets4$SampleName
#colnames(logCPM) <- targets3$SampleName

o <- order(lrt$table$PValue)
logCPM <- logCPM[o[1:998],]


logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")

library(pheatmap)
sum(colnames(logCPM)==targets4$SampleName)
targets5<-targets[!is.na(targets$AgeGroup2),]
sum(targets5$SampleName==colnames(logCPM))

targets5<-targets5[,c(4:6)]
row.names(targets5)<-targets4$SampleName
colnames(targets5)<-c("Gender","Stduy","Age Group")
pheatmap(logCPM,annotation_col = targets5,show_rownames=F,main="All DEgeen (998 genes) between 60+vs30- AM Human")

heatmap.2(logCPM, col=col.pan, Rowv=T, scale="none",
          trace="none", dendrogram="both", cexRow=0.8, cexCol=1, density.info="none",
          ylab="gene", labRow = F,margin=c(13,2), main=" \n1358 sig genes")


#calculate the DEgenes between all samples.
#output data the same DE gene list for sample with all ages.
gene.names<-row.names(logCPM)
x3<-read.csv("AM_6dataset_integrated_pesudo_bulk_counts.csv",header=T,row.names=1)
cellnames<-read.csv("colnames_6datasets.csv",header=T)
colnames(x3)<-colnames(cellnames)
x3<-x3[row.names(x3) %in% gene.names,]
targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(logCPM)==targets$SampleName)
write.csv(logCPM,file="AM_logCPMscale_integrated_aging3000_quest_heatmap.csv")
sum(colnames(logCPM)==targets$SampleName)
write.csv(targets,file="AM_targets_quest_heatmap.csv")

#################################
#plot 998 genes list for all the samples
#################################
gene_list<-read.csv("Pesudo_bulk_AM_DEgene60vs30_list_filtered_clean_021920.csv",header=T,row.names=1)
gene_list<-gene_list[gene_list$FDR<0.05,]


x3<-read.csv("AM_6dataset_integrated_pesudo_anchor_counts_clean03.csv",header=T,row.names=1)


colnames(x3)
cellnames<-read.csv("colnames_6datasets.csv",header=T)
colnames(x3)<-colnames(cellnames)

targets<-read.csv("AM_targets_quest_heatmap.csv",header=T)
targets$Sex[targets$Sex=="Male"]<-"male"


targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(x3)==targets$SampleName)


remove.names<-c(grep("_24h",colnames(x3)),grep("_12h",colnames(x3)),grep("_72h",colnames(x3)))
x3<-x3[,-remove.names]
targets <- targets[-remove.names,-1]
sum(colnames(x3)==targets$SampleName)


#remove negative value;
x3[x3<0]<-0
sum(x3<0)

library(ggplot2)
sum(x3)

targets4<-targets[,c(1,3)]
colnames(targets4)<-c("SampleName","condition")
x3<-x3[,colnames(x3) %in% targets4$SampleName]
sum(colnames(x3)==targets4$SampleName)

group <- factor(paste0(targets4$condition)) #define groups based on targets
library(edgeR)
samplenames <- colnames(x3)
y <- DGEList(counts=x3,group=group) #Define DGEList and design matrix
y <- calcNormFactors(y) #calculate normalization factors

y$samples$lib.size <- colSums(y$counts)
barplot(y$samples$lib.size*1e-6, names=1:38, ylab="Library size (millions)", xlab="sample")
abline(h=4, col="red", lty=2, lwd=2)
#TMM normalization 
y <- calcNormFactors(y)
y$samples
#Vizualize effect of TMM normalization (inspect all samples)
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#
logCPM <- cpm(y, prior.count = 2, log=TRUE)
rownames(logCPM) <- rownames(y$counts)
colnames(logCPM) <- targets4$SampleName
#colnames(logCPM) <- targets3$SampleName

logCPM <- logCPM[row.names(logCPM) %in% row.names(gene_list),]
dim(logCPM)

logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")

library(pheatmap)
sum(colnames(logCPM)==targets4$SampleName)
sum(targets$SampleName==colnames(logCPM))

targets5<-targets[,c(3:6)]
row.names(targets5)<-targets4$SampleName
colnames(targets5)<-c("Age Group1","Gender","Study","Age Group2")
pheatmap(logCPM,annotation_col = targets5,show_rownames=F,main="All DEgeen (998 genes) between all AM Human")


targets5<-targets[,c(2:6)]

write.csv(logCPM,file="Human_AM_logCPM_pseudo_bulk_heatmap_plotting_rmcluster0n3.csv")
write.csv(targets5,file="Human_AM_logCPM_pseudo_bulk_heatmap_plotting_rmcluster0n3_meta.csv")



################################
#running for cleaning sample
################################

#######
#model V6 clean twice
#######
setwd("/Users/zrt781/Desktop/doublehit_paper/scRNA_seq_analysis")

x3<-read.csv("Human_AM_clean_integrated_counts_022720.csv",header=T,row.names=1)



targets<-read.csv("AM_targets_quest_heatmap.csv",header=T)
targets$Sex[targets$Sex=="Male"]<-"male"

targets<-targets[targets$SampleName %in% colnames(x3),]
targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(x3)==targets$SampleName)




#remove negative value;
x3[x3<0]<-0
#x3<-round(x3)
sum(x3<0)
#x3<-round(log2(x3+1))

targets4<-targets[!is.na(targets$AgeGroup2),]
targets4<-targets4[,c(2,7)]
colnames(targets4)<-c("SampleName","condition")
x3<-x3[,colnames(x3) %in% targets4$SampleName]
sum(colnames(x3)==targets4$SampleName)

#group <- factor(paste0(targets4$condition)) #define groups based on targets
group <- factor(paste0(targets4$condition)) #define groups based on targets

samplenames <- colnames(x3)
y <- DGEList(counts=x3,group=group) #Define DGEList and design matrix
y <- calcNormFactors(y) #calculate normalization factors

#Filtering
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
#Recompute library size
y$samples$lib.size <- colSums(y$counts)
#Vizualize library size
barplot(y$samples$lib.size*1e-6, names=1:17, ylab="Library size (millions)", xlab="sample")
abline(h=4, col="red", lty=2, lwd=2)
#TMM normalization 
y <- calcNormFactors(y)
y$samples
#Vizualize effect of TMM normalization (inspect all samples)
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
#EstimatingDispersion, GLM
#To estimate common dispersion and tagwise dispersions in one run (recommended):
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmFit(y, design, robust=TRUE) #Fit a negative binomial generalized log-linear model to counts

#AMD vs. AMR
SvsM <- makeContrasts(Age60more - Age30less, levels=design)
lrt <- glmLRT(fit, contrast=SvsM)

is.de <- decideTestsDGE(lrt)
summary(is.de) #check how many sig genes
#-1*Age30less 1*Age60more
#Down                        423
#NotSig                     1120
#Up                           66

dataH<-as.data.frame(topTags(lrt,n=14000))
write.csv(dataH,file="Pesudo_bulk_AM_DEgene60vs30_list_filtered_clean_022720.csv")



#HC heatmap
logCPM <- cpm(y, prior.count = 2, log=TRUE)
#write.csv(logCPM,file="AM_logCPM_integrated_aging3000.csv")
rownames(logCPM) <- rownames(y$counts)
colnames(logCPM) <- targets4$SampleName
#colnames(logCPM) <- targets3$SampleName

o <- order(lrt$table$PValue)
logCPM <- logCPM[o[1:489],]


logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")

library(pheatmap)
sum(colnames(logCPM)==targets4$SampleName)
targets5<-targets[!is.na(targets$AgeGroup2),]
sum(targets5$SampleName==colnames(logCPM))

targets5<-targets5[,c(5:7)]
row.names(targets5)<-targets4$SampleName

pheatmap(logCPM,annotation_col = targets5,show_rownames=F,main="All DEgeen (489 genes) between 60+vs30- AM Human")


#################################
#plot 489 genes list for all the samples
#################################
gene_list<-read.csv("Pesudo_bulk_AM_DEgene60vs30_list_filtered_clean_022720.csv",header=T)
gene_list<-gene_list[as.numeric(as.character(gene_list$FDR))<0.05,]


x3<-read.csv("Human_AM_clean_integrated_counts_022720.csv",header=T,row.names=1)


colnames(x3)
targets<-read.csv("AM_targets_quest_heatmap.csv",header=T)
targets$Sex[targets$Sex=="Male"]<-"male"

targets<-targets[targets$SampleName %in% colnames(x3),]
targets<-targets[match(colnames(x3),targets$SampleName),]
sum(colnames(x3)==targets$SampleName)

#remove negative value;
x3[x3<0]<-0
sum(x3<0)

library(ggplot2)
sum(x3)

targets4<-targets[,c(2,4)]
colnames(targets4)<-c("SampleName","condition")
x3<-x3[,colnames(x3) %in% targets4$SampleName]
sum(colnames(x3)==targets4$SampleName)

group <- factor(paste0(targets4$condition)) #define groups based on targets
library(edgeR)
samplenames <- colnames(x3)
y <- DGEList(counts=x3,group=group) #Define DGEList and design matrix
y <- calcNormFactors(y) #calculate normalization factors

y$samples$lib.size <- colSums(y$counts)
barplot(y$samples$lib.size*1e-6, names=1:38, ylab="Library size (millions)", xlab="sample")
abline(h=4, col="red", lty=2, lwd=2)
#TMM normalization 
y <- calcNormFactors(y)
y$samples
#Vizualize effect of TMM normalization (inspect all samples)
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#
logCPM <- cpm(y, prior.count = 2, log=TRUE)
rownames(logCPM) <- rownames(y$counts)
colnames(logCPM) <- targets4$SampleName
#colnames(logCPM) <- targets3$SampleName

logCPM <- logCPM[row.names(logCPM) %in% gene_list$X,]
dim(logCPM)

logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")

library(pheatmap)
sum(colnames(logCPM)==targets4$SampleName)
sum(targets$SampleName==colnames(logCPM))

targets5<-targets[,c(4:7)]
row.names(targets5)<-targets4$SampleName
colnames(targets5)<-c("Age Group1","Gender","Study","Age Group2")
pheatmap(logCPM,annotation_col = targets5,show_rownames=F,main="All DEgeen (489 genes) between all AM Human")



