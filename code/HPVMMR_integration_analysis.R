#HPVMMR integration 
#model V3
#Paul+Valenzi_Morse_Habermann+Madissoon+Raredon

#Raredon counts table 
#Add Raredon sample
setwd("/projects/b1038/Pulmonary/ZiyouRen/Doublehit_scRNA_analysis/Raredon")
load("Raredon_hum2.robj")
Raredon2<-subset(Raredon,ident=1)
DimPlot(Raredon2)
meta_hum2<-Raredon2@meta.data
x_hum2<-Raredon2@assays$RNA@counts
dim(x_hum2)

load("Raredon_222C.robj")
Raredon2<-subset(Raredon,ident=6)
DimPlot(Raredon2)
meta_222C<-Raredon2@meta.data
x_222C<-Raredon2@assays$RNA@counts
dim(x_222C)


load("Raredon_1372C.robj")
Raredon2<-subset(Raredon,ident=7)
DimPlot(Raredon2)
meta_1372C<-Raredon2@meta.data
x_1372C<-Raredon2@assays$RNA@counts
dim(x_1372C)



setwd("/projects/b1038/Pulmonary/ZiyouRen/Doublehit_scRNA_analysis/Raredon/098C")
load("Raredon_098C.robj")
Raredon2<-subset(Raredon,ident=c(3,4))
DimPlot(Raredon2)
meta_098C<-Raredon2@meta.data
x_098C<-Raredon2@assays$RNA@counts
dim(x_098C)


setwd("/projects/b1038/Pulmonary/ZiyouRen/Doublehit_scRNA_analysis/Raredon/133C")
load("Raredon_133C.robj")
Raredon2<-subset(Raredon,ident=6)
DimPlot(Raredon2)
meta_133C<-Raredon2@meta.data
x_133C<-Raredon2@assays$RNA@counts
dim(x_133C)



setwd("/projects/b1038/Pulmonary/ZiyouRen/Doublehit_scRNA_analysis/Raredon/Hum1")
load("Raredon_Hum1.robj")
Raredon2<-subset(Raredon,ident=3)
DimPlot(Raredon2)
meta_hum1<-Raredon2@meta.data
x_hum1<-Raredon2@assays$RNA@counts
dim(x_hum1)


sum(row.names(x_hum1)==row.names(x_098C))
sum(row.names(x_hum1)==row.names(x_133C))
sum(row.names(x_hum1)==row.names(x_1372C))
sum(row.names(x_hum1)==row.names(x_hum2))
sum(row.names(x_hum1)==row.names(x_222C))
Raredon_b<-cbind(Hum1,Hum2,X133C,X1372C,X222C,X098C)
Mega_meta<-rbind(meta_098C,meta_133C,meta_1372C,meta_222C,meta_hum1,meta_hum2)
counts<-cbind(x_098C,x_133C,x_1372C,x_222C,x_hum1,x_hum2)
dim(counts)
sum(colnames(counts)==Mega_meta$CellID)
write.csv(counts,file="Raredon_6samples_raw_count.csv")
write.csv(Mega_meta,file="Raredon_6samples_meta_data.csv")

################################################################################
#running V3 model
################################################################################
setwd("/projects/b1038/Pulmonary/ZiyouRen/Doublehit_scRNA_analysis")
Raredon<-read.csv("Raredon_6samples_raw_count.csv",row.names=1,header=T)
Raredon_meta<-read.csv("Raredon_6samples_meta_data.csv",header=T)

Raredon_meta2<-cbind(as.character(Raredon_meta$X),as.character(Raredon_meta$Study),Raredon_meta$Age, as.character(Raredon_meta$Sex), as.character(Raredon_meta$SampleID))
colnames(Raredon_meta2)<-c("CellID","Study","Age","Sex","SampleID")

AM.integrate2<-readRDS("AM_HPVMM_Seurat_SCtrans.rds")
counts<-AM.integrate@assays$RNA@counts
meta<-AM.integrate2@meta.data
meta<-cbind(as.character(meta$CellID),as.character(meta$Study),as.character(meta$Age), as.character(meta$Sex), as.character(meta$SampleID))
meta<-as.data.frame(meta)
colnames(meta)<-c("CellID","Study","Age","Sex","SampleID")
counts<-as.data.frame(counts)

gene.list<-intersect(row.names(counts),row.names(Raredon))
counts<-counts[row.names(counts) %in% gene.list,]
Raredon<-Raredon[row.names(Raredon) %in% gene.list,]

Raredon_meta2<-as.data.frame(Raredon_meta2)

Raredon_meta2$CellID<-colnames(Raredon)
meta<-as.data.frame(meta)

sum(meta$CellID==colnames(counts))

meta<-rbind(meta,Raredon_meta2)
counts<-cbind(counts,Raredon)
sum(meta$CellID==colnames(counts))
row.names(meta)<-meta$CellID
AM <- CreateSeuratObject(counts = counts, meta.data = meta, project = "HPVMMR_Human_AM")
write.csv(meta,"AM_HPVMMR_meta_data.csv")
write.csv(counts,"AM_HPVMMR_raw_counts.csv")

AM.list <- SplitObject(AM, split.by = "Study")
for (i in names(AM.list)) {
  AM.list[[i]] <- SCTransform(AM.list[[i]], verbose = FALSE)
}

AM.features <- SelectIntegrationFeatures(object.list = AM.list, nfeatures = 3000)
AM.list <- PrepSCTIntegration(object.list = AM.list, anchor.features = AM.features)
AM.anchors <- FindIntegrationAnchors(object.list = AM.list, normalization.method = "SCT", 
                                     anchor.features = AM.features)
AM.integrated <- IntegrateData(anchorset = AM.anchors, normalization.method = "SCT")

AM.integrated <- RunPCA(object = AM.integrated, verbose = FALSE)
AM.integrated <- RunTSNE(object = AM.integrated)
DimPlot(AM.integrated, split.by = "Study", ncol = 3)

Mega_meta<-AM.integrated@meta.data
Mega_meta$Age_group2[as.character(Mega_meta$Age)<55]<-"Age<55"
Mega_meta$Age_group2[as.character(Mega_meta$Age)>=55]<-"55+"
AM.integrated@meta.data$Age_group2<-Mega_meta$Age_group2
AM.integrated@active.ident<-as.factor(AM.integrated@meta.data$Age_group2)
names(AM.integrated@active.ident)<-AM.integrated@meta.data$CellID
DimPlot(AM.integrated, reduction = "tsne",label = TRUE)

sum(Mega_meta$CellID==meta$CellID)
meta<-meta[match(meta$CellID,Mega_meta$CellID),]
Mega_meta$Age<-meta$Age
AM.integrated@meta.data<-Mega_meta

AM.integrated@active.ident<-as.factor(AM.integrated@meta.data$Study)
AM.integrated@active.assay<-"integrated"
AM.integrated@assays$integrated
FeaturePlot(AM.integrated, slot="data",features=c("FABP4","MARCO", "MRC1", "C1QA", "PPARG", "SPP1", "F13A1", "CHI3L1", "CCL2", "CCL3", "CCL4", "CD1C", "JCHAIN", "CCR7", "CCL22"),ncol=4)
AM.integrated@assays$RNA@scale.data
ElbowPlot(AM.integrated, ndims = 100)

AM.integrated <- FindNeighbors(AM.integrated, reduction = "pca", dims = 1:30, nn.eps = 0.5)
AM.integrated <- FindClusters(AM.integrated, resolution = 0.2, n.start = 10)
AM.integrated <- FindClusters(AM.integrated)

DimPlot(AM.integrated, reduction = "tsne",label = TRUE)
DimPlot(AM.integrated, reduction = "tsne",label = TRUE,split.by = "Age_group2",ncol = 2)
DimPlot(AM.integrated, reduction = "tsne",label = TRUE,split.by = "Study",ncol = 3)
markers<-FindAllMarkers(AM.integrated)
AM.integrated@active.ident<-as.factor(AM.integrated@meta.data$integrated_snn_res.0.2)
names(AM.integrated@active.ident)<-AM.integrated@meta.data$CellID
DimPlot(AM.integrated)

saveRDS(AM.integrated,file="AM_6dataset_anchorIntegrate_norm_020520.rds")

library(dplyr)
barplotdata=AM.integrated@meta.data
barplotdata$ident = as.character(barplotdata$integrated_snn_res.0.2)
toplotdata <- barplotdata %>% group_by(ident,Age_group2) %>% count()
toplotdata2 <- barplotdata %>% group_by(Study) %>% count()

sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))


sums <- rep(sums$Sum,each=6)
sums2<-rep(sum(toplotdata2$n),each=6)

toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)

toplotdata2$sums <- sums2
toplotdata2$proportion <- (toplotdata2$n/toplotdata2$sums)

toplotdata2$ident<-c("TotalAvg","TotalAvg")
toplotdata3<-rbind(toplotdata,toplotdata2)
library(ggplot2)
colors=c("#8dd3c7", "#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462")
ggplot(data=toplotdata3,mapping = aes(x=ident,y=proportion,fill=Study)) + geom_col(inherit.aes = T,color="black") +theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)+ ylim(0,1)+theme_classic()

#Age_group
barplotdata=AM.integrated@meta.data
barplotdata$ident = as.character(barplotdata$integrated_snn_res.0.2)
toplotdata <- barplotdata %>% group_by(ident,Age_group2) %>% count()
toplotdata2 <- barplotdata %>% group_by(Age_group2) %>% count()

sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))


sums <- rep(sums$Sum,each=2)
sums2<-rep(sum(toplotdata2$n),each=2)

toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)

toplotdata2$sums <- sums2
toplotdata2$proportion <- (toplotdata2$n/toplotdata2$sums)

toplotdata2$ident<-c("TotalAvg","TotalAvg")
toplotdata3<-rbind(toplotdata,toplotdata2)
library(ggplot2)
colors=c("#FFFFFF","#A3A3A3")
ggplot(data=toplotdata3,mapping = aes(x=ident,y=proportion,fill=Age_group2)) + geom_col(inherit.aes = T,color="black") +theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)+ ylim(0,1)+theme_classic()



#create pesudo bulk RNA sequencing
counts<-AverageExpression(AM.integrated,assays="RNA",slot="counts")
AM.integrated@active.ident<-as.factor(as.character(AM.integrated@meta.data$SampleID))
names(AM.integrated@active.ident)<-AM.integrated@meta.data$CellID
counts<-AverageExpression(AM.integrated,assays="integrated")
counts<-counts$RNA
sum(counts)

write.csv(counts,file="AM_6dataset_integrated_pesudo_bulk_RNA_counts.csv")

##########################################################################
#clean cluster 0 and 3 to re-analyze the data
##########################################################################
setwd("/projects/b1038/Pulmonary/ZiyouRen/Doublehit_scRNA_analysis/")
x<-readRDS(file="AM_6dataset_anchorIntegrate_norm_020520.rds")
x@active.assay<-"SCT"
FeaturePlot(x,features = c("FABP4","MRC1","PPARG","SPP1","CCL2","CCL3"),slot="data",ncol=3,)

FeaturePlot(x,features = c("FABP4","MRC1","PPARG","SPP1","CCL2","CCL3"),slot="data",ncol=3,)


DimPlot(x)
x2<-subset(x,ident=c(1,2,4,5,6,7,8))
AM.integrated<-x2
DimPlot(AM.integrated)
AM.integrated@active.ident<-as.factor(as.character(AM.integrated@meta.data$SampleID))
names(AM.integrated@active.ident)<-AM.integrated@meta.data$CellID
counts<-AverageExpression(AM.integrated,assays="RNA",slot="counts")
#AM.integrated@assays$integrated@scale.data
counts<-counts$integrated
counts<-counts$RNA

#reclustering
AM.integrated <- RunPCA(object = AM.integrated, verbose = FALSE)
AM.integrated <- RunTSNE(object = AM.integrated)
DimPlot(AM.integrated, split.by = "Study", ncol = 3)

AM.integrated <- FindNeighbors(AM.integrated, reduction = "pca", dims = 1:30, nn.eps = 0.5)
AM.integrated <- FindClusters(AM.integrated, resolution = 0.2, n.start = 10)
DimPlot(AM.integrated)

#age group
barplotdata=AM.integrated@meta.data
barplotdata$ident = as.character(barplotdata$integrated_snn_res.0.2)
toplotdata <- barplotdata %>% group_by(ident,Age_group2) %>% count()
toplotdata2 <- barplotdata %>% group_by(Age_group2) %>% count()

sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))


sums <- rep(sums$Sum,each=2)
sums2<-rep(sum(toplotdata2$n),each=2)

toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)

toplotdata2$sums <- sums2
toplotdata2$proportion <- (toplotdata2$n/toplotdata2$sums)

toplotdata2$ident<-c("TotalAvg","TotalAvg")
toplotdata3<-rbind(toplotdata,toplotdata2)
library(ggplot2)
colors=c("#FFFFFF","#A3A3A3")
ggplot(data=toplotdata3,mapping = aes(x=ident,y=proportion,fill=Age_group2)) + geom_col(inherit.aes = T,color="black") +theme(axis.text = element_text(angle=90,hjust = 1))+scale_fill_manual(values=colors)+ ylim(0,1)+theme_classic()
AM.integrated@active.assay<-"SCT"

FeaturePlot(AM.integrated,features=c("RPL32","RPL13"))

Markers<-FindAllMarkers(AM.integrated)
saveRDS(AM.integrated,file="Human_AM_6datasets_anchor_removing_cluster0n3.rds")
x<-readRDS("Human_AM_6datasets_anchor_removing_cluster0n3.rds")
DimPlot(x,split.by = "Study",ncol=3)
sum(counts)
write.csv(counts,file="AM_6dataset_integrated_pesudo_RNA_counts_clean_021720.csv")
names<-row.names(counts)[grep("RLP",row.names(counts))]
names1<-row.names(counts)[grepl("RPL",row.names(counts))]
names2<-row.names(counts)[grepl("RPS",row.names(counts))]
names<-union(names1,names2)
x@active.assay<-"SCT"
FeaturePlot(x, slot="data",features=c("FABP4","PPARG", "C15orf48", "CCL2", "SPP1"),ncol=3)

head(row.names(counts))
names<-as.data.frame(row.names(counts))

x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
x[["percent.RPL"]] <- PercentageFeatureSet(x, pattern = "^RPL")
x[["percent.RPS"]] <- PercentageFeatureSet(x, pattern = "^RPS")

VlnPlot(x, features = c("percent.RPL", "percent.RPS", "percent.mt"), ncol = 3)

###############################################
#further cleaning cluster 5 and 2 and 6
###############################################
x<-readRDS("Human_AM_6datasets_anchor_removing_cluster0n3.rds")
DimPlot(x)
x2<-subset(x,ident=c(0,1,3,4))
DimPlot(x2)
AM.integrated<-x2
DimPlot(AM.integrated)
AM.integrated@active.ident<-as.factor(as.character(AM.integrated@meta.data$SampleID))
names(AM.integrated@active.ident)<-AM.integrated@meta.data$CellID

counts<-AverageExpression(AM.integrated,assays="RNA",slot="counts")

#AM.integrated@assays$integrated@scale.data
counts<-counts$integrated
counts<-counts$RNA

#reclustering
AM.integrated <- RunPCA(object = AM.integrated, verbose = FALSE)
AM.integrated <- RunTSNE(object = AM.integrated)

AM.integrated <- FindNeighbors(AM.integrated, reduction = "pca", dims = 1:30, nn.eps = 0.5)
AM.integrated <- FindClusters(AM.integrated, resolution = 0.2, n.start = 10)
DimPlot(AM.integrated,split.by = "Study",ncol=3)
AM.integrated@active.assay<-"integrated"

AM.integrated@active.assay<-"SCT"
FeaturePlot(AM.integrated, slot="data",features=c("FABP4","PPARG", "C15orf48", "CCL2", "SPP1"),ncol=3)
saveRDS(AM.integrated,file="Human_AM_6datasets_clean_022720.rds")

AM.integrated<-readRDS("Human_AM_6datasets_clean_022720.rds")
barplotdata=AM.integrated@meta.data
barplotdata$ident = as.character(barplotdata$SampleID)
toplotdata <- barplotdata %>% group_by(ident) %>% count()
barplot(toplotdata$n,ylab="number of cells",xlab="Samples")
abline(100,0,col="red",cex=3)
names<-read.csv("colnames_6datasets.csv")
toplotdata$ident<-colnames(names)
abline(100,0,col="red",cex=3)
plotdata<-toplotdata$n
names(plotdata)<-colnames(names)
barplot(plotdata,ylab="number of cells",las=2)
remove.names<-c(grep("_24h",colnames(names)),grep("_12h",colnames(names)),grep("_72h",colnames(names)))
plotdata<-plotdata[-remove.names]
barplot(plotdata,ylab="number of cells",las=2,ylim=c(0,3000))
abline(100,0,col="red",cex=3)
title("Human_AM_6datasets_clean_cluster0n3_2n5n6")
AM.integrated@active.ident<-as.factor(as.character(AM.integrated@meta.data$SampleID))
names(AM.integrated@active.ident)<-AM.integrated@meta.data$CellID
counts<-AverageExpression(AM.integrated)
counts2<-counts$integrated
colnames(names)
colnames(counts2)
names(counts2)<-colnames(names)
remove.names<-c(grep("_24h",colnames(names)),grep("_12h",colnames(names)),grep("_72h",colnames(names)))
counts2<-counts2[,-remove.names]
toplotdata2<-toplotdata[toplotdata$n>10,]
counts3<-counts2[,colnames(counts2) %in% toplotdata2$ident]
write.csv(counts3,file="Human_AM_clean_integrated_counts_022720.csv")

x<-readRDS(file="AM_6dataset_anchorIntegrate_norm_020520.rds")
barplotdata=x@meta.data
barplotdata$ident = as.character(barplotdata$SampleID)
toplotdata <- barplotdata %>% group_by(ident) %>% count()
toplotdata$ident<-colnames(names)
plotdata<-toplotdata$n
names(plotdata)<-colnames(names)
barplot(plotdata,ylab="number of cells",las=2)
remove.names<-c(grep("_24h",colnames(names)),grep("_12h",colnames(names)),grep("_72h",colnames(names)))
plotdata<-plotdata[-remove.names]
barplot(plotdata,ylab="number of cells",las=2,ylim=c(0,3000))
abline(100,0,col="red",cex=3)
title("Human_AM_6datasets_initial")

#
toplotdata3<-toplotdata[toplotdata$n<=10,]
