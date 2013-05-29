current_directory=getwd()

data1_folder_name="dataset1"
dataset1_directory=paste(current_directory,data1_folder_name,sep="/")


library("affy")
library("AnnotationDbi")


#####################################################
############# LOAD DATA #############################
#####################################################


#Download raw data
source("obtainRawData.R")
obtainRawData("GSE28146",folder=data1_folder_name)


#Load phenodata (from dataset1.txt)
phdataset1<- read.table("dataset1.txt", sep="\t")
colnames(phdataset1)<-c("samples","age","diseaseStage","title")


#Filter control and severe, and store those .CEL names
severeIndex<-grep("^severe",as.character(phdataset1[,3]))
ControlIndex<-grep("^control",as.character(phdataset1[,3]))

d1names<-phdataset1[c(severeIndex,ControlIndex),1]
d1names<- as.vector(d1names)
d1names<-paste(d1names, ".CEL",sep = "")



#Load AffyBatch (load only .CEL of "severe" and "control" ) 
setwd(dataset1_directory)
AffyBatchObject = ReadAffy(filenames=d1names)
setwd(current_directory)


#Attach the phenodata
pData(AffyBatchObject)<-cbind(
			pData(AffyBatchObject),
			phdataset1[c(severeIndex,ControlIndex),])

#Check the data has been correctly assigned to each sample
pData(AffyBatchObject)


#####################################################
#############PREPROCESS###############################
#####################################################


####Check for microarrays potentially problematics


library("affyPLM")
AffyBatchPLM= fitPLM(AffyBatchObject)

#RLE plot
Mbox(AffyBatchPLM, main="RLE")
abline(h=0)

#NUSE plot
boxplot(AffyBatchPLM, main="NUSE")
abline(h=1)



###Check why normalization is needed

hist(AffyBatchObject)

library (simpleaffy)
#Execute qc function
QCstats = qc(AffyBatchObject)
#Display a plot of qc stats
plot(QCstats)



######Execute both preprocess and compare


library("affyPLM")

##MAS preprocess
AffyBatchMAS<-preprocess(AffyBatchObject,
	background.method="MAS",
	normalize.method="scaling"
)


##RMA  preprocess
AffyBatchRMA<-preprocess(AffyBatchObject,
	background.method="RMA.2",
	normalize.method="quantile"
)


hist(AffyBatchMAS)

hist(AffyBatchRMA)



#######Convert to ExpressionSet

ExpressionSet<- expresso(AffyBatchObject,
		bgcorrect.method="rma",
		normalize.method="quantiles",
		pmcorrect.method="pmonly",
		summary.method="medianpolish")


############################################################
######################FILTER GENES#########################
############################################################

library("genefilter")
ExpressionSet<-nsFilter(ExpressionSet,var.cutoff=0.6)$eset


#Filter by mean tests
N=50 #Number of genes we want to keep

varLabels(ExpressionSet)
#Two groups will be regarding diseaseStage


RowTestTable = rowttests(ExpressionSet, "diseaseStage")
order=order(RowTestTable$p.value)[1:N]
features=featureNames(ExpressionSet)[order] 

#In variable features here is stored the N genes which have more variance between
# control groul and alzheimer group


#Filter the ExpressionSet and keep only N genes
ExpressionSet=ExpressionSet[features,]



############################################################
################   MACHINE LEARNING   #####################
############################################################

#Scale ExpressionSet
library("matrixStats")
M<-exprs(ExpressionSet) #Obtain the matrix of ExpressionSet
esetScaled<-(M- rowMedians(M))/rowIQRs(M) #Scale the ExpressionSet
exprs(ExpressionSet)<-esetScaled


#Compute distance matrix (for later uses)
Distances=dist(exprs(ExpressionSet),method="euclidean")




#################CLustering

#We display PCA to have an intuition of the number of clusters
PCA<-prcomp(exprs(ExpressionSet))
plot(PCA$x [,1], PCA$x [,2])



#We decide to set the number of clusters as 2
CLUSTER_NUMBER=2
kmeans=kmeans(exprs(ExpressionSet),centers=CLUSTER_NUMBER,nstart=5)

#Save in a list the name of genes of each cluster
ClusterGenes<-list()
for (i in 1:CLUSTER_NUMBER){
	ClusterGenes[[i]]<-names(kmeans$cluster[kmeans$cluster==i])
}

ClusterGenes



##Hierarchichal clustering
HierarchicalClust=hclust(Distances,method="single")
plot(HierarchicalClust)



