#####################################################
############# IMPORT DATA ###########################
#####################################################


current_directory=getwd()

data_folder_name="data"
dataset_directory=paste(current_directory,data_folder_name,sep="/")


library("affy")
library("AnnotationDbi")


#Download raw data
source("obtainRawData.R")
obtainRawData("GSE28146",folder=data_folder_name)


#Load phenodata (from "phenodata.txt") and name columns
phdataset<- read.table("phenodata.txt", sep="\t")
colnames(phdataset)<-c("sampleNames","age","diseaseStage","title")


#Filter control and severe indexes 
#(those which start by "sever" or "control" in "diseaseStage" column )
severeIndex<-grep("^severe",as.character(phdataset[,"diseaseStage"]))
ControlIndex<-grep("^control",as.character(phdataset[,"diseaseStage"]))

#Obtain the names of samples that are "severe" or "control"
dnames<-phdataset[c(severeIndex,ControlIndex),"sampleNames"]
dnames<-paste(dnames, ".CEL",sep = "")



#Load AffyBatch (load only .CEL of "severe" and "control" ) 
setwd(dataset_directory)
AffyBatchObject = ReadAffy(filenames=dnames)
setwd(current_directory)


#Attach the phenodata
#(the original phenodata with sample names plus the phenodata of "phenodata.txt")
pData(AffyBatchObject)<-cbind(
			pData(AffyBatchObject),
			phdataset[c(severeIndex,ControlIndex),])

#Check the data has been correctly assigned to each sample
pData(AffyBatchObject)



#####################################################
#############PREPROCESS###############################
#####################################################


#####Check for microarrays potentially problematics####


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

par( mfrow = c( 1, 2 ) )
hist(AffyBatchMAS)
title("Histogram with MAS")

hist(AffyBatchRMA)
title("Histogram with RMA")


#######Convert to ExpressionSet

ExpressionSet<- expresso(AffyBatchObject,
		bgcorrect.method="rma",
		normalize.method="quantiles",
		pmcorrect.method="pmonly",
		summary.method="medianpolish")

###Also it is possible to execute RMA directly: 
### ExpressionSet= rma(AffyBatchObject)








############################################################
###################### FILTER GENES ########################
############################################################

####Filter by variance####
library("genefilter")
ExpressionSet<-nsFilter(ExpressionSet,var.cutoff=0.6)$eset


####Filter by mean tests####
N=50 #Number of genes we want to keep. Keep 50 genes

varLabels(ExpressionSet)
#Two groups will be created depending on diseaseStage

# Perform a mean difference hypothesis test 
# between two groups for each gene
RowTestTable = rowttests(ExpressionSet, "diseaseStage")

# Obtain the genes with the lowest p-value
order=order(RowTestTable$p.value)[1:N]
features=featureNames(ExpressionSet)[order] 


# Display a plot volcano.
# We highlight with a blue dot
# the 50 points with the lowest pvalue
plot(RowTestTable$dm,
	-log10(RowTestTable$p.value),
	xlab="mean differences (in log-ratio)",
	ylab=expression(-log[10]~(p-value)))


points(	RowTestTable$dm[order], 
		-log10(RowTestTable$p.value)[order], 
		pch=18, col="blue")




#Order RowTestTable by p-value
RowTestTable<-RowTestTable[order(RowTestTable$p.value),]


# Filter the ExpressionSet and keep only obtained
# N genes (50 in this case)
ExpressionSet=ExpressionSet[features,]



############################################################
################     CLASSIFICATION    #####################
############################################################

#####Scale ExpressionSet###
library("matrixStats")
M<-exprs(ExpressionSet) #Obtain the matrix of ExpressionSet
esetScaled<-(M- rowMedians(M))/rowIQRs(M) #Scale the ExpressionSet
exprs(ExpressionSet)<-esetScaled




library("MLInterfaces")

fvalidation<-xvalSpec("LOO")
f<-formula(paste("diseaseStage", "~ .")) #f is "diseaseStage ~ ."

##### KNN ####
ClassifierKNN = MLearn(f, 
		data=ExpressionSet,
		.method=knnI(k=1),
		trainInd=fvalidation)


# Obtain confusion matrix
confusionMatrixKNN<-confuMat(ClassifierKNN ) 





### Random forest ####
ClassifierRForest =MLearn(f,
		data=ExpressionSet,
		.method=randomForestI,
		ntree=100,
		trainInd=fvalidation)


# Obtain confusion matrix
confusionMatrixRForest<-confuMat(ClassifierRForest) 


#### Display confusion matrices###
confusionMatrixRForest<-
		confusionMatrixRForest[c("control","severe stage"),
						c("control","severe stage")]
confusionMatrixRForest



## Display both confusion matrices
confusionMatrixKNN
confusionMatrixRForest







############################################################
##################     CLUSTERING    #######################
############################################################



#####Compute distance matrix #####
# between samples, and also between probes
DistanceSamples=dist(t(exprs(ExpressionSet)),method="euclidean")
DistanceProbes=dist(exprs(ExpressionSet),method="euclidean")




######Clustering of samples#############
#We display PCA to have an intuition of 
#the number of clusters of samples
PCAsamples<-prcomp(t(exprs(ExpressionSet)))
plot(PCAsamples$x [,1], PCAsamples$x [,2])


HierarchicalClustSamples=hclust(DistanceSamples,method="single")
plot(HierarchicalClustSamples)



######Clustering of probes##################
#We display PCA to have an intuition of the number of clusters
PCAprobes<-prcomp(exprs(ExpressionSet))
plot(PCAprobes$x [,1], PCAprobes$x [,2])

text(PCAprobes$x[,1],
	PCAprobes$x[,2]-0.2,
	labels=featureNames(ExpressionSet))




##Hierarchical clustering of probes####
HierarchicalClustProbes=hclust(DistanceProbes,method="single")
plot(HierarchicalClustProbes)



####K-means####
#We decide to set the number of clusters as 2
CLUSTER_NUMBER=2
kmeans=kmeans(exprs(ExpressionSet),centers=CLUSTER_NUMBER,nstart=5)

#Save in a list the name of genes of each cluster
ClusterGenes<-list()
for (i in 1:CLUSTER_NUMBER){
	ClusterGenes[[i]]<-names(kmeans$cluster[kmeans$cluster==i])
}


ClusterGenes








############################################################
#############     GENE INFORMATION    ######################
############################################################




annotation(ExpressionSet)

library("hgu133plus2.db")

#Entrez genes data.frame
genesDB<-merge(
	toTable(hgu133plus2ENTREZID),
	toTable(hgu133plus2GENENAME))






# In 'features' is stored the 50 probes
# with smallest p-value (ordered)

#Obtain 10 features
features10<-features[1:10]

#Find the probes in the table
genesDB[ genesDB$probe_id %in% features10,]



