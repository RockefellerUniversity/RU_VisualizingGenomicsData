## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(DESeq2)
library(tximport)
library(org.Mm.eg.db)
library(goseq)
library(rtracklayer)
library(vsn)
library(ggplot2)
library(pheatmap)
options(ucscChromosomeNames=FALSE)


## ----eval=F-------------------------------------------------------------------
## setwd("/PathToMyDownload/VisualizingGenomicsData/viz_course/presentations/Slides")
## # e.g. setwd("~/Downloads/VisualizingGenomicsData/viz_course/presentations/Slides")


## ----a,eval=F-----------------------------------------------------------------
## 


## ----eval=FALSE,echo=TRUE-----------------------------------------------------
## library(Rsamtools)
## bamFilesToCount <- c("Sorted_Heart_1.bam","Sorted_Heart_2.bam",
##                      "Sorted_Kidney_1.bam","Sorted_Kidney_2.bam",
##                      "Sorted_Liver_1.bam","Sorted_Liver_2.bam")
## myBams <- BamFileList(bamFilesToCount,yieldSize = 10000)
## 


## ----eval=FALSE,echo=TRUE-----------------------------------------------------
## library(TxDb.Mmusculus.UCSC.mm10.knownGene)
## library(GenomicAlignments)
## geneExons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,by="gene")
## geneCounts <- summarizeOverlaps(geneExons,myBams,
##                                     ignore.strand = TRUE)
## colData(geneCounts)$tissue <- c("Heart", "Heart","Kidney","Kidney","Liver","Liver")


## ----gC1,eval=TRUE,echo=FALSE,cache=TRUE,warning=FALSE,message=FALSE----------
load("data/gC_TissueFull.RData")
geneCounts <- geneCounts_Tissue
colData(geneCounts)$Tissue <- c("Heart", "Heart","Kidney","Kidney","Liver","Liver") 
geneCounts <- geneCounts[rowSums(assay(geneCounts)) > quantile(rowSums(assay(geneCounts)),0.4)]


## ----gD,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1",warning=FALSE,message=FALSE----
dds <- DESeqDataSet(geneCounts,design = ~Tissue)
dds <- DESeq(dds)
heartVsLiver <- results(dds,c("Tissue","Heart","Liver"))
heartVskidney <- results(dds,c("Tissue","Heart","Kidney"))
LiverVskidney <- results(dds,c("Tissue","Liver","Kidney"))


## ----gDA,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gd",warning=FALSE,message=FALSE----
dds2 <- DESeq(dds,test="LRT",reduced = ~1)
AllChanges <- results(dds2)


## ----gnT,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gDA",warning=FALSE,message=FALSE----
normLog2Counts <- normTransform(dds)
normLog2Counts


## ----gnTM,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gnT",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
matrixOfNorm <- assay(normLog2Counts)
boxplot(matrixOfNorm)


## ----gnTMP,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gnTM",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
library(vsn)
vsn::meanSdPlot(matrixOfNorm)


## ----gRL,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gD",warning=FALSE,message=FALSE----
rlogTissue <- rlog(dds)
rlogTissue


## ----gRLM,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gRL",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
rlogMatrix <- assay(rlogTissue)
vsn::meanSdPlot(rlogMatrix)


## ----gPCA,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gRLM",warning=FALSE,message=FALSE,fig.height=3,fig.width=7----
plotPCA(rlogTissue,
        intgroup="Tissue",
        ntop=nrow(rlogTissue))


## ----gPCA3,eval=TRUE,echo=FALSE,cache=TRUE,dependson="gRLM",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
plotPCA(rlogTissue,
        intgroup="Tissue",
        ntop=nrow(rlogTissue))


## ----gPRcomp,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gRLM",warning=FALSE,message=FALSE----
pcRes <- prcomp(t(rlogMatrix))
class(pcRes)
pcRes$x[1:2,]


## ----gPRcosmp,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
plot(pcRes$x,
     col=colData(rlogTissue)$Tissue,
     pch=20,
     cex=2)
legend("top",legend = c("Heart","Kidney","Liver"),
       fill=unique(colData(rlogTissue)$Tissue))


## ----gPRloading,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE----
pcRes$rotation[1:5,1:4]


## ----gPRload2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE----
PC2markers <- sort(pcRes$rotation[,2],decreasing = FALSE)[1:100]
PC2markers[1:10]


## ----gPRcompRot,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
PC2_hVsl <- heartVsLiver$stat[rownames(heartVsLiver) %in% names(PC2markers)]
PC2_hVsk <- heartVskidney$stat[rownames(heartVskidney) %in% names(PC2markers)]
PC2_LVsk <- LiverVskidney$stat[rownames(LiverVskidney) %in% names(PC2markers)]



## ----gPRcompRot2,eval=FALSE,echo=FALSE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE----
## PC1markers <- sort(pcRes$rotation[,1],decreasing = FALSE)[1:100]
## PC1_hVsl <- heartVsLiver$stat[rownames(heartVsLiver) %in% names(PC1markers)]
## PC1_hVsk <- heartVskidney$stat[rownames(heartVskidney) %in% names(PC1markers)]
## PC1_LVsk <- LiverVskidney$stat[rownames(LiverVskidney) %in% names(PC1markers)]
## boxplot(PC1_hVsl,PC1_hVsk,PC1_LVsk,names=c("HeartVsLiver","HeartVsKidney","LiverVsKidney"))


## ----gPRcompRot3,eval=FALSE,echo=FALSE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE----
## PC1markers <- sort(pcRes$rotation[,1],decreasing = TRUE)[1:100]
## PC1_hVsl <- heartVsLiver$stat[rownames(heartVsLiver) %in% names(PC1markers)]
## PC1_hVsk <- heartVskidney$stat[rownames(heartVskidney) %in% names(PC1markers)]
## PC1_LVsk <- LiverVskidney$stat[rownames(LiverVskidney) %in% names(PC1markers)]
## boxplot(PC1_hVsl,PC1_hVsk,PC1_LVsk,names=c("HeartVsLiver","HeartVsKidney","LiverVsKidney"))


## ----gPRcompRotf,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcompRot",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
boxplot(PC2_hVsl,PC2_hVsk,PC2_LVsk,
        names=c("HeartVsLiver","HeartVsKidney","LiverVsKidney"))


## ----gSampleDista,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE----
sampleCor <- cor(rlogMatrix)
sampleCor


## ----gSampleDistb,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gSampleDista",warning=FALSE,message=FALSE----
library(pheatmap)
sampleDists <- as.dist(1-cor(rlogMatrix))
sampleDistMatrix <- as.matrix(sampleDists)



## ----gSampleDistc,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gSampleDistb",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)


## ----gSampleDistca,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gSampleDistb",warning=FALSE,message=FALSE,fig.height=3,fig.width=7----
library(RColorBrewer)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255)
plot(1:255,rep(1,255),
     col=colors,pch=20,cex=20,ann=FALSE,
     yaxt="n")



## ----gSampleDistd,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gSampleDistca",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color = colors)


## ----gSampleDiste,eval=FALSE,echo=TRUE,cache=TRUE,dependson="gSampleDistc",warning=FALSE,message=FALSE,fig.height=3,fig.width=7----
## annoCol <- as.data.frame(colData(dds))
## pheatmap(sampleDistMatrix,
##          clustering_distance_rows=sampleDists,
##          clustering_distance_cols=sampleDists,
##          color = colors,annotation_col = annoCo)


## ----gSampleDistl,eval=TRUE,echo=FALSE,cache=TRUE,dependson="gSampleDistc",warning=FALSE,message=FALSE,fig.height=3,fig.width=7----
annoCol <- as.data.frame(colData(dds))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color = colors,annotation_col = annoCol[,1,drop=FALSE])


## ----he1,eval=TRUE,echo=TRUE,cache=TRUE,dependson=c("gRLM","gDA"),warning=FALSE,message=FALSE----
sigChanges <- rownames(AllChanges)[AllChanges$padj < 0.01 & !is.na(AllChanges$padj)]
sigMat <- rlogMatrix[rownames(rlogMatrix) %in% sigChanges,]
nrow(sigMat)


## ----he2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE,dependson=c("he1"),fig.height=3,fig.width=7----
library(pheatmap)
pheatmap(sigMat,
         scale="row",
         show_rownames = FALSE)


## ----he2a,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPRcomp",warning=FALSE,message=FALSE,dependson=c("he1"),fig.height=6,fig.width=7----
library(pheatmap)
pheatmap(sigMat,
         scale="row",
         show_rownames = FALSE)


## ----km1,eval=FALSE,echo=TRUE,cache=TRUE,dependson="he2",warning=FALSE,message=FALSE----
## library(pheatmap)
## set.seed(153)
## k <-   pheatmap(sigMat,
##            scale="row",kmeans_k = 7)


## ----km2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="km1",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
library(pheatmap)
set.seed(153)
k <-   pheatmap(sigMat,
           scale="row",kmeans_k = 7)


## ----km3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="km2",warning=FALSE,message=FALSE----
names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
clusterDF[1:10,,drop=FALSE]


## ----km4,eval=FALSE,echo=TRUE,cache=TRUE,dependson="km3",warning=FALSE,message=FALSE----
## OrderByCluster <- sigMat[order(clusterDF$Cluster),]
## 
## pheatmap(OrderByCluster,
##            scale="row",annotation_row = clusterDF,
##            show_rownames = FALSE,cluster_rows = FALSE)
## 


## ----km4r,eval=TRUE,echo=TRUE,cache=TRUE,dependson="km3",warning=FALSE,message=FALSE,fig.height=5,fig.width=7----
OrderByCluster <- sigMat[order(clusterDF$Cluster),]

pheatmap(OrderByCluster,
           scale="row",annotation_row = clusterDF,
           show_rownames = FALSE,cluster_rows = FALSE)



## ----kms5,eval=FALSE,echo=TRUE,cache=FALSE,dependson="km4",warning=FALSE,message=FALSE----
## library(NbClust)
## rowScaledMat <- t(scale(t(sigMat)))
## clusterNum <- NbClust(rowScaledMat,distance = "euclidean",
##           min.nc = 2, max.nc = 12,
##           method = "kmeans", index ="silhouette")
## 
## clusterNum$Best.nc


## ----include=FALSE------------------------------------------------------------
load("data/ClusterNum.RData")
clusterNum$Best.nc


## ----kmsh,eval=TRUE,echo=TRUE,cache=FALSE,dependson="km4",warning=FALSE,message=FALSE----
clusterNum$Best.partition[1:10]
orderedCluster <- sort(clusterNum$Best.partition)
sigMat <- sigMat[match(names(orderedCluster),rownames(sigMat)),]


## ----kmsha,eval=TRUE,echo=TRUE,cache=FALSE,dependson="km4",warning=FALSE,message=FALSE----
pheatmap(sigMat,
           scale="row",annotation_row = clusterDF,
           show_rownames = FALSE,cluster_rows = FALSE)


## ----km5,eval=TRUE,echo=TRUE,cache=TRUE,dependson="km4",warning=FALSE,message=FALSE----
heartSpecific <- rownames(clusterDF[clusterDF$Cluster == 1,,drop=FALSE])
heartSpecific[1:10]


## ----km6,eval=TRUE,echo=TRUE,cache=TRUE,dependson="km5",warning=FALSE,message=FALSE----
bckGround <- rownames(AllChanges)[!is.na(AllChanges$padj)]
heartSpecific <- bckGround %in% heartSpecific
names(heartSpecific) <- bckGround
heartSpecific[1:10]


## ----km7,eval=TRUE,echo=TRUE,cache=TRUE,dependson="km6",warning=FALSE,message=FALSE,fig.height=5,fig.width=7----
library(goseq)
toTest <- nullp(heartSpecific,"mm10",id = "knownGene")



## ----km8,eval=TRUE,echo=TRUE,cache=TRUE,dependson="km7",warning=FALSE,message=FALSE----
heartRes <- goseq(toTest,"mm10","knownGene",test.cats = "GO:BP")
heartRes[1:5,]


## ----gPCACM,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPCAC",warning=FALSE,message=FALSE,fig.height=4,fig.width=7----
plotPCA(rlogTissue,
        intgroup="Tissue",
        ntop=nrow(rlogTissue))


## ----gPCACk,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gPCACM",warning=FALSE,message=FALSE----
PC1_rnk <- sort(pcRes$rotation[,1],decreasing = TRUE)
PC1_mat <- sigMat[match(names(PC1_rnk),rownames(sigMat),nomatch = 0),]
PC1_mat[1:3,]


## ----gPCACkl,eval=FALSE,echo=TRUE,cache=TRUE,dependson=c("gPCACk","gSampleDistl"),warning=FALSE,message=FALSE,fig.height=5,fig.width=7----
## pheatmap(PC1_mat,
##          scale="row",
##          cluster_rows=FALSE,
##          show_rownames = FALSE,annotation_col = annoCol
##          )


## ----gPCACk333,eval=TRUE,echo=FALSE,cache=TRUE,dependson=c("gPCACk","gSampleDistl"),warning=FALSE,message=FALSE,fig.height=5,fig.width=7----
pheatmap(PC1_mat,
         scale="row",
         cluster_rows=FALSE,
         show_rownames = FALSE,annotation_col = annoCol[,1,drop=FALSE]
         )

