## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(DESeq2)
library(tximport)
library(org.Mm.eg.db)
library(goseq)
library(DEXSeq)
library(limma)
library(rtracklayer)
"%in%" <- BiocGenerics::"%in%"


## ----eval=F-------------------------------------------------------------------
## setwd("/PathToMyDownload/VisualizingGenomicsData/viz_course/presentations/Slides")
## # e.g. setwd("~/Downloads/VisualizingGenomicsData/viz_course/presentations/Slides")


## ---- echo=T,eval=F-----------------------------------------------------------
## ## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite("Gviz")
## 


## ---- echo=T------------------------------------------------------------------
library(Gviz) 
genomeAxis <- GenomeAxisTrack(name="MyAxis") 
genomeAxis 


## ---- echo=T,eval=F,fig.width=23,fig.height=3---------------------------------
## plotTracks(genomeAxis,from=100,to=10100)


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
plotTracks(genomeAxis,from=100,to=10100,cex=3) 


## ---- echo=T,eval=F,fig.width=23,fig.height=3---------------------------------
## plotTracks(genomeAxis,from=100,to=10100,
##            add53=T,add35=T)


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
plotTracks(genomeAxis,from=100,to=10100, 
           add53=T,add35=T,cex=3) 


## ---- echo=T,eval=F,fig.width=23,fig.height=3---------------------------------
## plotTracks(genomeAxis,from=100,to=10100,
##            littleTicks = TRUE)


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
plotTracks(genomeAxis,from=100,to=10100, 
           littleTicks = TRUE,cex=3) 


## ---- echo=T,eval=F,fig.width=23,fig.height=3---------------------------------
## plotTracks(genomeAxis,from=100,to=10100,
##            labelPos="below")


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
plotTracks(genomeAxis,from=100,to=10100, 
           labelPos="below",cex=3) 


## ---- echo=T,eval=F,fig.width=10,fig.height=5---------------------------------
## plotTracks(genomeAxis,from=100,to=10100,
##            scale=1,labelPos="below")


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
plotTracks(genomeAxis,from=100,to=10100, 
           scale=1,labelPos="below",cex=3) 


## ---- echo=T,eval=F,fig.width=10,fig.height=5---------------------------------
## plotTracks(genomeAxis,from=100,to=10100,
##            scale=0.3)


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
plotTracks(genomeAxis,from=100,to=10100, 
           scale=0.3,cex=3) 


## ---- echo=T,eval=F,fig.width=10,fig.height=5---------------------------------
## plotTracks(genomeAxis,from=100,to=10100,
##            scale=2500)


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
plotTracks(genomeAxis,from=100,to=10100, 
           scale=2500,cex=3) 


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
library(IRanges) 
regionsOfInterest <- IRanges(start=c(140,5140),end=c(2540,7540)) 
names(regionsOfInterest) <- c("ROI_1","ROI_2") 
regionsOfInterest 


## ---- echo=T,eval=F,fig.width=10,fig.height=5---------------------------------
## genomeAxis <- GenomeAxisTrack(name="MyAxis",
##                               range = regionsOfInterest)
## plotTracks(genomeAxis,from=100,to=10100)


## ---- echo=F,fig.width=23,fig.height=3----------------------------------------
genomeAxis <- GenomeAxisTrack(name="MyAxis", 
                              range = regionsOfInterest) 
plotTracks(genomeAxis,from=100,to=10100,cex=3) 


## ---- echo=F,fig.width=23,fig.height=5----------------------------------------
 
plotTracks(genomeAxis,from=100,to=10100, 
           range=regionsOfInterest, 
           showId=T,cex=3,col.id="black") 


## ---- echo=T,eval=F,fig.width=10,fig.height=5---------------------------------
## 
## plotTracks(genomeAxis,from=100,to=10100,
##            range=regionsOfInterest,
##            showId=T)


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
mcols(regionsOfInterest) <- data.frame(Sample1=c(30,20),
                                       Sample2=c(20,200)) 
regionsOfInterest <- GRanges(seqnames="chr5",
                             ranges = regionsOfInterest) 
regionsOfInterest 


## ---- echo=T,eval=F,fig.width=10,fig.height=5---------------------------------
## dataROI <- DataTrack(regionsOfInterest)
## plotTracks(dataROI)


## ---- echo=F,fig.width=23,fig.height=5----------------------------------------
dataROI <- DataTrack(regionsOfInterest) 
plotTracks(dataROI,cex=3) 


## ---- echo=T,fig.width=10,fig.height=5,message=FALSE,warning=FALSE------------
library(rtracklayer) 
allChromosomeCoverage <- import.bw("data/small_Sorted_SRR568129.bw",
                                   as="GRanges") 
class(allChromosomeCoverage)


## ---- echo=T,eval=TRUE,fig.width=10,fig.height=5,message=FALSE,warning=FALSE----
allChromosomeCoverage 


## ---- echo=FALSE,eval=FALSE,fig.width=10,fig.height=5,message=FALSE,warning=FALSE----
## require(BiocGenerics)
## "%in%" <- BiocGenerics::"%in%"
## allChromosomeCoverage[!seqnames(allChromosomeCoverage) %in% "chrM"]


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
accDT <- DataTrack(allChromosomeCoverage,chomosome="chr5") 
accDT 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(accDT, 
           from=134887451,to=134888111, 
           chromosome="chr5") 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(accDT, 
           from=134887451,to=134888111, 
           chromosome="chr5",type="l") 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(accDT, 
           from=134887451,to=134888111, 
           chromosome="chr5",type="smooth") 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(accDT, 
           from=134887451,to=134888111, 
           chromosome="chr5",type="h") 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(accDT, 
           from=134887451,to=134888111, 
           chromosome="chr5",type="mountain") 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(accDT, 
           from=134887451,to=134888111, 
           chromosome="chr5",type="heatmap") 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(accDT, 
           from=134887451,to=134888111, 
           chromosome="chr5", 
           col="red",cex=4) 


## ---- echo=T,fig.width=25,fig.height=5----------------------------------------
plotTracks(c(accDT,genomeAxis), 
           from=134887451,to=134888111, 
           chromosome="chr5" 
           ) 


## ---- echo=T,fig.width=25,fig.height=5----------------------------------------
plotTracks(c(genomeAxis,accDT), 
           from=134887451,to=134888111, 
           chromosome="chr5" 
           ) 


## ---- echo=T,fig.width=25,fig.height=5----------------------------------------
plotTracks(c(genomeAxis,accDT), 
           from=134887451,to=134888111, 
           chromosome="chr5", 
           sizes=c(0.5,1) 
           ) 


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
toGroup <- GRanges(seqnames="chr5", 
        IRanges( 
          start=c(10,500,550,2000,2500), 
          end=c(300,800,850,2300,2800) 
        )) 
names(toGroup) <- seq(1,5) 
toGroup 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
 
annoT <- AnnotationTrack(toGroup, 
                group = c("Ann1","Ann1", 
                          "Ann2", 
                          "Ann3","Ann3")) 
 
plotTracks(annoT) 
 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
 
plotTracks(annoT,groupAnnotation = "group") 
 


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
strand(toGroup) 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
strand(toGroup) <- c("+","+","*","-","-") 
annoT <- AnnotationTrack(toGroup, 
                group = c("Ann1","Ann1", 
                          "Ann2", 
                          "Ann3","Ann3")) 
 
plotTracks(annoT, groupingAnnotation="group") 


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------


## ---- echo=F,fig.width=25,fig.height=5----------------------------------------
toGroup <- GRanges(seqnames="chr5", 
        IRanges( 
          start=c(100,100,500,700,2000,2500), 
          end=c(300,300,800,1050,2300,2800) 
        )) 
names(toGroup) <- seq(1,6) 
 
#toGroup 
 
strand(toGroup) <- c("*","*","*","*","*","*") 
annoT <- AnnotationTrack(toGroup, 
                group = c("Ann1", 
                          "Ann2", 
                          "Ann1", 
                          "Ann2", 
                          "Ann3", 
                          "Ann3")) 


## ---- echo=T,fig.width=25,fig.height=5----------------------------------------
plotTracks(annoT, groupingAnnotation="group",stacking="squish") 


## ---- echo=T,fig.width=25,fig.height=5----------------------------------------
plotTracks(annoT, groupingAnnotation="group",stacking="dense") 


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
feature(annoT) 


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
feature(annoT) <- c(rep("Good",4),rep("Bad",2)) 
feature(annoT) 


## ---- echo=T,fig.width=25,fig.height=5----------------------------------------
plotTracks(annoT, featureAnnotation = "feature", 
           groupAnnotation = "group", 
           Good="Blue",Bad="Red") 


## ---- echo=T,fig.width=10,fig.height=5----------------------------------------
data(geneModels) 
geneModels[1,]


## ---- echo=F,fig.width=10,fig.height=5----------------------------------------
data(geneModels) 
head(geneModels) 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
grtrack <- GeneRegionTrack(geneModels, genome = "hg19", 
                           chromosome = "chr7", 
                           name = "smallRegions") 
plotTracks(grtrack) 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(grtrack) 


## ---- echo=T,fig.width=23,fig.height=5----------------------------------------
plotTracks(grtrack,transcriptAnnotation="gene") 


## ---- echo=T,fig.width=23,fig.height=8----------------------------------------
plotTracks(grtrack,transcriptAnnotation="transcript") 


## ---- echo=T,fig.width=15,fig.height=5----------------------------------------
plotTracks(grtrack,transcriptAnnotation="symbol") 


## ---- echo=T,fig.width=15,fig.height=5----------------------------------------
plotTracks(grtrack,exonAnnotation="exon",
           from=26677490,to=26686889,cex=0.5) 


## ---- echo=T,fig.width=20,fig.height=6----------------------------------------
plotTracks(grtrack, stacking="dense") 


## ---- echo=T,fig.width=12,fig.height=3----------------------------------------
plotTracks(grtrack, collapseTranscripts=T, 
           transcriptAnnotation = "symbol") 


## ---- echo=T,fig.width=12,fig.height=3----------------------------------------
plotTracks(grtrack, collapseTranscripts=T, 
           transcriptAnnotation = "symbol", 
           shape="arrow") 


## ---- echo=T,fig.width=12,fig.height=3----------------------------------------
plotTracks(grtrack, collapseTranscripts="longest", 
           transcriptAnnotation = "symbol") 


## ---- echo=T,fig.width=12,fig.height=3----------------------------------------
plotTracks(grtrack, collapseTranscripts="meta", 
           transcriptAnnotation = "symbol") 


## ---- echo=TRUE---------------------------------------------------------------
 
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
txdb 


## ---- echo=TRUE---------------------------------------------------------------
 
customFromTxDb <- GeneRegionTrack(txdb,chromosome="chr7") 
head(customFromTxDb) 


## ----echo=T,fig.width=12,fig.height=3-----------------------------------------
 
plotTracks(customFromTxDb, 
           from=26591341,to=27034958, 
           transcriptAnnotation="gene") 


## ---- echo=TRUE,fig.width=12,fig.height=3,message=FALSE,warning=FALSE---------
library(GenomicFeatures) 
ensembleGTF <- "data/hg19.gtf.gz"
txdbFromGFF <- makeTxDbFromGFF(file = ensembleGTF) 
customFromTxDb <- GeneRegionTrack(txdbFromGFF,chromosome="chr7") 
plotTracks(customFromTxDb, 
           from=26591341,to=27034958, 
           transcriptAnnotation="gene") 

