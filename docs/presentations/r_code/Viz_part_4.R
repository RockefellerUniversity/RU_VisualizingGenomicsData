## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicAlignments)
library(DESeq2)
library(tximport)
library(org.Mm.eg.db)
library(goseq)
library(pheatmap)
library(rtracklayer)
library(MotifDb)
library(Biostrings)
library(seqLogo)
library(BSgenome.Mmusculus.UCSC.mm10)
options(ucscChromosomeNames=FALSE)


## ----setup2, include=FALSE,eval=FALSE-----------------------------------------
## knitr::opts_chunk$set(echo = TRUE)
## AsSlides <- TRUE
## library(TxDb.Mmusculus.UCSC.mm10.knownGene)
## library(GenomicAlignments)
## library(DESeq2)
## library(tximport)
## library(org.Mm.eg.db)
## library(goseq)
## library(DEXSeq)
## library(limma)
## library(rtracklayer)
## library(Gviz)
## options(ucscChromosomeNames=FALSE)
## 
## 
## 
## library(soGGi)
## library(rtracklayer)
## 
## library(soGGi)
## library(rtracklayer)
## library(TxDb.Mmusculus.UCSC.mm10.knownGene)
## dsds <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
## names(dsds) <- NULL
## csc41 <- regionPlot(bamFile = "~/Downloads/ENCFF239XXP.bigWig",testRanges = dsds,style = "percentOfRegion",format = "bigwig")
## plotRegion(csc41)
## 
## myMatt <- matt[!dwdw,]
## plot <- pheatmap::pheatmap(myMatt,scale="none",kmeans_k = 2,cluster_cols = FALSE,useRaster=TRUE)
## library(gplots)
## library(pheatmap)
## 
## pheatmap(myMatt[order(plot$kmeans$cluster),][1:20,],cluster_rows = FALSE,cluster_cols =FALSE,scale="none",
##           ,useRaster=FALSE)
## 
## 
## heatmap.2(myMatt[order(plot$kmeans$cluster),], Rowv=FALSE, Colv=FALSE, useRaster=TRUE, labRow=NA, symbreaks=TRUE, breaks=51,
##            col="redblue", trace="none", key=FALSE)
## 


## ----eval=F-------------------------------------------------------------------
## setwd("/PathToMyDownload/VisualizingGenomicsData/viz_course/presentations/Slides")
## # e.g. setwd("~/Downloads/VisualizingGenomicsData/viz_course/presentations/Slides")


## ----a,eval=F-----------------------------------------------------------------
## 


## ----aa1,eval=TRUE,echo=TRUE,cache=TRUE---------------------------------------
peakFiles <- dir("data/CTCFpeaks/",pattern="*.peaks",
                 full.names = TRUE)
macsPeaks_GR <- list()
for(i in 1:length(peakFiles)){
  macsPeaks_DF <- read.delim(peakFiles[i],
                                  comment.char="#")
     macsPeaks_GR[[i]] <- GRanges(
     seqnames=macsPeaks_DF[,"chr"],
     IRanges(macsPeaks_DF[,"start"],
             macsPeaks_DF[,"end"]
     )
  )
}
macsPeaks_GRL <- GRangesList(macsPeaks_GR)
names(macsPeaks_GRL) <- c("Ch12_1","Ch12_2","Mel_1","Mel_2")


## ----aa2,eval=TRUE,echo=TRUE,cache=TRUE,dependson="aa1"-----------------------
allPeaksSet_nR <- reduce(unlist(macsPeaks_GRL))
overlap <- list()
for(i in 1:length(macsPeaks_GRL)){
  overlap[[i]] <- allPeaksSet_nR %over% macsPeaks_GRL[[i]]
}
overlapMatrix <- do.call(cbind,overlap)
colnames(overlapMatrix) <- names(macsPeaks_GRL)
mcols(allPeaksSet_nR) <- overlapMatrix
allPeaksSet_nR[1,]


## ----aa3,eval=TRUE,echo=TRUE,cache=TRUE,dependson="aa2"-----------------------
HC_Peaks <- allPeaksSet_nR[
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("Ch12_1","Ch12_2")])) >= 2 |
  rowSums(as.data.frame(
    mcols(allPeaksSet_nR)[,c("Mel_1","Mel_2")])) >= 2  
  ]
HC_Peaks[1:2,]


## ----eval=F, echo=T, warning=FALSE--------------------------------------------
## library(Rsamtools)
## library(GenomicAlignments)
## 
## bams <- c("~/Projects/Results/chipseq/testRun/BAMs/Sorted_CTCF_Ch12_1.bam",
##           "~/Projects/Results/chipseq/testRun/BAMs/Sorted_CTCF_Ch12_2.bam",
##           "~/Projects/Results/chipseq/testRun/BAMs/Sorted_CTCF_Mel_1.bam",
##           "~/Projects/Results/chipseq/testRun/BAMs/Sorted_CTCF_Mel_2.bam")
## bamFL <- BamFileList(bams,yieldSize = 5000000)
## myCTCFCounts <- summarizeOverlaps(HC_Peaks,
##                               reads = bamFL,
##                               ignore.strand = TRUE)
## colData(myCTCFCounts)$Group <- c("Ch12","Ch12","Mel","Mel")


## ----eval=F, echo=T, warning=FALSE--------------------------------------------
## library(DESeq2)
## deseqCTCF <- DESeqDataSet(myCTCFCounts,design = ~ Group)
## deseqCTCF <- DESeq(deseqCTCF)
## CTCF_Ch12MinusMel <- results(deseqCTCF,
##                         contrast = c("CellLine","Mel","Ch12"),
##                         format="GRanges")
## 


## ----eval=F, echo=T, warning=FALSE--------------------------------------------
## CTCFrlog <- rlog(deseqCTCF)
## CTCFrlog


## ----load,eval=T,echo=F,warning=FALSE,include=FALSE,cache=TRUE,message=FALSE----
load("data/rloggedData.RData")
CTCFrlog <- toPlot
myRes <- read.delim("data/Antibody_CTCF___Group_CTCF_Ch12_minus_CTCF_MelDEG.xls",sep="\t")
dvev <- matrix(unlist(strsplit(as.vector(myRes[,1]),"_")),ncol=4,byrow=T)
CTCF_Ch12MinusMel <- GRanges(seqnames=dvev[,2],IRanges(as.numeric(dvev[,3]),
                                                       as.numeric(dvev[,4])))
mcols(CTCF_Ch12MinusMel) <- as.data.frame(myRes[,-c(1,2,9:16)])
names(CTCF_Ch12MinusMel) <- myRes[,1]
CTCFrlog


## ----p1,eval=T,echo=T, warning=FALSE,dependson="load",cache=TRUE--------------
sampleDists <- as.dist(1-cor(assay(CTCFrlog)))
sampleDistMatrix <- as.matrix(sampleDists)


## ----p1b,eval=T,echo=T, warning=FALSE,dependson="load",cache=TRUE,fig.height=4,fig.width=7,dependson="p1"----
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)


## ----p2,eval=T,echo=T, warning=FALSE,dependson="load",cache=TRUE--------------
plotPCA(CTCFrlog,intgroup="Group")
 


## ----p3,eval=T,echo=T, warning=FALSE,dependson="load",cache=TRUE--------------
plotPCA(CTCFrlog,intgroup="Group",
                  ntop=nrow(CTCFrlog))


## ----p3a,eval=T,echo=T, warning=FALSE,dependson="load,cache=TRUE"-------------
myData <- plotPCA(CTCFrlog,intgroup="Group",
                  ntop=nrow(CTCFrlog),
                  returnData=TRUE)
myData


## ----p3b,eval=T,echo=T, warning=FALSE,dependson="p3a",cache=TRUE,fig.height=3,fig.width=7----
myData$FragmentLength <- c(129,136,133,125)
library(ggplot2)
ggplot(myData,aes(x=PC1,y=PC2,
                  colour=Group,size=FragmentLength))+
  geom_point()+
  scale_size_continuous(range=c(2,10))


## ----p4,eval=T,echo=T,warning=FALSE,dependson="load",cache=TRUE---------------
pcRes <- prcomp(t(assay(CTCFrlog)))
RankedPC1 <- rownames(pcRes$rotation)[order(pcRes$rotation[,1],
                                            decreasing=T)]
RankedPC1[1:3]


## ----p5,eval=T,echo=T,warning=FALSE,dependson=c("p4","load"),cache=TRUE,fig.height=4,fig.width=7----
library(pheatmap)
rlogMat <- assay(CTCFrlog)
Diff <-  names(CTCF_Ch12MinusMel)[CTCF_Ch12MinusMel$padj < 0.05 & 
                                      !is.na(CTCF_Ch12MinusMel$padj) &
                                    abs(CTCF_Ch12MinusMel$log2FoldChange) > 3]

sigMat <- rlogMat[rownames(rlogMat) %in% Diff,]
pheatmap(sigMat,scale="row",show_rownames = FALSE)


## ----p6,eval=T,echo=T,warning=FALSE,dependson="load",cache=TRUE---------------
#library(devtools)
#install_github("ThomasCarroll/soGGi")
library(soGGi)


## ----p6a,eval=T,echo=T,warning=FALSE,dependson="load",cache=TRUE--------------
UpInMel <-  names(CTCF_Ch12MinusMel)[CTCF_Ch12MinusMel$padj < 0.05 & 
                                      !is.na(CTCF_Ch12MinusMel$padj) &
                                    CTCF_Ch12MinusMel$log2FoldChange < -3]

DownInMel <-  names(CTCF_Ch12MinusMel)[CTCF_Ch12MinusMel$padj < 0.05 & 
                                      !is.na(CTCF_Ch12MinusMel$padj) &
                                    abs(CTCF_Ch12MinusMel$log2FoldChange) > 3]

SigVec  <- ifelse(names(CTCF_Ch12MinusMel) %in% UpInMel,"UpInMel",
             ifelse(names(CTCF_Ch12MinusMel) %in% DownInMel,"DownInMel",
                    "Other"))
CTCF_Ch12MinusMel$Sig <- SigVec
CTCF_Ch12MinusMel$name <- names(CTCF_Ch12MinusMel)


## ----p6b,eval=T,echo=T,warning=FALSE,dependson=c("load","p6a"),cache=TRUE-----
grPs <- quantile(CTCF_Ch12MinusMel$baseMean,c(0.33,0.66))
high <-  names(CTCF_Ch12MinusMel)[CTCF_Ch12MinusMel$baseMean > grPs[2]]
low <-  names(CTCF_Ch12MinusMel)[CTCF_Ch12MinusMel$baseMean < grPs[2]]

Strength  <- ifelse(names(CTCF_Ch12MinusMel) %in% high,"high",
             ifelse(names(CTCF_Ch12MinusMel) %in% low,"low",
                    "")
             )
CTCF_Ch12MinusMel$Strength <- Strength


## ----p6c,eval=T,echo=T,warning=FALSE,dependson=c("load","p6a","p6b"),cache=TRUE----

CTCF_Ch12MinusMel[1:2,]


## ----p7,eval=F,echo=T,warning=FALSE-------------------------------------------
## rownames(CTCF_Ch12MinusMel) <- NULL
## 
## CTCFbigWig <- "data/Sorted_CTCF_Ch12_1Normalised.bw"
## CTCF_ch12 <- regionPlot(CTCFbigWig,
##                         testRanges = CTCF_Ch12MinusMel,
##                         style = "point",
##                         format = "bigwig",
##                         distanceAround = 750)


## ----load2,eval=T,echo=F,warning=FALSE----------------------------------------
#save(CTCF_ch12,file="data//CTCF_ch12_1_P_small.RData")
load("data/CTCF_ch12_1_P_small.RData")


## ----p8,eval=F,echo=T,warning=FALSE,cache=TRUE--------------------------------
## p <- plotRegion(CTCF_ch12)
## p


## ----p9,eval=F,echo=T,warning=FALSE,cache=TRUE--------------------------------
## p+theme_bw()+
##   geom_point(colour="red")


## ----p10,eval=F,echo=T,warning=FALSE,cache=TRUE,dependson="load2"-------------
## plotRegion(CTCF_ch12,
##            summariseBy = "Strength")


## ----p11,eval=F,echo=T,warning=FALSE,cache=TRUE-------------------------------
## rownames(CTCF_ch12) <- mcols(CTCF_ch12)$name
## plotRegion(CTCF_ch12,
##            gts=list(high=high,low=low))


## ----p12,eval=F,echo=T,warning=FALSE,cache=TRUE-------------------------------
## plotRegion(CTCF_ch12,
##            gts=list(high=high,low=low),
##            summariseBy = "name")


## ----p16,eval=F,echo=T,warning=FALSE,cache=TRUE-------------------------------
## library(TxDb.Mmusculus.UCSC.mm10.knownGene)
## dsds <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
## TSS <- resize(dsds,500,fix = "start")
## TTS <- resize(dsds,500,fix = "end")


## ----p17,eval=F,echo=T,warning=FALSE,cache=TRUE-------------------------------
## plotRegion(CTCF_ch12,
##            gts=list(TSS=TSS,
##                     TTS=TTS))
## 


## ----p13,eval=F,echo=T,warning=FALSE,cache=TRUE-------------------------------
## plotRegion(CTCF_ch12,
##            summariseBy = "Strength",
##            colourBy = "Group")


## ----p14,eval=F,echo=T,warning=FALSE,cache=TRUE-------------------------------
## plotRegion(CTCF_ch12,
##            summariseBy = "Strength",
##            colourBy = "Group",
##            groupBy = "Group")


## ----p15,eval=F,echo=T,warning=FALSE,cache=TRUE-------------------------------
## p <- plotRegion(CTCF_ch12,
##                  summariseBy = "Strength",
##                  colourBy = "Group",
##                  groupBy = "Group")
## cols <- c("high" = "darkred", "low" = "darkblue")
## p + scale_colour_manual(values = cols)+
##   theme_minimal()
## 


## ----p18,eval=F,echo=T,warning=FALSE------------------------------------------
## CTCFbigWig <- "data//Sorted_CTCF_MEL_1Normalised.bw"
## CTCF_mel <- regionPlot(CTCFbigWig,
##                         testRanges = CTCF_Ch12MinusMel,
##                         style = "point",
##                         format = "bigwig",
##                         distanceAround = 750)


## ----load3,eval=T,echo=F,warning=FALSE----------------------------------------
#save(CTCF_mel,file="data//CTCF_mel_1_P_small.RData")
load("data//CTCF_mel_1_P_small.RData")
library(soGGi)


## ----p19,eval=T,echo=T,warning=FALSE,message=FALSE,dependson=c("load2","load3")----
CTCF <- c(CTCF_mel,CTCF_ch12)
CTCF


## ----p20,eval=F,echo=T,warning=FALSE------------------------------------------
## plotRegion(CTCF,
##            colourBy = "Sample")


## ----p21,eval=F,echo=T,warning=FALSE------------------------------------------
## plotRegion(CTCF,
##            groupBy = "Sample")


## ----p22,eval=F,echo=T,warning=FALSE------------------------------------------
## plotRegion(CTCF,
##            summariseBy = "Sig",
##            groupBy = "Group",
##            colourBy="Sample")


## ----p23,eval=T,echo=T,warning=FALSE,cache=T,message=FALSE,warning=FALSE,fig.height=3,fig.width=5----
library(MotifDb)
library(Biostrings)
library(seqLogo)
library(BSgenome.Mmusculus.UCSC.mm10)

CTCFm <- query(MotifDb, c("CTCF"))
ctcfMotif <- CTCFm[[1]]
seqLogo(ctcfMotif)


## ----p24,eval=T,echo=T,warning=FALSE,cache=T,dependson=c("p23")---------------
myRes <- matchPWM(ctcfMotif,
                  BSgenome.Mmusculus.UCSC.mm10[["chr19"]])
CTCFmotifs <- GRanges("chr19",ranges(myRes))
CTCFrle <- coverage(CTCFmotifs)
CTCFrle[[1]]


## ----p25,eval=F,echo=T,warning=FALSE,dependson=c("p2"),message=FALSE,warning=FALSE----
## motifPlot <- regionPlot(CTCFrle,
##                CTCF_Ch12MinusMel,
##                format = "rlelist")
## plotRegion(motifPlot)


## ----p26,eval=F,echo=T,warning=FALSE------------------------------------------
## CTCFSig19 <- CTCF_mel[seqnames(CTCF_mel) %in% "chr19",]
## rownames(motifPlot) <- rownames(CTCFSig19) <- NULL
## CTCFall <- c(motifPlot,CTCFSig19)
## p <- plotRegion(CTCFall)
## p


## ----p27,eval=F,echo=T,warning=FALSE------------------------------------------
## p <- plotRegion(CTCFall,freeScale=TRUE)
## p


## ----p28,eval=T,echo=T,warning=FALSE,cache=TRUE-------------------------------
dsds <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
names(dsds) <- NULL
dsds


## ----p29,eval=F,echo=T,warning=FALSE------------------------------------------
## 
## pol2Ser2 <- regionPlot(bamFile = "~/Downloads/ENCFF239XXP.bigWig",
##                     testRanges = dsds,style = "percentOfRegion",
##                     format = "bigwig")
## 


## ----load4,eval=T,echo=F,warning=FALSE----------------------------------------
#save(pol2Ser2,file="data//pol2Ser2_1_P.RData")
load("data//pol2Ser2_1_P.RData")
#rm(CTCF)
#rm(CTCF_ch12)
#rm(CTCF_mel)


## ----p30,eval=F,echo=T,warning=FALSE------------------------------------------
## plotRegion(pol2Ser2)


## ----p31,eval=F,echo=T,warning=FALSE------------------------------------------
## CTCF_gene <- regionPlot(CTCFbigWig,
##                     testRanges = dsds,
##                     style = "percentOfRegion",
##                     format = "bigwig")


## ----load5,eval=T,echo=F,warning=FALSE----------------------------------------
#save(CTCF_gene,file="data//CTCF_gene_P.RData")
load("data//CTCF_gene_P.RData")
#rm(CTCF)
#rm(CTCF_ch12)
#rm(CTCF_mel)


## ----p32,eval=F,echo=F,warning=FALSE,include=FALSE----------------------------
## rownames(CTCF_gene) <- rownames(pol2Ser2) <- NULL
## 
## plotRegion(c(CTCF_gene,pol2Ser2[-c(1:5),]),
##            colourBy = "Sample",
##            freeScale = TRUE,
##            groupBy="Sample")


## ----p33,eval=F,echo=T,warning=FALSE,include=TRUE-----------------------------
## rownames(CTCF_gene) <- rownames(pol2Ser2) <- NULL
## 
## plotRegion(c(CTCF_gene,pol2Ser2),
##            colourBy = "Sample",
##            freeScale = TRUE,
##            groupBy="Sample")


## ----p34,eval=F,echo=T,warning=FALSE------------------------------------------
## CTCF_Mel_Up <- CTCF[mcols(CTCF)$name %in% UpInMel,]


## ----p35,eval=F,echo=T,warning=FALSE------------------------------------------
## k1 <- plotHeatmap(CTCF_Mel_Up[[1]],
##                   col = colorRampPalette(blues9)(100),
##                   maxValue=15)


## ----p36,eval=F,echo=T,warning=FALSE------------------------------------------
## k2 <- plotHeatmap(CTCF_Mel_Up[[2]],
##                   col = colorRampPalette(blues9)(100),
##                   maxValue=15)

