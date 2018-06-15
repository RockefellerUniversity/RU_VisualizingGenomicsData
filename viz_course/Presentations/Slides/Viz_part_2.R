## ----setup, include=FALSE------------------------------------------------
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
library(Gviz)
options(ucscChromosomeNames=FALSE)

## ----eval=F--------------------------------------------------------------
## setwd("/PathToMyDownload/VisualizingGenomicsData/viz_course/presentations/Slides")
## # e.g. setwd("~/Downloads/VisualizingGenomicsData/viz_course/presentations/Slides")

## ----bbunny1, echo=T,eval=T,fig.width=5,fig.height=5,fig.align='center'----
library(Gviz) 
genomeAxis <- GenomeAxisTrack(name="MyAxis") 
plotTracks(genomeAxis,from=100,to=10100) 

## ----bbunny2, echo=T,eval=T,fig.width=10,fig.height=3,dependson="bbunny1",warning=FALSE----
library(rtracklayer) 
allChromosomeCoverage <- import.bw("../../Data/activatedReads.bw",
                                   as="GRanges") 
accDT <- DataTrack(allChromosomeCoverage,chomosome="chr17") 
plotTracks(c(accDT,genomeAxis), 
           from=47504051,to=47600688, 
           chromosome="chr17",type="hist") 

## ----bbunny3, echo=T,eval=T,fig.width=10,fig.height=3,dependson="bbunny2"----

library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
customFromTxDb <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                  chromosome="chr17") 
plotTracks(c(accDT,customFromTxDb,genomeAxis), 
           from=47504051,to=47600688, 
           chromosome="chr17",type="hist") 

## ---- echo=TRUE----------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19) 
BSgenome.Hsapiens.UCSC.hg19[["chr7"]] 

## ---- echo=TRUE,fig.width=20,fig.height=3--------------------------------
sTrack <- SequenceTrack(Hsapiens) 
plotTracks(sTrack,from=134887024,to=134887074, 
           chromosome = "chr7",cex=2.5) 

## ---- echo=T,fig.width=20,fig.height=3-----------------------------------
dsSet <- DNAStringSet(Hsapiens[["chr7"]]) 
names(dsSet) <- "chr7" 
sTrack <- SequenceTrack(dsSet) 
plotTracks(sTrack,from=134887024,to=134887074, 
           chromosome = "chr7",cex=2.5) 

## ---- echo=F,eval=F------------------------------------------------------
## dsSet <- DNAStringSet(Hsapiens[["chr7"]])
## tempSet <- DNAStringSet(dsSet[[1]][134887024:134887074])
## names(tempSet) <- "chr7"
## writeXStringSet(tempSet,file="../../Data/chr7Short.fa")
## sTrack <- SequenceTrack("../../Data/chr7Short.fa")
## plotTracks(sTrack,from=1,to=50,
##            chromosome = "chr7")

## ---- echo=T,eval=T,fig.width=20,fig.height=3----------------------------
sTrack <- SequenceTrack("../../Data/chr7Short.fa") 
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7",cex=3) 

## ---- echo=T,eval=T,fig.width=20,fig.height=3----------------------------
sTrack <- SequenceTrack("../../Data/chr7Short.fa") 
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7",complement=T,cex=3) 

## ---- echo=T,eval=T,fig.width=20,fig.height=3----------------------------
sTrack <- SequenceTrack("../../Data/chr7Short.fa") 
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7",complement=F, 
           add53=T,cex=2.5) 

## ---- echo=T,eval=T,fig.width=20,fig.height=2----------------------------
sTrack <- SequenceTrack("../../Data/chr7Short.fa") 
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7",complement=T, 
           add53=T,cex=2.5) 

## ---- echo=T,collapse=T,fig.width=20,fig.height=2------------------------
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7",cex=2.5) 
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7", 
           cex=5) 

## ---- echo=T,collapse=T,fig.width=20,fig.height=2------------------------
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7",cex=1.1,noLetters=TRUE) 

## ---- echo=T,collapse=T,fig.width=20,fig.height=2------------------------
colForLetters <- c(A = "darkgrey", C = "red",
                    T = "darkgrey",G = "darkgrey")
plotTracks(sTrack,from=1,to=50, 
           chromosome = "chr7",cex=3,
           fontcolor=colForLetters) 

## ---- echo=T,fig.width=20,fig.height=5-----------------------------------
   peakReads <- AlignmentsTrack("../../Data/small_Sorted_SRR568129.bam") 
   peakReads 

## ---- echo=T,fig.width=20,fig.height=5-----------------------------------
   plotTracks(peakReads, 
              chromosome="chr5", 
              from=135312577, 
              to=135314146) 

## ---- echo=T,fig.width=20,fig.height=5-----------------------------------
   plotTracks(peakReads, 
              chromosome="chr5", 
              from=135312577, 
              to=135314146) 

## ---- echo=T,fig.width=20,fig.height=5-----------------------------------
   plotTracks(peakReads, 
              chromosome="chr5", 
              from=135312577, 
              to=135314146, 
              type="pileup") 

## ---- echo=T,fig.width=20,fig.height=5-----------------------------------
   plotTracks(peakReads, 
              chromosome="chr5", 
              from=135312577, 
              to=135314146, 
              type="coverage") 

## ---- echo=T,fig.width=20,fig.height=4-----------------------------------
   plotTracks(peakReads, 
              chromosome="chr5", 
              from=135312577, 
              to=135314146, 
              type=c("pileup","coverage")) 

## ----aa, echo=T,fig.width=20,fig.height=5,cache=T------------------------
 
heartReads <- AlignmentsTrack("../../Data/heart.bodyMap.bam", 
                           isPaired = TRUE) 
liverReads <- AlignmentsTrack("../../Data/liver.bodyMap.bam",  
                           isPaired = TRUE) 
 
liverReads 

## ----bb, echo=T,fig.width=20,fig.height=5,cache=T,dependson="aa"---------
plotTracks(c(heartReads,liverReads), 
           chromosome="chr12", 
           from=98986825,to=98997877) 
 

## ----cc, echo=T,fig.width=20,fig.height=5,cache=T,dependson="bb"---------
plotTracks(c(heartReads,liverReads), 
           chromosome="chr12", 
           from=98986825, 
           to=98997877, 
           type=c("coverage","sashimi")) 
 

## ----dd, echo=T,fig.width=18,fig.height=5,cache=T,dependson="cc"---------
plotTracks(c(liverReads), 
           chromosome="chr12", 
           from=98986825,to=98997877, 
           col.gap="Red",col.mate="Blue") 

## ----ee, echo=T,fig.width=18,fig.height=5,cache=T,dependson="dd"---------
plotTracks(c(liverReads), 
           chromosome="chr12", 
           from=98986825,to=98997877, 
           lty.gap=2,lty.mate=1) 
 

## ----aaas, echo=F,fig.width=18,fig.height=5,cache=T----------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
customFromTxDb <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                  chromosome="chr6") 

## ----ff, echo=T,fig.width=18,fig.height=5,cache=T,dependson="ee"---------
library(BSgenome.Mmusculus.UCSC.mm10)
sTrack <- SequenceTrack(BSgenome.Mmusculus.UCSC.mm10) 
activatedReads <- AlignmentsTrack("../../Data/activatedSNPread.bam", 
                           isPaired = TRUE, 
                           referenceSequence=sTrack) 

## ----gg, echo=T,fig.width=18,fig.height=5,cache=T,dependson=c("ff","aaas")----
plotTracks(c(activatedReads,customFromTxDb), 
           chromosome="chr6", 
           from=124815373,to=124815412) 
  

## ----hh, echo=T,fig.width=18,fig.height=5,cache=T,dependson="gg"---------
plotTracks(c(activatedReads,sTrack), 
           chromosome="chr6", 
           from=124815373,to=124815412,cex=2) 

## ----a, echo=T,fig.width=18,fig.height=5,cache=T-------------------------
bgrTrack <- BiomartGeneRegionTrack(genome="hg19", 
                                   start=26591341, 
                                   end=27034958, 
                                   chromosome = "chr7", 
                                   name="ENSEMBL") 

## ----b, echo=T,fig.width=18,fig.height=5,cache=T,dependson="a"-----------
plotTracks(bgrTrack) 

## ----c, echo=T,eval=FALSE,fig.width=10,fig.height=5,cache=T,dependson="b",warning=FALSE,message=FALSE----
## library(biomaRt)
## martList <- listMarts()
## mart = useMart("ENSEMBL_MART_ENSEMBL")
## dataList <- listDatasets(mart)
## mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
## filterList <- listFilters(mart)

## ----d, echo=T,eval=FALSE,fig.width=20,fig.height=5,cache=T,dependson="c"----
## bgrTrack <- BiomartGeneRegionTrack(genome="hg19",
##                                    start=26591341,
##                                    end=27034958,
##                                    chromosome = "chr7",
##                                    name="ENSEMBL",
##                                   filter=list(source="ensembl_havana"))

## ----dd1, echo=F,eval=T,cache=F------------------------------------------
load("../../Data/ensembl_havana_encode.RData") 

## ----e, echo=T,eval=FALSE,fig.width=18,fig.height=5,cache=F--------------
## plotTracks(bgrTrack)

## ----a1, echo=T,fig.width=10,fig.height=5,cache=T------------------------
library(rtracklayer) 
session <- browserSession() 
genome(session) <- "hg19" 
trackNames(session) 

## ----b1, echo=T,fig.width=10,fig.height=5,cache=T,dependson="a1"---------
query <- ucscTableQuery(session, "Ensembl Genes", 
                        GRangesForUCSCGenome("hg19", "chr7", 
                                             IRanges(26591341,27034958))) 
tableNames(query) 


## ----c1, echo=T,eval=F,fig.width=10,fig.height=5,cache=T,dependson="b1"----
## ucscTrack <- UcscTrack(genome = "hg19",
##                        chromosome = "chr7",
##                        track = "ensGene",
##                        from = 26591341,
##                        to = 27034958,
##                        trackType = "GeneRegionTrack",
##                        rstarts = "exonStarts",
##                        rends = "exonEnds",
##                        gene ="name",
##                        symbol = "name2",
##                        transcript = "name",
##                        strand = "strand"
## )
## 

## ----d1, echo=F,eval=T,cache=T,dependson="c1"----------------------------
load("../../Data/ensGene_UCSC.RData") 

## ----e1, echo=T,fig.width=18,fig.height=5,cache=T,dependson="d1"---------
plotTracks(c(bgrTrack,ucscTrack), 
           from = 26591341,to = 27034958) 

## ----f1, eval=F,echo=T,fig.width=10,fig.height=5,cache=T,dependson="e1"----
## conservationTrack <- UcscTrack(genome = "hg19", chromosome = "chr5",track = "Conservation", table = "phyloP100wayAll",from = 135313003, to = 135313570, trackType = "DataTrack",start = "start", end = "end", data = "score",type = "hist", window = "auto", col.histogram = "darkblue",fill.histogram = "darkblue", ylim = c(-3.7, 4),name = "Conservation")
## 

## ----f11, eval=F,echo=T,fig.width=10,fig.height=5,cache=T,dependson="e1"----
## conservationTrack <- UcscTrack(genome = "hg19",
##                                chromosome = "chr5",
##                                track = "Conservation",
##                                table = "phyloP100wayAll",
##                                from = 135313003,to = 135313570,
##                                trackType = "DataTrack",
##                                start = "start", end = "end",
##                                data = "score",type = "hist", window = "auto",
##                                col.histogram = "darkblue",fill.histogram = "darkblue",
##                                ylim = c(-3.7, 4),name = "Conservation")
## 

## ----g1, echo=F,eval=T,cache=T,dependson="f1"----------------------------
load("../../Data/conservation.RData") 

## ----h1, echo=F,cache=T,dependson="g1"-----------------------------------
genomeAxis <- GenomeAxisTrack(name="MyAxis",scale=250) 

## ----i1, echo=T,fig.width=18,fig.height=5,cache=T,dependson="h1"---------
plotTracks(c(conservationTrack,peakReads,genomeAxis), 
           from=135313003, 
           to=135313570, 
           chromosome = "chr5", 
           type = c("hist","coverage"), 
           sizes = c(1,1,0.2), 
           cex=2) 

