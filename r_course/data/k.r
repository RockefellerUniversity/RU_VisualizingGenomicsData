

---
  ## Plotting ChIPseq profiles.
  
  
  ```{r gD,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1",warning=FALSE,message=FALSE}
library(soGGi)
```


---
  ## Plotting ChIPseq profiles.
  
  
  ```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
myRes <- read.delim("../../Data/Antibody_CTCF___Group_CTCF_Ch12_minus_CTCF_MelDEG.xls",sep="\t")
dvev <- matrix(unlist(strsplit(as.vector(myRes[,1]),"_")),ncol=4,byrow=T)
CTCF_Ch12MinusMel <- GRanges(seqnames=dvev[,2],IRanges(as.numeric(dvev[,3]),
                                                       as.numeric(dvev[,4])))
mcols(CTCF_Ch12MinusMel) <- as.data.frame(myRes[,-c(1,2,9:16)])

# load("../../Data/dds.RData")
# library(DESeq2)
# myRes <- results(dds,contrast = c("Group","CTCF_Ch12","CTCF_Mel"),
#                  format = "GRanges")

```


---
  ## Plotting ChIPseq profiles.
  
  
  ```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
CTCF_Ch12MinusMel <- CTCF_Ch12MinusMel[order(CTCF_Ch12MinusMel$baseMean,decreasing = TRUE),]
# load("../../Data/dds.RData")
# library(DESeq2)
# myRes <- results(dds,contrast = c("Group","CTCF_Ch12","CTCF_Mel"),
#                  format = "GRanges")
CTCF_ch12 <- regionPlot(bamFile = "../../Data/Sorted_CTCF_Ch12_1Normalised.bw",testRanges = CTCF_Ch12MinusMel,
                        style = "point",format = "bigwig")
CTCF_ch12 <- CTCF_profile
save(CTCF_ch12,file="../../Data/CTCF_ch12.RData")
load(file="../../Data/CTCF_ch12.RData")
plotRegion(CTCF_ch12)
```



```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
CTCF_mel <- regionPlot(bamFile = "../../Data/Sorted_CTCF_MEL_1Normalised.bw",testRanges = CTCF_Ch12MinusMel,
                       style = "point",format = "bigwig")
plotRegion(CTCF_mel)
CTCF_all <- c(CTCF_mel,CTCF_ch12)
plotRegion(CTCF_all,groupBy = "Sample")

save(CTCF_all,file="../../Data/CTCF_ch12.RData")


```


```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}

plotRegion(CTCF_all,groupBy = "Sample",
           colourBy = "Sample")


```



```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}

plotRegion(CTCF_all,
           colourBy = "Sample")


```


```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}

MelUp <- CTCF_Ch12MinusMel[CTCF_Ch12MinusMel$padj < 0.05 & 
                             !is.na(CTCF_Ch12MinusMel$padj) & 
                             CTCF_Ch12MinusMel$log2FoldChange < -4]
Ch12Up <- CTCF_Ch12MinusMel[CTCF_Ch12MinusMel$padj < 0.05 & 
                              !is.na(CTCF_Ch12MinusMel$padj) & 
                              CTCF_Ch12MinusMel$log2FoldChange > 4]
plotRegion(CTCF_all,gts = list(MelUp=MelUp,Ch12Up=Ch12Up),
           groupBy = "Group",colourBy="Sample")


```


```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
myTop <- CTCF_Ch12MinusMel[order(CTCF_Ch12MinusMel$stat),][1:100]
CTCF_Mel_Up <- CTCF_all[rowRanges(CTCF_all) %over% myTop,]
k1 <- plotHeatmap(CTCF_Mel_Up[[1]])
k2 <- plotHeatmap(CTCF_Mel_Up[[2]])


CTCF_Mel_Up <- CTCF_all[rowRanges(CTCF_all) %over% MelUp,]

k1 <- plotHeatmap(CTCF_Mel_Up[[1]],rowScale = FALSE)
k2 <- plotHeatmap(CTCF_Mel_Up[[2]],rowScale = FALSE)


temp <- cbind(k1,k2)
pheatmap::pheatmap(temp,
                   cluster_rows = FALSE,cluster_cols =FALSE,scale="row",
                   ,useRaster=TRUE,gaps_col = 100)
```



---
  
  
  ##Exercises 
  
  
  Time for exercises! [Link here](../../Exercises/Viz_part3_exercises.html) 

---
  ##Solutions 
  
  
  Time for solutions! [Link here](../../Answers/Viz_part3_answers.html) 





---
  ## Plotting ChIPseq profiles.
  
  
  ```{r gD,eval=TRUE,echo=TRUE,cache=TRUE,dependson="gC1",warning=FALSE,message=FALSE}
library(soGGi)
```


---
  ## Plotting ChIPseq profiles.
  
  
  ```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
myRes <- read.delim("../../Data/Antibody_CTCF___Group_CTCF_Ch12_minus_CTCF_MelDEG.xls",sep="\t")
dvev <- matrix(unlist(strsplit(as.vector(myRes[,1]),"_")),ncol=4,byrow=T)
CTCF_Ch12MinusMel <- GRanges(seqnames=dvev[,2],IRanges(as.numeric(dvev[,3]),
                                                       as.numeric(dvev[,4])))
mcols(CTCF_Ch12MinusMel) <- as.data.frame(myRes[,-c(1,2,9:16)])

# load("../../Data/dds.RData")
# library(DESeq2)
# myRes <- results(dds,contrast = c("Group","CTCF_Ch12","CTCF_Mel"),
#                  format = "GRanges")

```


---
  ## Plotting ChIPseq profiles.
  
  
  ```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
CTCF_Ch12MinusMel <- CTCF_Ch12MinusMel[order(CTCF_Ch12MinusMel$baseMean,decreasing = TRUE),]
# load("../../Data/dds.RData")
# library(DESeq2)
# myRes <- results(dds,contrast = c("Group","CTCF_Ch12","CTCF_Mel"),
#                  format = "GRanges")
CTCF_ch12 <- regionPlot(bamFile = "../../Data/Sorted_CTCF_Ch12_1Normalised.bw",testRanges = CTCF_Ch12MinusMel,
                        style = "point",format = "bigwig")
CTCF_ch12 <- CTCF_profile
save(CTCF_ch12,file="../../Data/CTCF_ch12.RData")
load(file="../../Data/CTCF_ch12.RData")
plotRegion(CTCF_ch12)
```



```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
CTCF_mel <- regionPlot(bamFile = "../../Data/Sorted_CTCF_MEL_1Normalised.bw",testRanges = CTCF_Ch12MinusMel,
                       style = "point",format = "bigwig")
plotRegion(CTCF_mel)
CTCF_all <- c(CTCF_mel,CTCF_ch12)
plotRegion(CTCF_all,groupBy = "Sample")

save(CTCF_all,file="../../Data/CTCF_ch12.RData")


```


```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}

plotRegion(CTCF_all,groupBy = "Sample",
           colourBy = "Sample")


```



```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}

plotRegion(CTCF_all,
           colourBy = "Sample")


```


```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}

MelUp <- CTCF_Ch12MinusMel[CTCF_Ch12MinusMel$padj < 0.05 & 
                             !is.na(CTCF_Ch12MinusMel$padj) & 
                             CTCF_Ch12MinusMel$log2FoldChange < -4]
Ch12Up <- CTCF_Ch12MinusMel[CTCF_Ch12MinusMel$padj < 0.05 & 
                              !is.na(CTCF_Ch12MinusMel$padj) & 
                              CTCF_Ch12MinusMel$log2FoldChange > 4]
plotRegion(CTCF_all,gts = list(MelUp=MelUp,Ch12Up=Ch12Up),
           groupBy = "Group",colourBy="Sample")


```


```{r y_a,eval=TRUE,echo=TRUE,cache=TRUE,warning=FALSE,message=FALSE}
myTop <- CTCF_Ch12MinusMel[order(CTCF_Ch12MinusMel$stat),][1:100]
CTCF_Mel_Up <- CTCF_all[rowRanges(CTCF_all) %over% myTop,]
k1 <- plotHeatmap(CTCF_Mel_Up[[1]])
k2 <- plotHeatmap(CTCF_Mel_Up[[2]])


CTCF_Mel_Up <- CTCF_all[rowRanges(CTCF_all) %over% MelUp,]

k1 <- plotHeatmap(CTCF_Mel_Up[[1]],rowScale = FALSE)
k2 <- plotHeatmap(CTCF_Mel_Up[[2]],rowScale = FALSE)


temp <- cbind(k1,k2)
pheatmap::pheatmap(temp,
                   cluster_rows = FALSE,cluster_cols =FALSE,scale="row",
                   ,useRaster=TRUE,gaps_col = 100)



profile=CTCF_Mel_Up[[1]]
bins=100
col=colorRampPalette(blues9)(100)
rowScale=FALSE
orderPosition=NULL
orderBy="maxAtPosition"
matt <- assay(profile)
dwdw <- apply(matt,1,function(x)any(is.na(x)))
matt <- matt[!dwdw,]
if(is.null(orderPosition)){
  if(is.null(bins)){
    orderPosition <- unique(c(floor(bins/2),ceiling(bins/2)))
  }else{
    orderPosition <- unique(c(floor(bins/2),ceiling(bins/2)))
  }
}
# if(rowScale == TRUE){
#   cols <- colorRampPalette(brewer.pal(9,"Blues"),bias=1)(100)
# }else{
#   cols <- colorRampPalette(brewer.pal(9,"Blues"),bias=10)(100)
# }

if(!is.null(bins)){
  binsize <- floor(ncol(assay(profile))/bins)
  binremainner <- ncol(matt)%%bins
  mat <- matrix(nrow=nrow(matt),ncol=bins)
  firstIndex <- 0
  endIndex <- floor(binsize/2)+binsize
  mat[,1] <- rowMeans(matt[,firstIndex:endIndex])
  for(i in 2:(bins-1)){
    firstIndex <- endIndex+1
    endIndex <- endIndex+binsize
    mat[,i] <- rowMeans(matt[,firstIndex:endIndex])
  }
  
  mat[,i+1] <- rowMeans(matt[,endIndex:ncol(matt)])
  matt <- mat
}else{
  message("No binning of matrix done")
}
if(orderBy=="maxAtPosition"){
  if(length(orderPosition) == 1){
    matt <- matt[order(matt[,orderPosition],decreasing=TRUE),]
  }else{
    matt <- matt[order(rowMeans(matt[,min(orderPosition):max(orderPosition)]),decreasing=TRUE),]
  } 
}
if(rowScale==TRUE){
  matt <- t(scale(t(matt),center=TRUE,scale=TRUE))
}

qs <- quantile(matt,c(0,0.1,0.9,1))  
breaks=seq(qs[2],qs[3],length.out = round(length(col)/2))
breaks2=seq(qs[3],qs[4],length.out = round(length(col)/2)+1)
breaks <- c(breaks,breaks2)
layout(matrix(data=c(1,2), nrow=1, ncol=2),
       widths=c(4,1), heights=c(1,1))
matttoPlot <- matt[rev(1:nrow(matt)),]
image(t(matttoPlot),useRaster=TRUE,
      xaxt='n',yaxt="n",col=col,breaks=breaks)

par(mar = c(3,2.5,2.5,2))
image(1, 1:length(col),
      matrix(data=seq(min(matt,na.rm = T),max(matt,na.rm=T),
                      length=length(col)),ncol=length(col),
             nrow=1),
      xlab="",ylab="",xaxt="n",las=2,col=col)
return(matt)
```



---
  
  
  ##Exercises 
  
  
  Time for exercises! [Link here](../../Exercises/Viz_part3_exercises.html) 

---
  ##Solutions 
  
  
  Time for solutions! [Link here](../../Answers/Viz_part3_answers.html) 





library(soGGi)
library(rtracklayer)

library(soGGi)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
dsds <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
names(dsds) <- NULL
csc41 <- regionPlot(bamFile = "~/Downloads/ENCFF239XXP.bigWig",testRanges = dsds,style = "percentOfRegion",format = "bigwig")
plotRegion(csc41)

polSer2 <- csc41
save(polSer2,file="../RU_VisualizingGenomicsData2/viz_course/Data/polSer2.RData")

plotRegion(csc41,outliers = 0.01)
library(RColorBrewer)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(blueColours)(50)
matt <- plotHeatmap(csc41,col = colors,orderPosition = c(60:70))
dwdw <- apply(matt,1,function(x)any(is.na(x)))
myMatt <- matt[!dwdw,]
plot <- pheatmap::pheatmap(myMatt,scale="none",kmeans_k = 2,cluster_cols = FALSE,useRaster=TRUE)
library(gplots)
library(pheatmap) 

dwd <-pheatmap(myMatt,cluster_rows = FALSE,cluster_cols =FALSE,scale="none",
               ,useRaster=TRUE,kmeans_k = 4)

pdf("../RU_VisualizingGenomicsData2/viz_course/imgs/pol2Ser2_raster.pdf")
heatmap.2(myMatt[order(dwd$kmeans$cluster),], Rowv=FALSE, Colv=FALSE, useRaster=TRUE, labRow=NA, labCol=NA,symbreaks=TRUE, breaks=51,
          trace="none", key=FALSE)
dev.off()



download.file("http://simonsoftware.se/other/xkcd.ttf",dest="xkcd.ttf", mode="wb")
system("mkdir ~/.fonts")
system("cp xkcd.ttf  ~/.fonts")
font_import(paths = " ~/.fonts",pattern = "xkcd", prompt=TRUE)
fonts()
fonttable()
if(.Platform$OS.type != "unix") {
 ## Register fonts for Windows bitmap output
loadfonts(device="win")
 } else {
  loadfonts()
}




