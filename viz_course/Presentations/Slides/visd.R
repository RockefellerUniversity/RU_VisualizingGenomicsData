

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


