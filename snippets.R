library(DESeq2)
library(BiocParallel)
library(parallel)

# here starting with provided read counts matrix
# adopt to your needs
inData <- read.table("/Users/alexander.bindeus/Desktop/RNAseq/GSE70094_CountTable_HTSeq_At_raw.txt",header=TRUE)

# removing Ci treated ones
# inData <- inData[,!(colnames(inData) %in% c("minusP10Ci_Exp1_1","minusP10Ci_Exp1_2","minusP10Ci_Exp1_3"))]
# inData <- inData[,!(colnames(inData) %in% c("minusP10Ci_Exp2_1","minusP10Ci_Exp2_2","minusP10Ci_Exp2_3"))]
# inData <- inData[,!(colnames(inData) %in% c("minusP24Ci_Exp2_1","minusP24Ci_Exp2_2","minusP24Ci_Exp2_3"))]

# get colnames to use as conditions
libs <- colnames(inData)
libs <- substr(libs,1,nchar(libs)-2)
conditions <- as.factor(libs) 

sample<-as.matrix(inData)

coldata <- data.frame(row.names=colnames(sample), conditions)
dds <- DESeqDataSetFromMatrix(countData=sample, colData=coldata, design=~conditions)

numCores <- detectCores(logical = TRUE) # use all available cores
register(MulticoreParam(numCores))
system.time(dds <- DESeq(dds, parallel=TRUE)) # now got mu and cooks 

# alpha = padj (FDR) threshold (def: 0.1), lfcThreshold = log2 fold change threshold (def. 0)
# system.time(res <- results(dds, alpha = 0.05, lfcThreshold = 0, parallel=TRUE))
# summary(res)

system.time(rld <- rlogTransformation(dds))
head(assay(rld))
hist(assay(rld))

##########################

library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(conditions))]

###### Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[conditions], RowSideColors=mycols[conditions],
          margin=c(10, 10), main="Sample Distance Matrix")

require(ggplot2)
library(d3heatmap)

d3heatmap(sampleDists, #forheatmap,
          #cellnote = forheatmap,  # same data set for cell labels
          main = "CSI", # heat map title
          notecol="black",      # change font color of cell labels to black
          #density.info="none",  # turns off density plot inside color legend
          key = TRUE,
          #labCol = emptylabels, # myylabels,
          #labRow = emptylabels, #myylabels,
          #labRow = myylabels,
          na.rm = TRUE,
          symkey=FALSE,         # colorkey not sym around 0
          revC = TRUE,
          colors = "Blues",
          dendrogram="none",     # only draw a row dendrogram
          srtCol = 45,
          height = 900,
          width = 900
)
#############################

#############################

myPCA <- plotPCA(rld, intgroup="conditions",returnData=TRUE)
percentVar <- round(100*attr(myPCA,"percentVar"))
tooltips<-paste(myPCA$name)
library(scatterD3)

finCol <- c(rep(mycols,3))

scatterD3(x=myPCA$PC1, y=myPCA$PC2,xlab=paste0("PC1: ",percentVar[1],"% variance"),ylab = paste0("PC2: ", percentVar[2],"% variance"),
         tooltip_text = tooltips,
         #col_var = mycols
         lab = myPCA$name
          )

############################

res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="/Users/alexander.bindeus/Desktop/RNAseq/deseq_res_all_2505.csv")

### counts plot ####
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="conditions",
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=conditions, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))
  
############################

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
maplot(resdata, main="MA Plot")

############################

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
#png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))

#############################
