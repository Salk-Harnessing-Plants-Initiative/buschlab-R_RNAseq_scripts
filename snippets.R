library(DESeq2)

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
dds

# takes approx. 450 sec
system.time(dds <- DESeq(dds))

# takes approx. 390 sec 
system.time(rld <- rlogTransformation(dds))
head(assay(rld))
hist(assay(rld))