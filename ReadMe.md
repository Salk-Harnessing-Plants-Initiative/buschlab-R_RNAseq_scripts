collection of various R snippets for RNA diffseq analyses

terminology:

FDR:  False discovery rate

GLM:  Generalized linear model

HTS:  High-throughput sequencing

LFC:  Logarithmic fold change

MAP:  Maximum a posteriori

MLE:  Maximum-likelihood estimate

RNA-seq:  RNA sequencing

VST:  Variance-stabilizing transformation

BCV:  biological coefficient of variation 

logCPM: log2 counts per million, normalized for library sizes

PCA:  principal components analysis 

lfcSE: log2 fold change Standard Error
stat: Wald test stat

The lfcSE gives the standard error of the log2FoldChange. For the Wald test, stat is the Wald statistic: 
the log2FoldChange divided by lfcSE, which is compared to a standard Normal distribution to generate
a two-tailed pvalue. For the likelihood ratio test (LRT), stat is the difference in deviance between
the reduced model and the full model, which is compared to a chi-squared distribution to generate a pvalue.

baseMean = row means of normalized counts

baseMean: The base mean (i.e., mean of the counts divided by the size factors) for the counts for both conditions

foldChange: ratio meanB/meanA

log2FoldChange: log2 of the fold change

pval: p value for rejecting the null hypothesis 'meanA==meanB'

padj: adjusted p values (adjusted with 'p.adjust( pval, method="BH")')

DEsequ2: function DESeq runs the following functions in order:
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)  (prior just nbinomTest)

ad estimateSizeFactors:
"Each column is divided by the geometric means of the rows. The median (or, ir requested, another location estimator)
of these ratios (skipping the genes with a geometric mean of zero) is used as the size factor for this column."
