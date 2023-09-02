options(mc.cores = 4L)
graphics.off()
library(tximport)
library(DESeq2)
library(tidyverse)
#library(readr)
#library(dplyr)
#library(magrittr)
#library(ggplot2)
library(ggpubr)
library(biomaRt)
library(data.table)
library(pheatmap)

setwd("C:/Users/Evangelos Mourkas/OneDrive - Nexus365/Documents/Oxford_Files/12.GitHub_preparation_scripts/2.Differential_expression_using_DeSeq2/2.Neisseria")

samples <- read_csv("RNAseq_samples.csv", show_col_types = FALSE)
head(samples)

gene_map <- read.csv("gene_map2.csv", col.names = c('Names','Genes'))
head(gene_map)

# Import a file that specifies which transcripts are associated with which genes. Skip if transcripts already have gene information #
# See salmon transcriptome index #
# tx2gene_map <- read_tsv("path/to/.tsv") #
# head(tx2gene_map) #

samples = samples[1:6,]
samples$condition = factor(rep(c('Sample_1', 'Sample_2'), each =3), levels = c('Sample_1', 'Sample_2'))
# Make a vector list with the quant.sf files from the different samples #
# Make the tx2gene file manually for bacteria #
txi <- tximport(files = samples$quant_file, type = "salmon", tx2gene = gene_map, ignoreTxVersion = TRUE)
summary(txi)
# Count data #
# Number of reads in each of these samples that match to transcripts that belong to these genes #
txi[['counts']]
head(txi$counts)
colnames(txi$counts) <- samples$sample


# Normalization and transformation of read counts #
# DeSeq2 will: 1. model the raw counts, using normalization factors (size factors) to account for differences in library depth
# 2. estimate the gene-wise dispersion and shrink these estimates to generate more accurate estimates of dispersion to model the counts
# 3. fit the negative binomial model and perform hypothesis testing using the Wald test or Likelihood Ratio Test

# Transform txi object into DeSeq2 readable file #
# A design formula tells the statistical software the known sources of variation to control for, as well as the factor of interest to test for during differential expression testing #
samples = as.data.frame(samples)
condition = c('Sample_1', 'Sample_2')
condition = rep(condition, each=3)
condition = factor(condition)

### 3 steps to DESeq2 analysis ###
## 1) estimate size factors (normalization process) ##

dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = samples, 
                                design = ~ condition)

# Check normalization after Deseq analysis #
counts(dds)[1:6, 1:6]
txi$counts[1:6, 1:6]
boxplot(txi$counts)


dds = estimateSizeFactors(dds)
normalizationFactors(dds)
counts(dds, normalized = TRUE)[1:6,1:6]
boxplot(counts(dds, normalized=TRUE))

# Variance stabilized transformation on the count data, while controlling for library size of samples #
# logarithmic normalization #
vst <- varianceStabilizingTransformation(dds)
boxplot(assay(vst))

# Apply PCA analysis to the data #
plotPCA(vst, intgroup= 'condition') + theme_bw()

# Gives back those variance stabilized transformed counts
# Distance matrix calculation using euclidean (best method) between our samples
d = assay(vst)
d = t(d)
d = dist(d)

# Clustering using complete clustering
h = hclust(d)
plot(h)

k = kmeans(t(assay(vst)), centers = 2)

## 2) estimate dispersions of variance - distribution to see whether genes are differential expressed ##
# Uses variance across the whole experiment to submit disperse estimates to shrinkage #
# imperial bayesian shrinkage #

dds = estimateDispersions(dds)
#png("Dispersion.png", width = 1500, height = 1300)
plotDispEsts(dds)
#dev.off()
# Notes: Mean and variance in RNASeq data are not independent from one another. 
# The variance is not continuous with the mean of the count data
# As the mean of normalized counts continuous your variance goes down, so your dispersion reduces as your mean counts goes up #
# The black points (raw estimate of dispersion) are the genewise estimates of dispersion (for each gene for the 6 observations the dispersion has been estimated)
# Mean of the counts and the mean of the dispersion (mean dispersion - red points)
# THe bayesian shrinkage estimation is applied - the black points are shrank towards that red fitted disperse estimates (blue points)

## 3) apply statistics (Wald Test) across the negative binomial distributions
dds = nbinomWaldTest(dds)


## apply all these 3 in one line
# dds = DESeq(dds)

# Extract information from our data
results_table = results(dds, alpha = 0.05)
summary(results_table)
View(as.data.frame(results_table))


# The reference condition is defined by alphabetical order unless you set it up as a level after the factor (levels=c(',')) #
# results_table is a DESeq DataFrame not a data.frame!
# make a copy of results_table as a data.frame format
results_df = as.data.frame(results_table)

# Notes: baseMean - the average normalized count for each gene across all samples in the experiment
# if it's 0 then the gene is not expressed. When is it below cutoff (mean count < 0) then they are filtered out from multiple correction testing, 
# some genes are also filtered because they are outliers (failed the cutoff filter - test for whether a gene has an outlier in each set of observations)
# gene example
#setwd("C:/Users/Evangelos/Documents/BATH_UNI/9.RNA_Expression/2.Test2/DeSeq2/1.N188_1_VS_N222_1")

#png("outlier_example.png", width = 1500, height = 1300)
plotCounts(dds, gene='abcZ', intgroup = 'condition')
#dev.off()

# Remove those outliers and filtered genes out of our dataset
sum(complete.cases(results_df))
filter_df1 = results_df[complete.cases(results_df),]
dim(filter_df1)
View(filter_df1)

# log2FoldChange - Difference between the two conditions in our differentially expressed test as a log ratio
# We use log2 because: fold changes are quoted in log space
# Practical reasons: 
# lfcSE: standard error of the log2foldchange
# Wald test: 1. stat - pvalue (raw pvalue) the raw output of this statistical test
# 2. padj - (adjusted p-value) - those tests are dependend on one another for each gene. Account for the 5% of our observations being wrong multiplied for 1709 genes.  

# Filter results with padj < 0.05 , log2FoldChange > 1 < -1
filter_df1$padj < 0.05
filter_df2 = filter_df1[filter_df1$padj < 0.05,]
dim(filter_df2)

# controlling for effect size because of low number of replicates, so adding log2foldchange filter sequentially
# turn all values to absolute values and then filter out the ones outside the magnitute of log2foldchange < -1 and > 1
abs(filter_df2$log2FoldChange) > 1
filter_df3 = filter_df2[abs(filter_df2$log2FoldChange) > 1,]
dim(filter_df3)

View(filter_df3)

### VISUALIZATION ###
# 1. DESEQ2 visualize all of the genes #
#png("dds_MA.png", width = 1500, height = 1300)
plotMA(dds)
#dev.off()
# Comparison of the mean normalized count (basemean) and fold change
# Illustrates the high level of variance in low counts. Genes with low counts tend to high fold change. Adjusted pvalue < 0.1 are blue (we filtered for lower)
# Genes of low count do not get to be robusted. The noise

# Volcano plot
# The log of a number between 1 and 0 then you get a negative number. Log base 10 negative. so transforming into a positive.
# The lowest the pvalue the higher up the point will be.

filter_df1$test = filter_df1$padj < 0.05 & abs(filter_df1$log2FoldChange) > 1
# converting rownames as an extra column for the heatmap
test_df <- tibble::rownames_to_column(filter_df1, "Gene")
write.csv(as.data.frame(test_df), "filtered_results.csv")

#png("volcano_plot.png", width = 1500, height = 1300)
ggplot(test_df, aes(x=log2FoldChange, y=-log10(padj), name = Gene)) + 
  geom_point(aes(colour=test), size = 1) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept = 1, colour = 'darkgreen', linetype = 2) + 
  geom_vline(xintercept = -1, colour = 'darkgreen', linetype = 2) + 
  geom_hline(yintercept = -log10(0.05), colour = 'blue', linetype = 2) + 
  theme_bw() + 
  theme(legend.position = 'top')
dev.off()

#png("ggmaplot.png", width = 1500, height = 1300)
ggmaplot(filter_df1, fdr = 0.05, fc = 2, genenames= NULL, size = 1, alpha = 1, 
         font.label = c(9, "bold"), palette = c("#B31B21", "#1454AC", "darkgray"),
         top = 30, select.top.method = c("padj", "fc"), label.rectangle = TRUE, label.select = NULL, main = "N188.1 vs N222.1", xlab = "Log2 mean expression",
         ylab = "Log2 fold change", ggtheme = theme_classic(), legend = "top", font.main= "bold", font.legend = "bold",
)
dev.off()


anno_df3 <- tibble::rownames_to_column(filter_df3, "Gene")
dim(anno_df3)

# Set up differentially expressed genes - degs
degs = anno_df3$Gene

vst_mat = assay(vst)

data_for_hm = vst_mat[degs,]
dim(data_for_hm)


#generate a quick and dirty heatmap
heatmap(data_for_hm)

#Generate a good heatmap
#png("pheatmap.png", width = 1500, height = 1300)
pheatmap(data_for_hm, fontsize_row = 4, scale = 'row', cutree_cols = 2, cutree_rows = 2, border_color = NA)
dev.off()
#Upload a dataset that contains the top 20 genes differentially expressed genes with the lowest padj values



