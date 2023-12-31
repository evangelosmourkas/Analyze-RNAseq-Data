# Analyze RNAseq data using DESeq2
This repository contains the required script to analyze RNAseq data, with a focus on differential gene expression using DESeq2.
The script has been provided as an Rscript and can be run locally in Rstudio or Visual Studio Core.

## Files Summary
The files include data for two strains (three replicates for each strain), each strain representing a different condition. Briefly, illumina sequence reads were mapped to a reference transcriptome index using Salmon. Gene expression was calculated in units of transcripts per million (TPM) by Salmon, a relative abundance measure used for downstream analysis with DESeq2. Each directory has the following:
* **_Salmon_TPM_data_**: directory with data from the two strains used in the analysis
* **_Differential_gene_expression.R_**: rscript for the analysis
* **_RNAseq_samples.csv_**: .csv file containing the strain setup information used as input
* **_gene_map.csv_**: .csv file containing gene information as inferred by the reference transcriptome to be used as input

# Differential Gene Expression script
## R packages
* DESeq2
* ggplot2
* tximport
* tidyverse
* ggpubr
* biomaRt
* data.table
* pheatmap

## Set working directory
```
setwd("/Path/to/working/directory")
```
## Import and read files
```
samples <- read_csv("RNAseq_samples.csv", show_col_types = FALSE)
head(samples)

gene_map <- read.csv("gene_map2.csv", col.names = c('Names','Genes'))
head(gene_map)
```
## Specify samples and conditions
```
samples = samples[1:6,]
samples$condition = factor(rep(c('Sample_1', 'Sample_2'), each =3), levels = c('Sample_1', 'Sample_2'))
```
Note: The reference condition is defined by alphabetical order unless you set it up as a level after the factor (levels=c(',')).
## Import the tx2gene file with gene names information
```
txi <- tximport(files = samples$quant_file, type = "salmon", tx2gene = gene_map, ignoreTxVersion = TRUE)
summary(txi)
```
## Count data
```
txi[['counts']]
head(txi$counts)
colnames(txi$counts) <- samples$sample
```
## Normalization and transformation of read counts
DeSeq2 will: 
* Model the raw counts, using normalization factors (size factors) to account for differences in library depth
* Estimate the gene-wise dispersion and shrink these estimates to generate more accurate estimates of dispersion to model the counts
* Fit the negative binomial model and perform hypothesis testing using the Wald test or Likelihood Ratio Test

## Transform txi object into DeSeq2 readable file
```
samples = as.data.frame(samples)
condition = c('Sample_1', 'Sample_2')
condition = rep(condition, each=3)
condition = factor(condition)
```
## Normalization - estimate size factors
```
dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = samples, 
                                design = ~ condition)
```
## Check counts
```
counts(dds)[1:6, 1:6]
txi$counts[1:6, 1:6]
boxplot(txi$counts)
```
![boxplots](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/97932d71-5610-4b45-b288-bdb0632be0ab)
## Check counts after normalization
```
dds = estimateSizeFactors(dds)
normalizationFactors(dds)
counts(dds, normalized = TRUE)[1:6,1:6]
boxplot(counts(dds, normalized=TRUE))
```
![boxplots_norm](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/560b96ca-aad0-43af-aa21-faf45bd84b67)
## Variance stabilized transformation on the count data, while controlling for library size of samples
```
vst <- varianceStabilizingTransformation(dds)
boxplot(assay(vst))
```
![boxplots_lognorm](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/b5268289-e795-446a-bf08-bb7783d0a985)
## PCA analysis 
```
plotPCA(vst, intgroup= 'condition') + theme_bw()
```
![PCA](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/61a3ab53-3776-4d76-8a21-0a3b61f58224)
## Distance matrix calculation using euclidean (best method) between samples
```
d = assay(vst)
d = t(d)
d = dist(d)
```
## Clustering using complete clustering
```
h = hclust(d)
plot(h)
k = kmeans(t(assay(vst)), centers = 2)
```
![hclustering](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/38f1e210-db43-4730-ab44-7c9c72c95812)
## Estimate dispersions of variance - distribution to see whether genes are differential expressed
Uses variance across the whole experiment to submit disperse estimates to shrinkage - imperial bayesian shrinkage
```
dds = estimateDispersions(dds)
plotDispEsts(dds)
```
![Dispersion](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/e56b2947-dd30-4618-9a65-751bd83e7e7f)
Note: Mean and variance in RNASeq data are not independent from one another. The variance is not continuous with the mean of the count data. As the mean of normalized counts continuous, the variance goes down, so the dispersion reduces as the mean counts goes up. The black points (raw estimate of dispersion) are the genewise estimates of dispersion (for each gene for the 6 observations the dispersion has been estimated). Mean of the counts and the mean of the dispersion (mean dispersion - red points). The bayesian shrinkage estimation is applied - the black points are shrank towards that red fitted disperse estimates (blue points).
## Apply Wald test across the negative binomial distributions
```
dds = nbinomWaldTest(dds)
```
## Extracting information from the data
```
results_table = results(dds, alpha = 0.05)
summary(results_table)
View(as.data.frame(results_table))
results_df = as.data.frame(results_table)
```
## Removing outliers and filtered genes out of the dataset
```
sum(complete.cases(results_df))
filter_df1 = results_df[complete.cases(results_df),]
dim(filter_df1)
View(filter_df1)
```
Note: baseMean - the average normalized count for each gene across all samples in the experiment. If it's 0 then the gene is not expressed. When is it below cutoff (mean count < 0) then they are filtered out from multiple correction testing. Some genes are also filtered because they are outliers (failed the cutoff filter - test for whether a gene has an outlier in each set of observations).
## Visualize counts for individual genes - example
```
plotCounts(dds, gene='NEIS1437', intgroup = 'condition')
```
![Gene_example](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/915d8bd9-a8bd-442d-bea5-025a7dbba192)
## Filter results with padj < 0.05 , log2FoldChange > 1 < -1
```
filter_df1$test = filter_df1$padj < 0.05 & abs(filter_df1$log2FoldChange) > 1
```
## Converting rownames as an extra column for the heatmap
```
test_df <- tibble::rownames_to_column(filter_df1, "Gene")
write.csv(as.data.frame(test_df), "filtered_results.csv")
```
# Visualization 
## Volcano plot
```
ggplot(test_df, aes(x=log2FoldChange, y=-log10(padj), name = Gene)) + 
  geom_point(aes(colour=test), size = 1) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept = 1, colour = 'darkgreen', linetype = 2) + 
  geom_vline(xintercept = -1, colour = 'darkgreen', linetype = 2) + 
  geom_hline(yintercept = -log10(0.05), colour = 'blue', linetype = 2) + 
  theme_bw() + 
  theme(legend.position = 'top')
dev.off()
```
![Volcano_plot](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/17295186-e173-4341-b152-a3e6f640bf14)
## Alternative to volcano plot with gene names
```
ggmaplot(filter_df1, fdr = 0.05, fc = 2, genenames= NULL, size = 1, alpha = 1, 
         font.label = c(9, "bold"), palette = c("#B31B21", "#1454AC", "darkgray"),
         top = 30, select.top.method = c("padj", "fc"), label.rectangle = TRUE, label.select = NULL, main = "N188.1 vs N222.1", xlab = "Log2 mean expression",
         ylab = "Log2 fold change", ggtheme = theme_classic(), legend = "top", font.main= "bold", font.legend = "bold",
)
dev.off()
```
![MA_plot](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/407dfa23-ec8a-424f-b7a8-f40e419c6fe6)
## Set up differentially expressed genes - Heatmap visualization
```
anno_df3 <- tibble::rownames_to_column(filter_df3, "Gene")
dim(anno_df3)
degs = anno_df3$Gene
vst_mat = assay(vst)
data_for_hm = vst_mat[degs,]
dim(data_for_hm)

pheatmap(data_for_hm, fontsize_row = 4, scale = 'row', cutree_cols = 2, cutree_rows = 2, border_color = NA)
dev.off()
```
![Heatmap](https://github.com/evangelosmourkas/Analyze-RNAseq-Data/assets/73548463/35750064-71f6-4c2f-86c2-c4b74ab3e152)
# How to cite
Monteith, W., Pascoe, B., Mourkas, E., Clark, J., Hakim, M., Hitchings, M.D., McCarthy, N., Yahara, K., Asakura, H., Sheppard, S.K. (2023), **Contrasting genes conferring short and long-term biofilm adaptation in _Listeria_** 
Microbial Genomics (in press). bioRxiv 2023.06.22.546149. https://doi.org/10.1101/2023.06.22.546149
