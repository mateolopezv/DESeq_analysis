#Differential Gene Expression Analysis

#Installing packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")

#Packages

library(DESeq2)
library(tidyverse)
library(airway)

#Data

data(airway)
airway

#Load csv
counts_data <- read.csv("C:/Users/ignac/Desktop/Airway/counts_data.csv")
colData <- read.csv("C:/Users/ignac/Desktop/Airway/sample_info.csv")

#Naming columns in colData
names(colData) <- c("NA", "cellLine", "dexamethasone")
colData <- select(colData, cellLine, dexamethasone)

#Columns in counts data match rows in colData
all(colnames(counts_data) %in% rownames(colData))

#Same order?
all(colnames(counts_data) == rownames(colData))

#DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~dexamethasone)
View(dds)

#Filtering samples with low counts (>= 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#Set factors levels
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

#Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

#Explore results
summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

#Contrast
resultsNames(dds)

#This only has two levels: treated and untreated -> If you have +2 you can use contrast
results(dds, contrast = c("dexamethasone", "treated 4hs", "untreated")) #e.g.

#MA plot
plotMA(res)

#Volcano Plot
#Colors
old.pal <- palette(c("#00BFFF", "#FF3030"))
#Plot
plot(res$log2FoldChange, -log10(res$padj), main = "Treated vs Untreated", xlab = "Log2(FC)", ylab="-log10(padj)")
with(subset(res, padj < 0.05 & abs(log2FoldChange) >= 1),
     points(log2FoldChange, -log10(padj), col=(sign(log2FoldChange)+6)/2))
legend("topright", title = paste("padj < 0.05"), legend = c("Up", "Down"), col = 1:2, pch = 20)

#Heatmap
normalized_counts <- counts(dds, normalized = TRUE)
top_hits <- res[order(res$padj),][1:10,]
top_hits <- row.names(top_hits)
cal_z_score <- function(x) {(x - mean(x))/sd(x)}
zscore_all <- t(apply(normalized_counts,1,cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset)
