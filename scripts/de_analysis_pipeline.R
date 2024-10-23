.libPaths('C:/r_packages')

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(here)
  library(LSD)
  library(RColorBrewer)
  library(SummarizedExperiment)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

source(here("scripts/volcano_plot.R"))
pal <- brewer.pal(8,"Dark2")

#---------------------


WT40 <- as.matrix(read.table('Data/WT40vsMock.counts', header = T))

WT40.meta <- as.data.frame(read.csv('Data/WT40_meta.csv', header = T))

WT40.meta <- WT40.meta %>% remove_rownames %>% column_to_rownames(var="Names")


all(rownames(WT40.meta) %in% colnames(WT40))

wt40.dds <- DESeqDataSetFromMatrix(countData = WT40,
                              colData = WT40.meta,
                              design = ~ Condition)



meanSdPlot(log(assay(wt40.dds)[rowSums(assay(wt40.dds))>30,]))

wt40.dds <- DESeq2::DESeq(wt40.dds)
DESeq2::plotDispEsts(wt40.dds)

resultsNames(wt40.dds)
res.wt40 <- results(wt40.dds,name="Condition_untreated_vs_treated")

head(res.wt40)
summary(res.wt40)

hist(res.wt40$pvalue,breaks=seq(0,1,.01))


resSchurch.wt40 <- results(wt40.dds, name="Condition_untreated_vs_treated", lfcThreshold = 0.5, alpha = 0.01)
summary(resSchurch.wt40)


volcanoPlot(resSchurch.wt40)

DESeq2::plotMA(res.wt40, ylim = c(-5, 5))


mat <- assay(vsd)[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("Tissue", "Time")])
pheatmap(mat, annotation_col = df)

resOrdered.wt40 <- res.wt40[order(res.wt40$pvalue),]
head(resOrdered.wt40, n = 10L)

resOrdered.wt40Schurch <- resSchurch.wt40[order(resSchurch.wt40$pvalue),]
head(resOrdered.wt40Schurch,  n = 20L)



resSig.wt40 <- subset(resOrdered, padj < 0.1)





