# package loader ----
#.libPaths('C:/r_packages')

library(DESeq2)
# library(edgeR)
library(org.Bt.eg.db)
library("biomaRt")
library("AnnotationDbi")
library(tximport)
library(vsn)
library(SummarizedExperiment)

library(here)
library(LSD)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)
library(data.table)
source(here("scripts/volcano_plot.R"))
pal <- brewer.pal(11,"RdBu")

# count data and meta data loader ---------------------

dSH16 <- as.matrix(read.table('Data/dSH16vsMock.counts', header = T))

dSH16.meta <- as.data.frame(read.csv('Data/dSH16_meta.csv', header = T))

# reorder the conditions so they go control, condition
dSH16 <- dSH16[,sort(colnames(dSH16))]

dSH16.meta <- dSH16.meta %>% arrange(Names)

# removes the colum name of sample names
dSH16.meta <- dSH16.meta %>% remove_rownames %>% column_to_rownames(var="Names")

# filters out viral transcripts (only leaves stuff starting with ENS)
nrow(dSH16)
dSH16 <- dSH16[rownames(dSH16) %like% "ENS", ]
nrow(dSH16)

# check if count sample names and metdata sample names are the same
all(rownames(dSH16.meta) %in% colnames(dSH16))

write.csv(dSH16, file = "Data/dSH16.csv", row.names = TRUE)

# DESeq2 object creation ----------------------
dds <- DESeqDataSetFromMatrix(countData = dSH16,
                              colData = dSH16.meta,
                              design = ~ Condition)

dds$Condition <- relevel(dds$Condition, ref = "untreated")

# filering genes with no counts
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

# data transformation ----------
meanSdPlot(log(assay(dds)[rowSums(assay(dds))>30,]))

ggsave('MeanSdPlot-dds.png',
       path = 'results/dSH16')

# rlog transformation ----------
rld <- rlog(dds)
head(assay(rld))
meanSdPlot(log(assay(rld)[rowSums(assay(rld))>30,]))

ggsave('MeanSdPlot-rld.png',
       path = 'results/dSH16')

# variance stabilising transformation (VST) -----
vsd <- DESeq2::vst(dds)
meanSdPlot(log(assay(vsd)[rowSums(assay(vsd))>0,]))

ggsave('MeanSdPlot-vsd.png',
       path = 'results/dSH16')

# pca plot -------------
DESeq2::plotPCA(rld, intgroup = 'Condition')

ggsave('PCA-rld.png',
       path = 'results/dSH16')

DESeq2::plotPCA(vsd, intgroup = 'Condition')

ggsave('PCA-vsd.png',
       path = 'results/dSH16')

# sample distances vsd ------------
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# save it manually 'distance-matrix-vsd.png'

# sample distances rld -------------
sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# save it manually 'distance-matrix-rld.png'

#### Performing differential expression testing with DESeq2 ---------------
dds <- DESeq2::DESeq(dds)

# plot the dispersion estimation -----------
DESeq2::plotDispEsts(dds)

ggsave('dispersion-estimation.png',
       plot = DESeq2::plotDispEsts(dds),
       path = 'results/dSH16')

# DE results ---------
res <- results(dds)
head(res)
mcols(res, use.names = TRUE)
summary(res)
res <- res[order(res$padj),]
head(res)

hist(res$pvalue,breaks=seq(0,1,.01))

ggsave('p_value-distribution.png',
       plot = hist(res$pvalue,breaks=seq(0,1,.01)),
       path = 'results/dSH16')

# Following Schurch et al., (RNA, 2016) recommandations,  for 3 replicates per conditions:
resSchurch <- results(dds, lfcThreshold = 0.5, alpha = 0.01)
summary(resSchurch)

res.vgp <- results(dds, lfcThreshold = 1, alpha = 0.05)
summary(res.vgp)

# volcano plot --------
volcanoPlot(res.vgp)
ggsave('volcanoPlot-res.vgp.png',
       plot = volcanoPlot(res.vgp),
       path = 'results/dSH16')


volcanoPlot(resSchurch)
ggsave('volcanoPlot-resSchurch.png',
       volcanoPlot(resSchurch),
       path = 'results/dSH16')

# MA plot ----------
DESeq2::plotMA(res.vgp, ylim = c(-5, 5))
ggsave('MA-plot_res.vgp.png',
       DESeq2::plotMA(res.vgp, ylim = c(-5, 5)),
       path = 'results/dSH16')

DESeq2::plotMA(resSchurch, ylim = c(-5, 5))
ggsave('MA-plot_resSchurch.png',
       DESeq2::plotMA(resSchurch, ylim = c(-5, 5)),
       path = 'results/dSH16')

# count plots
topGene <- rownames(res.vgp)[which.min(res.vgp$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Condition"))

ggsave('top-gene-count.png',
       plotCounts(dds, gene = topGene, intgroup=c("Condition")),
       path = 'results/dSH16')

# adding annotation -----------
columns(org.Bt.eg.db)

# for normal results
## substr: data, start, stop
ens.str <- substr(rownames(res),1,18)
res.vgp$symbol <- mapIds(org.Bt.eg.db,
                         keys=ens.str,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
res.vgp$gene_name <- mapIds(org.Bt.eg.db,
                            keys=ens.str,
                            column="GENENAME",
                            keytype="ENSEMBL",
                            multiVals="first")
res.vgp$go <- mapIds(org.Bt.eg.db,
                     keys=ens.str,
                     column="GO",
                     keytype="ENSEMBL",
                     multiVals="first")
res.vgp$ontology <- mapIds(org.Bt.eg.db,
                           keys=ens.str,
                           column="ONTOLOGY",
                           keytype="ENSEMBL",
                           multiVals="first")

# resSchurch
resSchurch$symbol <- mapIds(org.Bt.eg.db,
                            keys=ens.str,
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
resSchurch$gene_name <- mapIds(org.Bt.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
resSchurch$go <- mapIds(org.Bt.eg.db,
                        keys=ens.str,
                        column="GO",
                        keytype="ENSEMBL",
                        multiVals="first")
resSchurch$ontology <- mapIds(org.Bt.eg.db,
                              keys=ens.str,
                              column="ONTOLOGY",
                              keytype="ENSEMBL",
                              multiVals="first")

# subset significant results
resSig <- subset(res.vgp, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

# subset significant results for resSchurch
resSigSchurch <- subset(resSchurch, padj < 0.01)
head(resSigSchurch[ order(resSigSchurch$log2FoldChange), ])
head(resSigSchurch[ order(resSigSchurch$log2FoldChange, decreasing = TRUE), ])

# heatmaps --------
# plots difference of mean row from vsd dataset of top genes
mat <- assay(vsd)[head(order(res.vgp$padj), 30), ]
mat <- mat - rowMeans(mat)
pheatmap(mat)

# save manually 'heatmap-top-30-sig-genes_vsd-mean-row-difference.png'

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)

matros  <- assay(vsd)[ topVarGenes, ]
matros  <- matros - rowMeans(matros)

pheatmap(matros)

# save manually 'heatmap-top-variable-genes_vsd-mean-row-difference.png',

# filter results based on p value
resOrdered <- res.vgp[order(res.vgp$pvalue),]

# saving it as .csv dataframe
resOrderedDF_top100 <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF_top100, file = "results/dSH16/results_top100.csv")

# only significant results
write.csv(resSig[order(resSig$pvalue),], file = "results/dSH16/sig_results.csv")

# resSchurch results
write.csv(resSigSchurch[order(resSigSchurch$pvalue),], file = "results/dSH16/Schurch_results.csv")

# save all genes
write.csv(res.vgp[order(res.vgp$pvalue),], file = "results/dSH16/all_results.csv")

# makinghtml report
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="results/dSH16/report")
publish(resSig, htmlRep)
url <- finish(htmlRep)
browseURL(url)
