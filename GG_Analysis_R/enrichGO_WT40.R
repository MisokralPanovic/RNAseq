if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")
BiocManager::install("clusterProfiler")
BiocManager::install("biomaRt")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("Biobase")
BiocManager::install("DO.bd")
#library(topGO)

BiocManager::install("enrichplot")
library(enrichplot)
install.packages("XML")

install.packages("ggplot2")
library(ggplot2)

install.packages("UpSetR")
install.packages("ggupset")
library(UpSetR)
library(ggupset)

mart=useMart('ensembl')
listDatasets(mart)
####

db= useMart('ENSEMBL_MART_ENSEMBL',dataset='btaurus_gene_ensembl', host="useast.ensembl.org", ensemblRedirect = FALSE)
WT40<-read.table(file = "all_results_WT40.csv", header = T, sep = ",", dec = ".")
WT40 <- dplyr::rename(WT40,ensembl_gene_id=X)
#WT40<- dplyr::filter(WT40,pvalue<0.05)
WT40<-dplyr::filter(WT40, log2FoldChange>1)
GeneID_WT40= getBM(attributes=c('entrezgene_id', 'ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=WT40$ensembl_gene_id, mart=db)


GeneID_WT40 <- dplyr::filter(GeneID_WT40,!is.na(entrezgene_id))
GeneID_WT40 <- dplyr::left_join(GeneID_WT40,WT40,by="ensembl_gene_id")

Input_WT40 <- GeneID_WT40$log2FoldChange
names(Input_WT40) <- GeneID_WT40$entrezgene_id
head(Input_WT40)



data(GeneID_WT40, package="DOSE")
deWT40 <- names(Input_WT40)[abs(Input_WT40) > 1] #Seuil FC fixed a 2, peut être modifié
ego_BP_WT40 <- enrichGO(gene = deWT40, OrgDb = "org.Hs.eg.db", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable=TRUE)
#ego_MF <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db", ont="MF",pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable=TRUE)
ego_CC2vs2 <- enrichGO(gene = de2vs2, OrgDb = "org.Hs.eg.db", ont="CC",pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable=TRUE)

head(ego_BP2vs2)
as.data.frame(ego_BP2vs2)

goplot(ego_BP2vs2)
#goplot(ego_MF)#, Aucun enrichissement significatif en MF
goplot(ego_CC2vs2)


plotGOgraph(ego_BP2vs2)
#plotGOgraph(ego_MF)
plotGOgraph(ego_CC2vs2)

barplot(ego_BP2vs2, showCategory=20)
dotplot(ego_BP2vs2, showCategory=30)


go2vs2 <- enrichGO(de2vs2, OrgDb = "org.Hs.eg.db", ont="all")

dotplot(go2vs2, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
## remove redundent GO terms
ego2_BP2vs2 <- simplify(ego_BP2vs2)
cnetplot(ego2_BP2vs2, foldChange=Input_geneList2vs2)

ego2_MF <- simplify(ego_MF)
#cnetplot(ego2_MF, foldChange=Input_geneList)

ego2_CC2vs2 <- simplify(ego_CC2vs2)
cnetplot(ego2_CC2vs2, foldChange=Input_geneList2vs2)

cnetplot(ego2_BP2vs2, foldChange=Input_geneList2vs2, circular = TRUE, colorEdge = TRUE)

#cnetplot(ego2_MF, foldChange=Input_geneList, circular = TRUE, colorEdge = TRUE)

cnetplot(ego2_CC2vs2, foldChange=Input_geneList2vs2, circular = TRUE, colorEdge = TRUE)

upsetplot(ego_BP2vs2)
#upsetplot(ego_MF)
upsetplot(ego_CC)
heatplot(ego2_BP2vs2)
#heatplot(ego2_MF)
heatplot(ego2_CC2vs2)

heatplot(ego2_BP2vs2, foldChange=Input_geneList2vs2)
#heatplot(ego2_MF, foldChange=Input_geneList)
heatplot(ego2_CC2vs2, foldChange=Input_geneList2vs2)


emapplot(ego2_BP2vs2)
#emapplot(ego2_MF)
emapplot(ego2_CC2vs2)
