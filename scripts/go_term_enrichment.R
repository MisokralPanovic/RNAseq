BiocManager::install("ViSEAGO")
library("ViSEAGO")

selection <- as.data.frame(read.csv('results/dSH40/sig_results.csv', header = T))
selection <- selection %>% dplyr::select(symbol) %>% filter(symbol != 'NA')

background <- as.data.frame(read.csv('results/dSH40/all_results.csv', header = T))
background <- background %>% dplyr::select(symbol) %>% filter(symbol != 'NA')


# load gene identifiers column 1) and corresponding statistical value (column 2)
table<-data.table::fread("table.txt")

# rank gene identifiers according statistical value
data.table::setorder(gene_table,padj)

# connect to Bioconductor
Bioconductor<-ViSEAGO::Bioconductor2GO()

# connect to EntrezGene
EntrezGene<-ViSEAGO::EntrezGene2GO()

# connect to Ensembl
Ensembl<-ViSEAGO::Ensembl2GO()

# connect to Uniprot-GOA
Uniprot<-ViSEAGO::Uniprot2GO()


# Display table of available organisms with Bioconductor
ViSEAGO::available_organisms(Bioconductor)

# Display table of available organisms with EntrezGene
ViSEAGO::available_organisms(EntrezGene)

# Display table of available organisms with Ensembl
ViSEAGO::available_organisms(Ensembl)

# Display table of available organisms with Uniprot
ViSEAGO::available_organisms(Uniprot)

# load GO annotations from Ensembl
myGENE2GO<-ViSEAGO::annotate(
  "btaurus_gene_ensembl",
  Ensembl
)

# perform fgseaMultilevel tests
BP<-ViSEAGO::runfgsea(
  geneSel=gene_table,
  ont="BP",
  gene2GO=myGENE2GO, 
  method ="fgseaMultilevel",
  params = list(
    scoreType = "pos",
    minSize=5
  )
)

