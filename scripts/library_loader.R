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


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("vsn")

source(here("src/R/volcanoPlot.R"))
