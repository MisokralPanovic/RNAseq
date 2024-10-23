.libPaths('C:/r_packages')
library(tidyverse)
library(data.table)

sig_WT40 <- as.data.frame(read.csv('results/WT40/sig_results.csv', header = T))
Schurch_WT40 <- as.data.frame(read.csv('results/WT40/Schurch_results.csv', header = T))

sig_dSH40 <- as.data.frame(read.csv('results/dSH40/sig_results.csv', header = T))
Schurch_dSH40 <- as.data.frame(read.csv('results/dSH40/Schurch_results.csv', header = T))

sig_WT16 <- as.data.frame(read.csv('results/WT16/sig_results.csv', header = T))
Schurch_WT16 <- as.data.frame(read.csv('results/WT16/Schurch_results.csv', header = T))

sig_dSH16 <- as.data.frame(read.csv('results/dSH16/sig_results.csv', header = T))
Schurch_dSH16 <- as.data.frame(read.csv('results/dSH16/Schurch_results.csv', header = T))


write.csv(sig_WT40$X, file = "results/WT40/ensembl_sig.csv")
write.csv(Schurch_WT40$X, file = "results/WT40/ensembl_list_Schurch.csv")
write.csv(sig_WT40$symbol, file = "results/WT40/gene_list_sig.csv")
write.csv(Schurch_WT40$symbol, file = "results/WT40/gene_list_Schurch.csv")

write.csv(sig_dSH40$X, file = "results/dSH40/ensembl_sig.csv")
write.csv(Schurch_dSH40$X, file = "results/dSH40/ensembl_Schurch.csv")
write.csv(sig_dSH40$symbol, file = "results/dSH40/gene_list_sig.csv")
write.csv(Schurch_dSH40$symbol, file = "results/dSH40/gene_list_Schurch.csv")

write.csv(sig_WT16$X, file = "results/WT16/ensembl_sig.csv")
write.csv(Schurch_WT16$X, file = "results/WT16/ensembl_Schurch.csv")
write.csv(sig_WT16$symbol, file = "results/WT16/gene_list_sig.csv")
write.csv(Schurch_WT16$symbol, file = "results/WT16/gene_list_Schurch.csv")

write.csv(sig_dSH16$X, file = "results/dSH16/ensembl_sig.csv")
write.csv(Schurch_dSH16$X, file = "results/dSH16/ensembl_Schurch.csv")
write.csv(sig_dSH16$symbol, file = "results/dSH16/gene_list_sig.csv")
write.csv(Schurch_dSH16$symbol, file = "results/dSH16/gene_list_Schurch.csv")
