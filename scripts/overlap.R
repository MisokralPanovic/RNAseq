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

overlap_vector <- c()
for (gene in sig_WT40$symbol) {
  if (gene %in% sig_dSH40$symbol) {
    print(gene)
    vector_out <- c(T)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
  else {
    vector_out <- c(F)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
}
sig_WT40 <- sig_WT40 %>% mutate(Overlap = overlap_vector)

overlap_vector <- c()
for (gene in sig_dSH40$symbol) {
  if (gene %in% sig_WT40$symbol) {
    print(gene)
    vector_out <- c(T)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
  else {
    vector_out <- c(F)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
}
sig_dSH40 <- sig_dSH40 %>% mutate(Overlap = overlap_vector)

overlap_vector <- c()
for (gene in Schurch_WT40$symbol) {
  if (gene %in% Schurch_dSH40$symbol) {
    print(gene)
    vector_out <- c(T)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
  else {
    vector_out <- c(F)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
}
Schurch_WT40 <- Schurch_WT40 %>% mutate(Overlap = overlap_vector)

overlap_vector <- c()
for (gene in Schurch_dSH40$symbol) {
  if (gene %in% Schurch_WT40$symbol) {
    print(gene)
    vector_out <- c(T)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
  else {
    vector_out <- c(F)
    overlap_vector <- c(overlap_vector,
                        vector_out)
  }
}
Schurch_dSH40 <- Schurch_dSH40 %>% mutate(Overlap = overlap_vector)

write.csv(sig_WT40, file = "results/WT40/overlap_sig.csv")
write.csv(sig_dSH40, file = "results/dSH40/overlap_Schurch.csv")
write.csv(Schurch_WT40, file = "results/WT40/overlap_sig.csv")
write.csv(Schurch_dSH40, file = "results/dSH40/overlap_Schurch.csv")
