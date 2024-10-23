results_table <- as.data.frame(read.csv('results/WT16/sig_results.csv', header = T))
results_table <- results_table %>% dplyr::select(symbol)

genes <- c('SLC16A1',
  
  'SLC1A1',
  
  'SLC25A28',
  
  'SLC25A30',
  
  'SLC26A4',
  
  'SLC2A12',
  
  'SLC6A9',
  
  'SLCO5A1',
  
  'SLC15A3')


for (gene in genes) {
 if (gene %in% results_table$symbol) {
   print(gene)
 }
}


dSH40_slc <- c( "SLC16A1",
                "SLC1A1",
               "SLC6A9",
               "SLC15A3")
dSH16_scl <- c("SLC16A1")
WT40_scl <- c("SLC16A1",
              "SLC1A1",
              "SLC6A9",
              "SLC15A3")
WT16_scl <- c()