.libPaths('C:/r_packages')
library(tidyverse)
library(data.table)


count(dsh40_overlap_Schurch %>% filter(Overlap == T))

dsh40_overlap_Schurch <- dsh40_overlap_Schurch %>% filter(symbol != 'NA')
wt40_overlap_Schurch <- wt40_overlap_Schurch %>% filter(symbol != 'NA')


library(ggvenn)

ggvenn(
  list(WT40 = wt40_overlap_Schurch$symbol,
       dSH40 = dsh40_overlap_Schurch$symbol),
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)



list(WT40 = wt40_overlap_Schurch$symbol,
     dSH40 = dsh40_overlap_Schurch$symbol)

intersect(wt40_overlap_Schurch$symbol,
          dsh40_overlap_Schurch$symbol)
