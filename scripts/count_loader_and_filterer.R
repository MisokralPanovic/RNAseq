.libPaths('C:/r_packages')
library(tidyverse)
library(data.table)

# load counts

dSH16 <- as.matrix(read.table('Data/dSH16vsMock.counts', header = T))
dSH40 <- as.matrix(read.table('Data/dSH40vsMock.counts', header = T))
WT16 <- as.matrix(read.table('Data/WT16vsMock.counts', header = T))
WT40 <- as.matrix(read.table('Data/WT40vsMock.counts', header = T))

# load count meta data

dSH16.meta <- as.data.frame(read.csv('Data/dSH16_meta.csv', header = T))
dSH16.meta <- dSH16.meta %>% remove_rownames %>% column_to_rownames(var="Names")

dSH40.meta <- as.data.frame(read.csv('Data/dSH40_meta.csv', header = T))
dSH40.meta <- dSH40.meta %>% remove_rownames %>% column_to_rownames(var="Names")

WT16.meta <- as.data.frame(read.csv('Data/WT16_meta.csv', header = T))
WT16.meta <- WT16.meta %>% remove_rownames %>% column_to_rownames(var="Names")

WT40.meta <- as.data.frame(read.csv('Data/WT40_meta.csv', header = T))
WT40.meta <- WT40.meta %>% remove_rownames %>% column_to_rownames(var="Names")


# chosing only rows that start with ENS

dSH16[rownames(dSH16) %like% "ENS", ]

nrow(dSH16[rownames(dSH16) %like% "ENS", ])
nrow(dSH16)

nrow(dSH40[rownames(dSH16) %like% "ENS", ])
nrow(dSH40)

nrow(WT16[rownames(dSH16) %like% "ENS", ])
nrow(WT16)

nrow(WT40[rownames(dSH16) %like% "ENS", ])
nrow(WT40)

