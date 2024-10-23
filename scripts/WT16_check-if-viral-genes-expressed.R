WT16 <- as.matrix(read.table('Data/WT16vsMock.counts', header = T))

# filters out viral transcripts (only leaves stuff starting with ENS)
nrow(WT16)
WT16 <- WT16[!rownames(WT16) %like% "ENS", ]
nrow(WT16)

WT16
