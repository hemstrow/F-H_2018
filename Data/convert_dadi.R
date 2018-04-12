library(snpR)
library(readr)

#import data
anc <- read.table("anc_flank.txt", header = F)
anc <- anc[-1,]
ref <- read.table("ref_flank.txt", header = F)
ref <- ref[-1,]
dat <- readr::read_delim("rand_gap_snps.txt", delim = "\t", col.names = T)

#bind together
x <- cbind(ref = ref, anc = anc, dat)
