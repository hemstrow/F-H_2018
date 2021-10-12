library(snpR); library(data.table)

anc_flank <- "~/monarch/github/F-H_2018/data/dadi_inputs/anc_flank.txt"
ref_flank <- "~/monarch/github/F-H_2018/data/dadi_inputs/anc_flank.txt"
idat <- "~/monarch/github/F-H_2018/data/dadi_inputs/rand_10kgap_snps.txt"
odat <- "~/monarch/github/F-H_2018/data/dadi_inputs/dadi_10kgap_snps.txt"

#import data
anc <- read.table(anc_flank, header = F, stringsAsFactors = F)
anc <- anc[-1,]
ref <- read.table(ref_flank, header = F, stringsAsFactors = F)
ref <- ref[-1,]
dat <- read.table(idat, stringsAsFactors = F, header = T)

#bind together
x <- cbind(ref = ref, anc = anc, dat, stringsAsFactors = F)

######################

#sort and reformat
#if there is a wierd NA column at the end, remove it.
if(all(is.na(x[,ncol(x)]))){
    x <- x[,-ncol(x)]
}

pops <- table(substr(colnames(x[,5:ncol(x)]), 1, 3))
pops <- list(names(pops), as.numeric(pops))

#reformat
dadi_snps <- format_snps(x, 4, "dadi", pop = pops)

data.table::fwrite(dadi_snps, odat, quote = F, sep = "\t", col.names = T, row.names = F)
