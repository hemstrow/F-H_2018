library(snpR);
input <- readRDS("results/paralogs/nomaf_paralog_fix_snpR.RDS")
meta <- read.csv("Raw_data/sample_metadata.csv")
sample.meta(input) <- meta

#============marianas==============
mariana <- subset_snpR_data(input, pop = c("GUA", "ROT", "SAI"))
mariana <- filter_snps(mariana, 0.05, maf_facets = "pop", min_loci = .75)
mariana <- calc_isolation_by_distance(mariana, x_y = c("long", "lat"))
get.snpR.stats(mariana, stats = "ibd") 

#============north america==============
noram <- subset_snpR_data(input, pop = "NAM")
noram <- filter_snps(noram, 0.05, min_loci = .75)
noram <- calc_isolation_by_distance(noram, x_y = c("long", "lat"))
get.snpR.stats(noram, stats = "ibd")

#============hawaii==============
haw <- subset_snpR_data(input, pop = "HAW")
haw <- filter_snps(haw, 0.05, min_loci = .75)
haw <- calc_isolation_by_distance(haw, x_y = c("long", "lat"))
get.snpR.stats(haw, stats = "ibd")

#============australia==============
aus <- subset_snpR_data(input, pop = c("VIC", "QLD", "NSW"))
aus <- filter_snps(aus, 0.05, min_loci = .75)
aus <- calc_isolation_by_distance(aus, x_y = c("long", "lat"))
get.snpR.stats(aus, stats = "ibd")



#==================filtered======================================================================
dat <- data.table::fread("data/dadi_inputs/rand_10kgap_snps.txt")
samp.met <- colnames(dat)[-c(1:2)]
samp.met <- data.frame(samp = samp.met, pop = substr(samp.met, 1, 3))
samp.met <- merge(samp.met, meta, by = "samp", sort = FALSE)
samp.met$pop.y <- NULL
colnames(samp.met)[2] <- "pop"

dat <- import.snpR.data(dat[,-c(1:2)], dat[,1:2], samp.met)
dat <- filter_snps(dat, hwe = 0.000001, hwe_facets = "pop")


#============marianas==============
mariana.f <- subset_snpR_data(dat, pop = c("GUA", "ROT", "SAI"))
mariana.f <- filter_snps(mariana.f, 0.05, maf_facets = "pop", min_loci = .75)
mariana.f <- calc_isolation_by_distance(mariana.f, x_y = c("long", "lat"))
get.snpR.stats(mariana.f, stats = "ibd") 

#============north america==============
noram.f <- subset_snpR_data(dat, pop = "NAM")
noram.f <- filter_snps(noram.f, 0.05, min_loci = .75)
noram.f <- calc_isolation_by_distance(noram.f, x_y = c("long", "lat"))
get.snpR.stats(noram.f, stats = "ibd")

#============hawaii==============
haw.f <- subset_snpR_data(dat, pop = "HAW")
haw.f <- filter_snps(haw.f, 0.05, min_loci = .75)
haw.f <- calc_isolation_by_distance(haw.f, x_y = c("long", "lat"))
get.snpR.stats(haw.f, stats = "ibd")

#============australia==============
aus.f <- subset_snpR_data(dat, pop = c("VIC", "QLD", "NSW"))
aus.f <- filter_snps(aus.f, 0.05, min_loci = .75)
aus.f <- calc_isolation_by_distance(aus.f, x_y = c("long", "lat"))
get.snpR.stats(aus.f, stats = "ibd")
