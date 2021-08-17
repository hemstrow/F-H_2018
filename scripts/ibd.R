library(snpR);
input <- readRDS("results/paralogs/nomaf_paralog_fix_snpR.RDS")
meta <- read.csv("Raw_data/sample_metadata.csv")
sample.meta(input) <- meta

#============marianas==============
mariana <- subset_snpR_data(input, pop = c("GUA", "ROT", "SAI"))
mariana <- filter_snps(mariana, 0.05, maf.facets = "pop")
mariana <- calc_isolation_by_distance(mariana, x_y = c("long", "lat"))
get.snpR.stats(mariana, stats = "ibd") 

#============north america==============
noram <- subset_snpR_data(input, pop = "NAM")
noram <- filter_snps(noram, 0.05)
noram <- calc_isolation_by_distance(noram, x_y = c("long", "lat"))
get.snpR.stats(noram, stats = "ibd")

#============hawaii==============
haw <- subset_snpR_data(input, pop = "HAW")
haw <- filter_snps(haw, 0.05)
haw <- calc_isolation_by_distance(haw, x_y = c("long", "lat"))
get.snpR.stats(haw, stats = "ibd")

#============australia==============
aus <- subset_snpR_data(input, pop = c("VIC", "QLD", "NSW"))
aus <- filter_snps(aus, 0.05)
aus <- calc_isolation_by_distance(aus, x_y = c("long", "lat"))
get.snpR.stats(aus, stats = "ibd")

#============overall==============
all <- filter_snps(input, 0.05, maf.facets = "pop")
all <- calc_isolation_by_distance(all, x_y = c("long", "lat"))
all <- calc_isolation_by_distance(all, "pop", x_y = c("long", "lat"))

get.snpR.stats(all, stats = "ibd")
