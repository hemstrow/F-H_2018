library(snpR);
input <- readRDS("results/paralogs/nomaf_paralog_fix_snpR.RDS")
input@sample.meta <- read.csv("Raw_data/sample_metadata.csv")
dats <- import.snpR.data(input[1:nrow(input), 1:ncol(input)], input@snp.meta, input@sample.meta)
rm(input)

#============marianas==============
mariana <- subset_snpR_data(dats, facets = "pop", subfacets = c("GUA", "ROT", "SAI"))
mariana <- filter_snps(mariana, 0.05, maf.facets = "pop")
mariana <- calc_isolation_by_distance(mariana, x_y = c("long", "lat"))
get.snpR.stats(mariana, type = "ibd")

#============north america==============
noram <- subset_snpR_data(dats, facets = "pop", subfacets = "NAM")
noram <- filter_snps(noram, 0.05)
noram <- calc_isolation_by_distance(noram, x_y = c("long", "lat"))
get.snpR.stats(noram, type = "ibd")

#============hawaii==============
haw <- subset_snpR_data(dats, facets = "pop", subfacets = "HAW")
haw <- filter_snps(haw, 0.05)
haw <- calc_isolation_by_distance(haw, x_y = c("long", "lat"))
get.snpR.stats(haw, type = "ibd")

#============australia==============
aus <- subset_snpR_data(dats, facets = "pop", subfacets = c("VIC", "QLD", "NSW"))
aus <- filter_snps(aus, 0.05)
aus <- calc_isolation_by_distance(aus, x_y = c("long", "lat"))
get.snpR.stats(aus, type = "ibd")

#============overall==============
all <- filter_snps(dats, 0.05, maf.facets = "pop")
all <- calc_isolation_by_distance(all, x_y = c("long", "lat"))
all <- calc_isolation_by_distance(all, "pop", x_y = c("long", "lat"))

get.snpR.stats(all, type = "ibd")
