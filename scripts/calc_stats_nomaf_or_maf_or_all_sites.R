library(snpR)

# note, using nomaf from ANGSD because using a maf filter when most samples are not from 
# North America and are fixed at sites, which will
# result in a low maf at those sites and their removal!
dat <- readRDS("Raw_data/nomaf_snps_snpR.RDS")

#================with E and W NA split=====================
# pi, ho, pa
# dat <- calc_pi(dat, "pop")
# dat <- calc_ho(dat, "pop")
# dat <- calc_private(dat, "pop")
# 
# basic.stats <- get.snpR.stats(dat, "pop")
# pdat <- data.frame(pi = tapply(basic.stats$pi, basic.stats$subfacet, mean, na.rm = T),
#                    ho = tapply(basic.stats$ho, basic.stats$subfacet, mean, na.rm = T))
# pdat$pop <- rownames(pdat)
# 
# ## plot
# pdatm <- reshape2::melt(pdat, id.vars = c("pop"))
# ggplot(pdatm, aes(x = pop, y = value, color = variable)) + geom_point()
# basic.stats.m <- reshape2::melt(basic.stats[,c(2,6,7)], id.vars = c("subfacet"))
# colnames(basic.stats.m)[1] <- "pop"
# p <- ggplot(basic.stats.m, aes(x = pop, y = value, color = variable)) + geom_boxplot() + theme_bw()
# p + scale_y_continuous(limits = c(0,0.1))

#===============with E and W joined========================
# merge
datj <- as.data.frame(dat)
datj.samp.meta <- dat@sample.meta
datj.samp.meta[datj.samp.meta$pop %in% c("WNA", "ENA"),]$pop <- "NAM"
datj <- import.snpR.data(datj, dat@snp.meta, datj.samp.meta)

# pi, ho, private, fst
datj <- calc_pi(datj, "pop")
datj <- calc_ho(datj, "pop")
datj <- calc_private(datj, "pop")
datj <- calc_pairwise_fst(datj, "pop")

## plot
# basic.statsj <- get.snpR.stats(datj, "pop")
# pdatj <- data.frame(pi = tapply(basic.statsj$pi, basic.statsj$subfacet, mean, na.rm = T),
#                    ho = tapply(basic.statsj$ho, basic.statsj$subfacet, mean, na.rm = T))
# pdatj$pop <- rownames(pdat)
# pdatjm <- reshape2::melt(pdatj, id.vars = c("pop"))
# 
# bp1 <- ggplot(pdatjm, aes(x = pop, y = value, color = variable)) + geom_point() + theme_bw()
# basic.statsj$pa <- tapply(basic.statsj$pa, basic.statsj$subfacet, sum, na.rm = T)


# with a maf
datfj <- filter_snps(datj, maf = 0.05, maf.facets = "pop")
datfj <- calc_pairwise_fst(datfj, "pop")

# save
saveRDS(list(nomaf = datj, maf = datfj), "NAM_merged_noMAF_or_MAF_stats_calced_snpRdata.RDS")

rm(datj, datfj, dat)
gc(); gc(); gc(); gc();


#================all sites, joined==========================
dat <- readRDS("Raw_data/allsites_snps_snpR.RDS")

datj <- as.data.frame(dat)
datj.samp.meta <- dat@sample.meta
datj.samp.meta[datj.samp.meta$pop %in% c("WNA", "ENA"),]$pop <- "NAM"
datj <- import.snpR.data(datj, dat@snp.meta, datj.samp.meta)

# pi, ho, private, fst
datj <- calc_pi(datj, "pop")
datj <- calc_ho(datj, "pop")
datj <- calc_private(datj, "pop")
datj <- calc_tajimas_d(datj, "pop.group", sigma = 50, step = 50)

# save
saveRDS(datj, "NAM_merged_all_sites_stats_calced_snpRdata.RDS")

