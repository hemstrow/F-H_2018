.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR);

#================all sites, joined==========================
dat <- readRDS("../Raw_data/allsites_snps_snpR.RDS")

# prep
datj <- as.data.frame(dat)
datj.samp.meta <- dat@sample.meta
datj.samp.meta[datj.samp.meta$pop %in% c("WNA", "ENA"),]$pop <- "NAM"
datj <- import.snpR.data(datj, dat@snp.meta, datj.samp.meta)

# pi, ho, private, fst
datj <- calc_pi(datj, "pop")
datj <- calc_ho(datj, "pop")
datj <- calc_private(datj, "pop")
datj <- calc_tajimas_d(datj, "pop.group", sigma = 50, step = 50)


# write stats
write.table(get.snpR.stats(datj, "pop"), "../results/all_sites_single_stats.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(get.snpR.stats(datj, "pop.group", "single.window"), "../results/all_sites_tajimasD.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# save
saveRDS(datj, "../results/NAM_merged_all_sites_stats_calced_snpRdata.RDS")

