.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR);

args <- commandArgs(TRUE)
pop <- as.character(args[1]) #one pop at a time, since the file is too large otherwise!

#================all sites, joined==========================
dat <- readRDS("../Raw_data/allsites_snps_snpR.RDS")

# prep
datj <- as.data.frame(dat)
datj.samp.meta <- dat@sample.meta
datj.samp.meta[datj.samp.meta$pop %in% c("WNA", "ENA"),]$pop <- "NAM"
datj <- import.snpR.data(datj, dat@snp.meta, datj.samp.meta)
datj <- subset_snpR_data(datj, facets = "pop", subfacets = pop)

# pi, ho, private, fst
datj <- calc_pi(datj, "pop")
datj <- calc_ho(datj, "pop")
datj <- calc_private(datj, "pop")
datj <- calc_tajimas_d(datj, "pop.group", sigma = 50, step = 50)


# write stats
write.table(get.snpR.stats(datj, "pop"), paste0("../results/all_sites_single_stats_", pop, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(get.snpR.stats(datj, "pop.group", "single.window"), paste0("../results/all_sites_tajimasD_", pop, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t")

# save
saveRDS(datj, paste0("../results/NAM_merged_all_sites_stats_calced_snpRdata_", pop, ".RDS"))

