.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR);

args <- commandArgs(TRUE)
input <- args[1]
pop <- args[2]
outfile <- gsub(".RDS", "", input)
dat <- readRDS(input)

dat <- subset_snpR_data(dat, facets = "pop", subfacets = pop)
dat <- filter_snps(dat, min_loci = 0.75, non_poly = F, re_run = FALSE)
dat <- calc_tajimas_d(dat, "group", sigma = 50, step = 50)
write.table(get.snpR.stats(dat, "group", "single.window"), paste0(outfile, "_", pop, "_tajimas_D.txt"), col.names = T, quote = F, sep = "\t")
saveRDS(dat, paste0(outfile, "_", pop "_processed_tajimasD.RDS"))
