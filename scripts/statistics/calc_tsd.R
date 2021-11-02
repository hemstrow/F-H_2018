.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/4.1", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR);

args <- commandArgs(TRUE)
input <- args[1]
pop <- args[2]
outfile <- gsub(".RDS", "", input)
dat <- readRDS(input)

dat <- dat[pop = pop]
dat <- filter_snps(dat, min_loci = 0.75, non_poly = F, re_run = FALSE, hwe = 0.000001, hwe_facets = "pop")
if(nrow(dat) > 0){
    dat <- calc_tajimas_d(dat, "group", sigma = 50, step = 50)
    write.table(get.snpR.stats(dat, "group", "tajimas_d")$weighted.means, paste0(outfile, "_", pop, "_tajimas_D.txt"), col.names = T, quote = F, sep = "\t")
    saveRDS(dat, paste0(outfile, "_", pop, "_processed_tajimasD.RDS"))
}
