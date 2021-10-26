.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/4.1", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR);

args <- commandArgs(TRUE)

if(length(args) != 1){
  cat("Usage: calc_stats.R input_snpRdata_file\n")
}
if(length(args) == 1){
  input <- args[1]

  outfile <- gsub(".RDS", "", input)

  dat <- readRDS(input)
  dat <- filter_snps(dat, hwe = 0.000001, hwe_facets = "pop", min_loci = 0.75)

  # pi, ho, private, fst, tajima's D
  dat <- calc_pi(dat, "pop")
  dat <- calc_ho(dat, "pop")
  dat <- calc_private(dat, "pop")
  dat <- calc_het_hom_ratio(dat)

  ss <- get.snpR.stats(dat, "pop", c("pi", "ho"))$weighted.means
  ss <- na.omit(ss)
  write.table(ss, paste0(outfile, "_processed_stats.txt"), sep = "\t", quote = F, col.names = TRUE, row.names = FALSE)
  
  saveRDS(dat, paste0(outfile, "_processed.RDS"))
}
