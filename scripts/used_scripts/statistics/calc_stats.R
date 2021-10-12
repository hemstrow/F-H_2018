.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR);

args <- commandArgs(TRUE)

if(length(args) != 3){
  cat("Usage: calc_stats.R input_file FST_TRUE_or_FALSE Tajima'sD_TRUE_or_FALSE Tajima'sD_pop\n")
}
if(length(args) == 3){
  input <- args[1]
  fst <- args[2]
  tajimas <- args[3]

  outfile <- gsub(".RDS", "", input)

  dat <- readRDS(input)
  dat <- filter_snps(dat, min_ind = 0.5, non_poly = F)

  # pi, ho, private, fst, tajima's D
  dat <- calc_pi(dat, "pop")
  dat <- calc_ho(dat, "pop")
  dat <- calc_private(dat, "pop")
  write.table(get.snpR.stats(dat, "pop"), paste0(outfile, "_basic_stats_ind_filtered.txt"), col.names = T, quote = F, sep = "\t")

  # het/hom ratio
  dat <- calc_het_hom_ratio(dat)
  write.table(dat@sample.stats, paste0(outfile, "_hom_het_ratio_ind_filtered.txt"), col.names = T, quote = F, sep = "\t")

  if(fst == TRUE){
    dat <- calc_pairwise_fst(dat, "pop")
    write.table(get.snpR.stats(dat, "pop", "pairwise"), paste0(outfile, "_fst_ind_filtered.txt"), col.names = T, quote = F, sep = "\t")
  }
  #if(tajimas == TRUE){
  #  dat <- calc_tajimas_d(dat, "pop.group", sigma = 50, step = 50)
  #  write.table(get.snpR.stats(dat, "pop", "single.window"), paste0(outfile, "tajimas_D.txt"), col.names = T, quote = F, sep = "\t")
  #}

  saveRDS(dat, paste0(outfile, "_ind_filtered_processed.RDS"))
}
