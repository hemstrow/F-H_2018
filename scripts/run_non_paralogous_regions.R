.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(dplyr); library(readr); source("scripts/convert_paralog_output_to_selected_regions.R")

#args <- commandArgs(TRUE)
#ref <- args[1]
ref <- "Raw_data/chr_lengths.txt"
result_dir <- "results/paralogs/"
result_pattern <- ".paralogs"
output <- "results/paralogs/selected_clean_regions.txt"

g.list <- find_non_paralogous_sections(result_dir, result_pattern, ref,
                                       1000)

# check
chrs <- gsub(":.+$", "", g.list)
bases <- gsub("^.+:", "", g.list)
sbase <- gsub("-", "", gsub("-.+", "", bases))
fbase <- gsub("^.+-", "", bases)
fbase[fbase == ""] <- Inf
fbase <- as.numeric(fbase)
sbase <- as.numeric(sbase)
any(sbase > fbase)

write.table(data.frame(p = g.list), file = output, col.names = F, row.names = F, quote = F)

