library(dplyr)
x <- readr::read_delim("~/monarch/github/F-H_2018/Raw_data/nomaf_raw_genotypes.geno", delim = "\t", col_names = F)
colnames(x)[1:2] <- c("chr", "position")

#function to get n random snps every gap bases
get_rand <- function(x, gap, n){
  cs <- seq(gap/2, max(x$position) + gap/2, gap)
  starts <- cs - gap/2
  ends <- cs + gap/2
  pos <- x$position
  
  lmat <- outer(pos, starts, function(pos, starts) pos >= starts)
  lmat <- lmat + outer(pos, ends, function(pos, ends) pos < ends)
  colnames(lmat) <- cs
  rownames(lmat) <- pos
  lmat <- ifelse(lmat == 2, TRUE, FALSE)
  
  wins <- which(t(lmat) == TRUE) %% length(cs)
  wins[wins == 0] <- length(cs)
  
  x <- cbind(window = wins, x)
  
  rands <- x %>% group_by(window) %>% sample_n(size = n)
  
  rands <- rands[,-1]
  
  return(rands)
}

#run this
rsnps <- plyr::ddply(x, "chr", .fun = get_rand, gap = 20000, n = 1, .progress = "text")

#save
data.table::fwrite(rsnps, "~/monarch/github/F-H_2018/Data/rand_gap_snps.txt", quote = FALSE, col.names = T, sep = "\t", row.names = F)
