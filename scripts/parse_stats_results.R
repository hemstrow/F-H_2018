stats.nomaf <- read.table("results/nomaf_single_stats.txt", header = T, stringsAsFactors = F)
pdatj <- data.frame(pi = tapply(stats.nomaf$pi, stats.nomaf$subfacet, mean, na.rm = T),
                    ho = tapply(stats.nomaf$ho, stats.nomaf$subfacet, mean, na.rm = T))
tapply(stats.nomaf$pa, stats.nomaf$subfacet, sum, na.rm = T)
pdatj$pop <- rownames(pdatj)
pdatj <- reshape2::melt(pdatj, id.vars = "pop")
library(ggplot2)
ggplot(pdatj, aes(x = pop, y = value, color = variable)) + geom_point() + theme_bw()




stats.maf <- read.table("results/maf_single_stats.txt", header = T, stringsAsFactors = F)
pdatj <- data.frame(pi = tapply(stats.maf$pi, stats.maf$subfacet, mean, na.rm = T),
                    ho = tapply(stats.maf$ho, stats.maf$subfacet, mean, na.rm = T))
tapply(stats.maf$pa, stats.maf$subfacet, sum, na.rm = T)
pdatj$pop <- rownames(pdatj)
pdatj <- reshape2::melt(pdatj, id.vars = "pop")
library(ggplot2)
ggplot(pdatj, aes(x = pop, y = value, color = variable)) + geom_point() + theme_bw()



pairwise.maf <- readr::read_delim("results/maf_pairwise_stats.txt", col_names = T, delim = "\t")
pdat <- data.frame(fst = tapply(pairwise.maf$fst, pairwise.maf$comparison, mean, na.rm = T))
pdat$comparison <- rownames(pdat)
pdat$pop_1 <- substr(pdat$comparison, 1, 3)
pdat$pop_2 <- substr(pdat$comparison, 5, 7)
all_pops <- unique(c(pdat$pop_1, pdat$pop_2))
all_pops <- sort(all_pops)
pdat$pop_1 <- factor(pdat$pop_1, levels = all_pops)
pdat$pop_2 <- factor(pdat$pop_2, levels = all_pops)
ggplot(pdat, aes(x = pop_1, y = pop_2, color = fst, fill = fst)) + geom_tile() + theme_bw() + 
  scale_color_viridis_c(option = "inferno") +
  scale_fill_viridis_c(option = "inferno")
