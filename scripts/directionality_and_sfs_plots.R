# remotes::install_github("hemstrow/snpR", ref = "dev")
library(snpR)

x <- read.table("data/dadi_inputs/dadi_10kgap_snps.txt", header = T, stringsAsFactors = F)

# the make_SFS function arguments are: input data file, populations, projections (the number of
# sample gene copies (2 per individual!) to make the sfs from, lower numbers will keep more snps, maximizing snps
# kept over individuals kept tends to produce better results as long as you don't get TOO small),
# and then a logical for if the sfs should be folded or not. We can't fold the sfs if we are going to
# use the directionality index

# NAM, HAW
sf <- make_SFS(x, c("NAM", "HAW"), c(100, 10), F)
plot_sfs(sf)
calc_directionality(sf)

# HAW, GUA
sf <- make_SFS(x, c("HAW", "GUA"), c(10, 15), F)
plot_sfs(sf)
calc_directionality(sf)

# GUA, ROT
sf <- make_SFS(x, c("GUA", "ROT"), c(15, 15), F)
plot_sfs(sf)
calc_directionality(sf)

# HAW, QLD
sf <- make_SFS(x, c("HAW", "QLD"), c(10, 10), F)
plot_sfs(sf)
calc_directionality(sf)

# HAW, NOR
sf <- make_SFS(x, c("HAW", "NOR"), c(10, 10), F)
plot_sfs(sf)
calc_directionality(sf)
