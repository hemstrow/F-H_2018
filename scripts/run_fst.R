library(snpR)
input <- readRDS("../results/paralogs/nomaf_paralog_fix_snpR.RDS")
meta <- read.csv("../Raw_data/sample_metadata.csv")
sample.meta(input) <- meta

input <- import.snpR.data(genotypes(input), snp.meta = snp.meta(input), sample.meta = sample.meta(input))
input <- input[pop = c("NAM", "HAW", "GUA", "ROT", "SAI", "QLD", "NSW", "VIC")]

#=========fst===========
input <- filter_snps(input, maf = 0.05, maf.facets = "pop", min_loci = 0.75)
input <- calc_pairwise_fst(input, facets = "pop", method = "wc", boot = 1000, boot_par = 24)

saveRDS(input, "../results/fst_out.RDS")
