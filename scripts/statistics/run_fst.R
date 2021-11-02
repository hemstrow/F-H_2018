library(snpR)

input <- data.table::fread("../../../data/dadi_inputs/rand_10kgap_snps.txt")
samp.met <- colnames(input)[-c(1:2)]
samp.met <- data.frame(ID = samp.met, pop = substr(samp.met, 1, 3))
input <- import.snpR.data(input[,-c(1:2)], input[,1:2], samp.met)

input <- input[pop = c("NAM", "HAW", "GUA", "ROT", "SAI", "QLD", "NSW", "VIC")]

#=========fst===========
input <- filter_snps(input, hwe = 0.000001, hwe_facets = "pop", min_loci = .75, maf = 0.05, maf_facets = "pop")
input <- calc_pairwise_fst(input, facets = "pop", method = "genepop", boot = 1000, boot_par = 24)

 saveRDS(input, "../../../results/fst_out.RDS")
# input <- readRDS("../results/fst_out.RDS")

fst <- get.snpR.stats(input, "pop", "fst")

# condense fst results
fst <- fst$fst.matrix
fst_matrix <- fst$pop$fst
rn <- fst_matrix$p1
fst_matrix$p1 <- NA
colnames(fst_matrix)[1] <- "GUA"
fst_matrix <- as.matrix(fst_matrix)
fst_matrix <- rbind(fst_matrix, rep(NA, 8))
rn <- c(rn, "VIC")
cn <- colnames(fst_matrix)
fst_matrix <- matrix(as.numeric(fst_matrix), nrow(fst_matrix))

p_matrix <- fst$pop$p
p_matrix$p1 <- NA
colnames(p_matrix)[1] <- "GUA"
p_matrix <- as.matrix(p_matrix)
p_matrix <- rbind(p_matrix, rep(NA, 8))
rownames(p_matrix) <- rn
p_matrix <- t(p_matrix)
p_matrix <- matrix(as.numeric(p_matrix), nrow(p_matrix))


fst_matrix[is.na(fst_matrix)] <- p_matrix[is.na(fst_matrix)]
fst_matrix <- round(fst_matrix, 4)
colnames(fst_matrix) <- cn
rownames(fst_matrix) <- rn


