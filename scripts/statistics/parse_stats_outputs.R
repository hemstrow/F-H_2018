.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/4.1", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(data.table); library(ggplot2); library(snpR);


wb <- openxlsx::createWorkbook(creator = "William Hemstrom")


openxlsx::addWorksheet(wb, "Basic Stats")
openxlsx::addWorksheet(wb, "Fst")


#==================tajima's D===============
# pull in data
flist <- list.files(path = "./results/paralogs/", "tajimas_D")
out <- vector("list", length(flist))
for(i in 1:length(flist)){
  out[[i]]  <- data.table::fread(paste0("./results/paralogs/", flist[i]))
  out[[i]]$pop <- substr(flist[i], 27, 29)
}
tsd <- dplyr::bind_rows(out)
tsd <- tsd[,-1]

# get averages
aves <- tsd[tsd$snp.subfacet == ".OVERALL_MEAN",]
colnames(aves)[which(colnames(aves) == "n")] <- "thetas_n"

# # plot
# ggplot(aves, aes(y = log10(weighted_mean_ws.theta), x = log10(weighted_mean_ts.theta), label = pop)) + geom_text() +
#   theme_bw() + geom_abline(intercept = 0, slope = 1)


#==============dadi, gapped===============


dat <- data.table::fread("data/dadi_inputs/rand_10kgap_snps.txt")
samp.met <- colnames(dat)[-c(1:2)]
samp.met <- data.frame(ID = samp.met, pop = substr(samp.met, 1, 3))
dat <- import.snpR.data(dat[,-c(1:2)], dat[,1:2], samp.met)
dat <- filter_snps(dat, hwe = 0.000001, hwe_facets = "pop", min_loci = .75, maf_facets = "pop")

# ratio <- nrow(dat)/302446 # ratio of included snps
# L <-  (1373747 - 86)*ratio # approx number of considered bases, using x/1373747 (number of total sequenced bases) = 9370 (number of snps after LD gapping)/302446 (number of total snps)

dat <- calc_het_hom_ratio(dat, "pop")
dat <- calc_pi(dat, "pop")
dat <- calc_ho(dat, "pop")

ss <- get.snpR.stats(dat, "pop", c("pi", "ho"))$weighted.means
ss <- na.omit(ss)
ss <- reshape2::melt(ss[,c(2, 5, 6)], id.vars = c("subfacet"))
colnames(ss) <- c("pop", "stat", "value")

# ss$value <- ss$value/L # correct for non-polymorphic, considered sites.

# # plot
# ggplot(ss, aes(x = pop, y = value, fill = stat)) + 
#   geom_bar(stat = "identity", position = position_dodge()) +
#   theme_bw() + scale_fill_viridis_d(option = "inferno", end = 0.5)

# het/hom
samp.s <- get.snpR.stats(dat, "pop", "het_hom_ratio")$sample
samp.s <- na.omit(samp.s)
# ggplot(samp.s, aes(x = pop, y = `Het/Hom`)) + geom_point() + theme_bw()


#=======================tables==============

# basic stats:
## merge
comb.table <- merge(aves, reshape2::dcast(ss, pop~stat, value.var = "value"), by = "pop")
he_ho_means <- tapply(samp.s$`Het/Hom`, samp.s$pop, mean)
comb.table$het_hom <- he_ho_means
comb.table <- comb.table[,-c(2:5)]
n <- table(sample.meta(dat)$pop)
comb.table$n <- n[match(comb.table$pop, names(n))]


## clean
comb.table[, .SDcols = colnames(comb.table)[-1], colnames(comb.table)[-1] := lapply(.SD, round, digits = 3)]
comb.table <- comb.table[,-c(2:3)]
comb.table <- as.data.frame(comb.table)
colnames(comb.table) <- c("Population", "D", expression(n[D]), expression(pi), expression(H[o]), "Het/Hom", expression(n[other]))
facet.order <- c("NAM", "HAW", "GUA", "ROT", "SAI", "SAM", "FIJ", "NCA", "NOR", "QLD", "NSW", "VIC", "NZL")

## rename
comb.table$Population <- factor(comb.table$Population, facet.order)
comb.table <- dplyr::arrange(comb.table, Population)
key <- data.frame(code = facet.order, name = c("North America", "Hawaii", "Guam",
                                               "Rota", "Saipan", "Samoa", "Fiji", "New Caledonia", 
                                               "Norfolk Island", "Queensland", "New South Wales", "Victoria",
                                               "New Zealand"))
key$key <- paste0(key$name, " (", key$name, ")")
comb.table$Population <- key$key[match(comb.table$Population, key$code)]
# comb.table
openxlsx::writeData(wb, "Basic Stats", x = comb.table, keepNA = T)




#===================fst======================
dat <- calc_pairwise_fst(dat, facets = "pop", method = "genepop", boot = 1000, boot_par = 24)

fst <- get.snpR.stats(dat, "pop", "fst")

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

openxlsx::writeData(wb, "Fst", x = fst_matrix, keepNA = T)

#===============save and return============
openxlsx::saveWorkbook(wb, "./results/paralogs/statistics.xlsx", overwrite = TRUE)
