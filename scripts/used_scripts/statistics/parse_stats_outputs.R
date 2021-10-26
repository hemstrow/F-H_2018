library(data.table); library(ggplot2); library(snpR);



#==================tajima's D===============
# pull in data
flist <- list.files(path = "./results/paralogs/", "tajimas_D")
out <- vector("list", length(flist))
for(i in 1:length(flist)){
  out[[i]]  <- data.table::fread(paste0("./results/paralogs/", flist[i]))
  out[[i]]$pop <- substr(flist[i], 27, 29)
}
tsd <- dplyr::bind_rows(out)

# purge off odd values
tsd[is.nan(tsd$D)]$D <- NA
tsd[is.infinite(tsd$D)]$D <- NA

# get averages
tab.D <- data.table::dcast(tsd, 1 ~ pop, value.var = "D", 
                           fun.aggregate = mean, na.rm = T)
tab.ws <- data.table::dcast(tsd, 1 ~ pop, value.var = "ws.theta", fun.aggregate = mean, na.rm = T)
tab.ts <- data.table::dcast(tsd, 1 ~ pop, value.var = "ts.theta", fun.aggregate = mean, na.rm = T)

# plot
## thetas
pdat.t <- data.frame(population = c(names(tab.ws), names(tab.ts)),
                   stat = c(rep("watersons_theta", length(tab.ws)), rep("tajimas_theta", length(tab.ts))),
                   value = c(as.numeric(tab.ws), as.numeric(tab.ts)))
pdat.t <- pdat.t[-(which(pdat.t$population == ".")),]
ggplot(pdat.t, aes(x = population, y = value, fill = stat)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() + scale_fill_viridis_d(option = "inferno", end = 0.5)

## D
pdat.D <- data.frame(population = names(tab.D), D = as.numeric(tab.D))
pdat.D <- pdat.D[-which(pdat.D$population == "."),]
ggplot(pdat.D, aes(x = population, y = D)) + geom_point() + theme_bw() +
  geom_hline(yintercept = 0, color = "red")


#==============take two, paralogs===============
input <- readRDS("results/paralogs/nomaf_paralog_fix_snpR.RDS")
dat <- filter_snps(input, min_loci = 0.75)
dat <- calc_het_hom_ratio(dat)
dat <- calc_pi(dat, "pop")
dat <- calc_ho(dat, "pop")

ss <- get.snpR.stats(dat, "pop")
colnames(ss)[2] <- "pop"
tab.pi <- data.table::dcast(ss, 1 ~ pop, value.var = "pi", fun.aggregate = mean, na.rm = T)
tab.ho <- data.table::dcast(ss, 1 ~ pop, value.var = "ho", fun.aggregate = mean, na.rm = T)

# plot
pdat.div <- data.frame(population = c(names(tab.pi), names(tab.ho)),
                   stat = c(rep("pi", length(tab.pi)), rep("ho", length(tab.ho))),
                   value = c(as.numeric(tab.pi), as.numeric(tab.ho)))
pdat.div <- pdat.div[-(which(pdat.div$population == "1")),]
ggplot(pdat.div, aes(x = population, y = value, fill = stat)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() + scale_fill_viridis_d(option = "inferno", end = 0.5)

# het/hom
samp.s <- dat@sample.stats
colnames(samp.s)[5] <- "Het_Hom"
ggplot(samp.s, aes(x = pop, y = Het_Hom)) + geom_point() + theme_bw()


# fst
dat.maf <- filter_snps(dat, maf = 0.05)
dat.maf <- calc_pairwise_fst(dat.maf, "pop", method = "Genepop")


#=======================tables==============

# basic stats:
comb.table <- merge(pdat.D, reshape2::dcast(pdat.div, population~stat), by = "population")
he_ho_means <- tapply(samp.s$Het_Hom, samp.s$pop, mean)
comb.table$het_hom <- he_ho_means
comb.table[,-1] <- round(comb.table[,-1], 3)
colnames(comb.table) <- c("population", "Tajima's D", expression(H[o]), "pi", "Het/Hom")

facet.order <- c("NAM", "HAW", "GUA", "ROT", "SAI", "SAM", "FIJ", "NCA", "NOR", "QLD", "NSW", "VIC", "NZL")
comb.table$population <- factor(comb.table$population, facet.order)
comb.table <- dplyr::arrange(comb.table, population)
#formattable(comb.table, align = c("l", rep("c", 4)))
comb.table

# fst
fst <- dat.maf[[2]]
fst <- data.frame(pop_1 = substr(fst[,1], 1, 3), pop_2 = substr(fst[,1], 5, 7), fst = fst[,2], stringsAsFactors = F)
opts <- unique(c(fst$pop_1, fst$pop_2))
opts <- opts[order(match(opts, facet.order))]
fst$pop_1 <- factor(fst$pop_1, opts, ordered = T)
fst$pop_2 <- factor(fst$pop_2, opts, ordered = T)
fst$fst <- round(fst$fst, 3)
for(i in 1:nrow(fst)){
  if(fst$pop_1[i] < fst$pop_2[i]){
    np1 <- fst$pop_2[i]
    np2 <- fst$pop_1[i]
    fst$pop_1[i] <- np1
    fst$pop_2[i] <- np2
  }
}

fst <- reshape2::dcast(fst, pop_1~pop_2)
fst[is.na(fst)] <- ""
colnames(fst)[1] <- "population"
fst
