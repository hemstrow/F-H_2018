library(snpR)

# read in distance data
dmeta <- read.table("raw_data/migration_distances.csv", sep = ",", header = T, stringsAsFactors = F)

# read in genotypes
genos <- readRDS("Raw_data/flt_snps.RDS")

# prepare a new sample meta file
pops <- colnames(genos)[-c(1:3)]
pops <- gsub("_.+", "", pops)

meta <- genos[,c(1:3)]

genos <- genos[,-c(1:3)]
genos <- genos[,which(pops %in% c("ENA", "WNA"))]

pops <- pops[which(pops %in% c("ENA", "WNA"))]

sampIDs <- colnames(genos)
sampIDs <- gsub(".+_", "", sampIDs)


sampmeta <- data.frame(pop = pops, sampleID = sampIDs, stringsAsFactors = F)
sampmeta$sampleID <- tolower(sampmeta$sampleID)
sampmeta$ord <- 1:nrow(sampmeta)

dmeta$Sample_ID <- tolower(dmeta$Sample_ID)
dmeta$Sample_ID <- gsub("ellwood", "goleta", dmeta$Sample_ID)

new.meta <- merge(sampmeta, dmeta, by.x = "sampleID", by.y = "Sample_ID")

new.meta <- new.meta[order(new.meta$ord),]
dat <- import.snpR.data(genos[,new.meta$ord], snp.meta = meta, sample.meta = new.meta, mDat = "NN")

sn <- format_snps(dat, "sn")

tsn <- t(sn[,-c(1:3)])

tsn <- cbind(new.meta, tsn)

rm(genos, dmeta, sampmeta)

# run BGLR
phenotypes <- new.meta$migration_distance
tsn <- tsn[,-c(1:5)]
colnames(tsn) <- paste0("m", 1:ncol(tsn)) # marker names
rownames(tsn) <- paste0("s", 1:nrow(tsn)) # ind IDS
rm(sn)
gc();gc();gc();gc();

# prepare ETA
# need to impute missing data, use a binomial draw off of the allele frequency or just stick in 2*af.
ETA <- list(list(X = tsn, model = "BayesB", saveEffects = T))


BGLR_mod <- BGLR::BGLR(y = phenotypes, ETA = ETA, nIter = 10000, burnIn = 1000, thin = 100)

# grab h2 estimate
B <- BGLR::readBinMat('ETA_1_b.bin')
h2 <- rep(NA,nrow(B))
varU <- h2
varE <- h2
for(i in 1:length(h2)){
  u <- ind.genos%*%B[i,]	
  varU[i] <- var(u)
  varE[i] <- var(phenotypes-u)
  h2[i] <- varU[i]/(varU[i] + varE[i])
}
h2 <- mean(h2)