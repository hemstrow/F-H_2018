setwd("Raw_data/")
#import bamlist and annotations
bamlist <- read.table("gbamlist.txt")
plate1 <- read.csv("plate1.csv", header = T)
plate2 <- read.csv("plate2_MFMB1.csv", header = T, skip = 1)
plate3 <- read.csv("plate3_MFMB2.csv", header = T, skip = 1)

DNA274 <- "ACAGTG"
DNA275 <- "CAGATC"

plate2 <- rbind(plate2[,1:3], setNames(plate2[,4:6], colnames(plate2)[1:3]))
plate3 <- rbind(plate3[,1:3], setNames(plate3[,4:6], colnames(plate3)[1:3]))


#add index 
barcodes <- read.table("RAD barcodes.txt", header = T)
plate2 <- merge(plate2, barcodes, by = "Well")
plate3 <- merge(plate3, barcodes, by = "Well")
plate1 <- merge(plate1, barcodes, by = "Well")


#get sample info
sampinfo <- data.frame(plate = substr(bamlist[,1], 8, 13), ind = substr(bamlist[,1], 20, 27), stringsAsFactors = FALSE)
sampinfo$plate[sampinfo$plate == "_split"] <- "plate1"



#add sample info to plates
plate2 <- merge(sampinfo[sampinfo$plate == DNA274,], plate2, by.x = "ind", by.y = "Index")
plate3 <- merge(sampinfo[sampinfo$plate == DNA275,], plate3, by.x = "ind", by.y = "Index")
plate1 <- merge(sampinfo[sampinfo$plate == "plate1",], plate1, by.x = "ind", by.y = "Index")

#combine plates in the correct order
combplates <- rbind(plate2, plate3, plate1)

#fix up some pop IDs
pops <- substr(combplates$Pop, 1, 3)
table(pops)
pops[pops == "Gua"] <- "GUA"
pops[pops == "New"] <- "NCA"
pops[pops == "Nor"] <- "NOR"
pops[pops == "NZ2"] <- "NZL"
pops[pops == "NZR"] <- "NZL"
pops[pops == "Rot"] <- "ROT"
pops[pops == "Sai"] <- "SAI"
pops[pops == "Sam"] <- "SAM"
table(pops)
combplates$Pop <- pops


#import raw genotypes.
raw_genos <- read.table("../../genotypes.geno", header = F, stringsAsFactors = F)

#associate sample IDs to genotypes. Works as long as the combplates are in the same order as the bamfile.
colnames(raw_genos) <- c("group", "position", paste0(combplates$Pop, "_", combplates$ID))

#add snp numbers
raw_genos <- cbind(snp = 1:nrow(raw_genos), raw_genos)

#sort by sample ID (and thus by population)
temp <- raw_genos[,4:ncol(raw_genos)]
raw_genos <- cbind(raw_genos[,1:3], temp[,order(colnames(temp))])


#pool SAM, NCA, FIJ; SAI, GUA; VIC, NSW, QLD, NZL; for filtering.


#grab pop info
pops <- substr(colnames(raw_genos)[4:ncol(raw_genos)], 1, 3)

#paste on pooling IDs
pops[pops %in% c("SAM", "NCA", "FIJ")] <- paste0("NPA_", pops[pops %in% c("SAM", "NCA", "FIJ")])
pops[pops %in% c("SAI", "GUA")] <- paste0("SAG_", pops[pops %in% c("SAI", "GUA")])
pops[pops %in% c("VIC", "NSW", "QLD", "NZL")] <- paste0("AUS_", pops[pops %in% c("VIC", "NSW", "QLD", "NZL")])
ncoln <- paste0(pops, "_", gsub(".+_", "", colnames(raw_genos)[-c(1:3)])) #new column names

#replace names
colnames(raw_genos) <- c(colnames(raw_genos)[1:3], ncoln)

#order
meta <- raw_genos[,1:3]
dat <- raw_genos[,4:ncol(raw_genos)]
dat <- dat[,-which(colnames(dat) == "HAW_9.1")] #remove the duplicate HAW sample
raw_genos <- cbind(meta, dat[,order(colnames(dat))])

#new pop IDs
l <- table(substr(colnames(raw_genos)[-c(1:3)], 1, 3))
l <- list(names(l), as.numeric(l))

#####################################################
#Filter and prepare data
library(snpR)

#all samples: for most analysis.
flt_genos <- filter_snps(raw_genos, 3, 0.05, 0.55, floor((ncol(raw_genos)-3)/2), pop = l) #just filtering snps.
#only well sequenced samples
flt_genos_clean <- filter_snps(raw_genos, 3, 0.05, 0.55, floor((ncol(raw_genos)-3)/2), .5, pop = l) #also remove crappy samples.

#reorder - full
meta <- flt_genos[,1:3]
dat <- flt_genos[,-c(1:3)]
colnames(dat)[grepl("^.+_.+_", colnames(dat))] <- gsub("^[A-Z]+_", "", colnames(dat)[grepl("^.+_.+_", colnames(dat))])
flt_genos <- cbind(meta, dat[,order(colnames(dat))])
#reorder - clean
meta <- flt_genos_clean[,1:3]
dat <- flt_genos_clean[,-c(1:3)]
colnames(dat)[grepl("^.+_.+_", colnames(dat))] <- gsub("^[A-Z]+_", "", colnames(dat)[grepl("^.+_.+_", colnames(dat))])
flt_genos_clean <- cbind(meta, dat[,order(colnames(dat))])

#check with a quick PCA
# pa <- format_snps(flt_genos_clean[1:10000,], 3, 7)
# pa <- cbind(pop = substr(pa$samp, 1, 3), pa)
# pca <- PCAfromPA(pa, 2)
#looks fine

#save
saveRDS(flt_genos, "flt_snps.RDS")
saveRDS(flt_genos, "flt_snps_clean.RDS")
