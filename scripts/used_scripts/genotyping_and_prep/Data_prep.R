.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR); library(readr)


######################################################################
# prep metadata
setwd("../github/F-H_2018/Raw_data/")


# import bamlist and annotations
bamlist <- read.table("gbamlist.txt")
plate1 <- read.csv("plate1.csv", header = T)
plate2 <- read.csv("plate2_MFMB1.csv", header = T, skip = 1)
plate3 <- read.csv("plate3_MFMB2.csv", header = T, skip = 1)

DNA274 <- "ACAGTG"
DNA275 <- "CAGATC"

plate2 <- rbind(plate2[,1:3], setNames(plate2[,4:6], colnames(plate2)[1:3]))
plate3 <- rbind(plate3[,1:3], setNames(plate3[,4:6], colnames(plate3)[1:3]))


# add index 
barcodes <- read.table("RAD barcodes.txt", header = T)
plate2 <- merge(plate2, barcodes, by = "Well")
plate3 <- merge(plate3, barcodes, by = "Well")
plate1 <- merge(plate1, barcodes, by = "Well")


# get sample info
sampinfo <- data.frame(plate = substr(bamlist[,1], 8, 13), ind = substr(bamlist[,1], 20, 27), stringsAsFactors = FALSE)
sampinfo$plate[sampinfo$plate == "_split"] <- "plate1"


# add sample info to plates
plate2 <- merge(sampinfo[sampinfo$plate == DNA274,], plate2, by.x = "ind", by.y = "Index")
plate3 <- merge(sampinfo[sampinfo$plate == DNA275,], plate3, by.x = "ind", by.y = "Index")
plate1 <- merge(sampinfo[sampinfo$plate == "plate1",], plate1, by.x = "ind", by.y = "Index")

# combine plates in the correct order
combplates <- rbind(plate2, plate3, plate1)

# fix up some pop IDs
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
pops[pops == "WNA"] <- "NAM"
pops[pops == "ENA"] <- "NAM"

table(pops)
combplates$Pop <- pops

setwd("../../../paralogue_check/")

#######################################################
# import data and associate metadata
cat("Importing data.\n")

# import raw genotypes.
raw_genos <- read.table("all_sites_paralogs.geno", sep = "\t", header = F, stringsAsFactors = F)
raw_genos <- raw_genos[,-ncol(raw_genos)] # strip NA column


cat("Read in ", nrow(raw_genos), " sites. Sorting data.\n")
# associate sample IDs to genotypes. Works as long as the combplates are in the same order as the bamfile.
colnames(raw_genos) <- c("group", "position", paste0(combplates$Pop, "_", combplates$ID))

# add snp numbers
raw_genos <- cbind(snp = 1:nrow(raw_genos), raw_genos)

# sort by sample ID (and thus by population)
temp <- raw_genos[,4:ncol(raw_genos)]
raw_genos <- cbind(raw_genos[,1:3], temp[,order(colnames(temp))])

 
# grab pop info
pops <- substr(colnames(raw_genos)[4:ncol(raw_genos)], 1, 3)

# order
meta <- raw_genos[,1:3]
dat <- raw_genos[,4:ncol(raw_genos)]
dat <- dat[,-which(colnames(dat) == "HAW_9.1")] #remove the duplicate HAW sample
raw_genos <- cbind(meta, dat[,order(colnames(dat))])
rm(dat, temp); gc(); gc();
 
# gather metadata and create snpRdata object
cat("Creating snpRdata object.\n")
sample.meta <- as.data.frame(cbind(samp = colnames(raw_genos)[-c(1:3)], pop = substr(colnames(raw_genos)[-c(1:3)], 1, 3)),
                              stringsAsFactors = F)
snp.meta <- raw_genos[,2:3]
cat("snp.meta rows: ", nrow(snp.meta), "\nraw geno rows:", nrow(raw_genos), "\n")

dat <- import.snpR.data(raw_genos[,-c(1:3)], snp.meta, sample.meta)
cat("dat rows:", nrow(dat), "\n")


#########################################################
# all sites
cat("Saving all sites data.\n")
saveRDS(dat, "allsites_paralog_fix_snpR.RDS")

#########################################################
# maf
cat("Filtering by minor allele frequency.\n")
dat_maf <- filter_snps(dat, maf = 0.05)
saveRDS(dat_maf, "maf_paralog_fix_snpR.RDS") 
rm(dat_maf)

#########################################################
# no maf, but only polymorphic
cat("Filtering out only non-polymorphic sites.\n")
dat_nomaf <- filter_snps(dat)
saveRDS(dat_nomaf, "nomaf_paralog_fix_snpR.RDS")

