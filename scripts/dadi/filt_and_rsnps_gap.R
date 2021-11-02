library(snpR); library(plyr); library(data.table)
#set parameters:

gap <- 10000 #set gap size
idat <- "~/monarch/github/F-H_2018/raw_data/nomaf_raw_genotypes.geno" #set input data, in this case just the raw genotypes but with metadata column names.
odat <- "~/monarch/github/F-H_2018/data/dadi_inputs/rand_10kgap_snps.txt" #set output name
hfhf <- 0.55 #high frequency heterozygote filter cutoff
min_ind <- 50 #min number of sequenced individuals per loci
meta.dir <- "../../raw_data/" #directory where metadata files are stored.


############setup data with metadata#####################

setwd(meta.dir)
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
pops[pops == "ENA"] <- "NAM"
pops[pops == "WNA"] <- "NAM" #pool ENA and WNA
combplates$Pop <- pops

#import raw genotypes.
raw_genos <- read.table(idat, header = F, stringsAsFactors = F)

#associate sample IDs to genotypes. Works as long as the combplates are in the same order as the bamfile.
colnames(raw_genos) <- c("group", "position", paste0(combplates$Pop, "_", combplates$ID))

#sort by sample ID (and thus by population)
temp <- raw_genos[,4:ncol(raw_genos)]
raw_genos <- cbind(raw_genos[,1:3], temp[,order(colnames(temp))])
raw_genos <- raw_genos[,-which(colnames(raw_genos) == "HAW_9.1")] #remove the duplicate HAW sample


#####################################################
#Filter and prepare data

#only filtering high frequency hets and poorly sequenced SNPs, for dadi.
dadi_genos <- filter_snps(raw_genos, 2, FALSE, hfhf, min_ind)

###############################
#rgap it
dadi_genos <- rgap_snps(dadi_genos, gap) #get snps seperated by a mininum gap length

##############################
#save
data.table::fwrite(dadi_genos, odat, quote = FALSE, col.names = T, sep = "\t", row.names = F)
