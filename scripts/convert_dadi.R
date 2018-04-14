library(snpR); library(readr)

#import data
anc <- read.table("../Data/anc_flank_20k.txt", header = F, stringsAsFactors = F)
anc <- anc[-1,]
ref <- read.table("../Data/ref_flank_20k.txt", header = F, stringsAsFactors = F)
ref <- ref[-1,]
dat <- readr::read_delim("../Data/rand_gap_snps_20k.txt", delim = "\t", col_names = T)

#bind together
x <- cbind(ref = ref, anc = anc, dat, stringsAsFactors = F)

##################
#get metadata
setwd("../Raw_data/")
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

colnames(x) <- c("ref", "anc", "chr", "position", paste0(combplates$Pop, "_", combplates$ID))

#merge ENA and WNA
colnames(x) <- gsub("WNA", "NAM", colnames(x))
colnames(x) <- gsub("ENA", "NAM", colnames(x))

######################

#sort and reformat
x <- x[,-ncol(x)]
pops <- table(substr(colnames(x[,5:ncol(x)]), 1, 3))
pops <- list(names(pops), as.numeric(pops))
meta <- x[,1:4]
dat <- x[,5:ncol(x)]
dat <- dat[,order(colnames(dat))]
x <- cbind(meta, dat)
#remove the fixed snps that are here for some reason... (low snp p value but not called in anything? wierd.)
fixed <- unlist(strsplit("6 205 452 878 898 977 1719 2060 2175 2357 2419 3162 3545 3766 3845 3897 3952 4004 4159 4293 4332 4792 4808 4985 5118 5186 5793 5919 6015 6123 6426 6636 6653 6677 6932 6993 7045 7090 7184 7219 7440 7464 8038 8138 8387 8579 8715 8885 9013 9151 9316 9506 9764 9839 9881",
                         split = " "))
fixed <- as.numeric(fixed)
x <- x[-fixed,]
#reformat
dadi_snps <- format_snps(x, 4, "dadi", pop = pops)

write.table(dadi_snps, "../Data/dadi_snps_20kgap.txt", quote = F, sep = "\t", col.names = T, row.names = F)
