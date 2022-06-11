setwd("Raw_data/")

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


table(pops)
combplates$Pop <- pops

loc_meta <- read.table("sample_metadata.csv", header = TRUE, sep = ",")

combplates$samp <- paste0(combplates$Pop, "_", combplates$ID)
combplates$samp <- gsub("ENA", "NAM", combplates$samp)
combplates$samp <- gsub("WNA", "NAM", combplates$samp)
combplates <- merge(combplates, loc_meta[,-2], by = "samp")

addtl.meta <- data.table::fread("Sample_Inventory_Mol_Ecol.txt", header = TRUE)
addtl.meta$samp <- paste0(addtl.meta$Pop, "_", addtl.meta$ID)
addtl.meta$samp <- gsub("ENA", "NAM", addtl.meta$samp)
addtl.meta$samp <- gsub("WNA", "NAM", addtl.meta$samp)
combplates <- merge(combplates, addtl.meta[,-c(1:4)], by = "samp", all = TRUE)
combplates <- combplates[-which(duplicated(combplates)),]
combplates$Month <- substr(combplates$Month, 1, 3)
combplates$Day <- as.character(combplates$Day)
combplates$Day[nchar(combplates$Day) == 1 & !is.na(combplates$Day)] <- 
  paste0("0", combplates$Day[nchar(combplates$Day) == 1 & !is.na(combplates$Day)])
combplates$collection_date <- paste0(combplates$Day, "-", combplates$Month, "-", combplates$Year)
combplates$collection_date <- gsub(".+NA-", "", combplates$collection_date)
combplates$collection_date <- gsub("NA-", "", combplates$collection_date)

#=================biosample============
source_tab <- data.frame(Pop = unique(combplates$Pop),
                         Details = c("USA: California, ",
                                     "Mexico: ",
                                     "USA: Hawaii, ",
                                     "Australia: New South Wales, ",
                                     "Samoa: ",
                                     "Guam: ",
                                     "New Caledonia: ",
                                     "Northern Mariana Islands: Rota, ",
                                     "Northern Mariana Islands: Saipan, ",
                                     "Australia: Norfolk Island, ",
                                     "New Zealand: ",
                                     "Australia: Victoria, ",
                                     "Australia: Queensland, ",
                                     "Fiji: "))


sample_name <- paste0(combplates$samp, "_", combplates$ind, "_", combplates$plate)
organism <- "Danaus plexippus"
isolation_source <- paste0(combplates$Full_Location, ": ", combplates$lat, ",", combplates$long)
geo_loc_name <- paste0(source_tab$Details[match(combplates$Pop, source_tab$Pop)], isolation_source)
tissue <- "leg"
collection_date <- combplates$collection_date
isolate <- combplates$samp
isolate[isolate == "HAW_9"] <- paste0("HAW_9_aliquot_", 1:2)


biosample <- data.frame(sample_name = sample_name,
                        organism = organism,
                        isolate = isolate,
                        isolation_source = isolation_source,
                        geo_loc_name = geo_loc_name,
                        tissue = tissue,
                        collection_date = collection_date)
write.table(biosample, "../data/BioSample_meta.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

#================SRA===============
files <- read.table("all_fastqs_ncbi_list.txt")
files$V1 <- gsub(".+/", "", files$V1)
patterns <- gsub("rad_split_out_200_", "", files$V1)
patterns <- gsub("rad_split_out_199_", "", patterns)
patterns <- gsub("rad_split_out_199b_", "", patterns)
patterns <- gsub("_R.+_", "_", patterns)
patterns <- gsub(".fastq", "", patterns)
patterns <- gsub("SOMM270_split_", "SOMM270_", patterns)
patterns <- unique(patterns)
filenames <- vector("list", length(patterns))
names(filenames) <- patterns
p1 <- gsub("_.+", "_", patterns)
p2 <- gsub(".+_", "_", patterns)
SRA <- data.frame(library_ID = patterns)
SRA <- cbind(SRA, matrix("", nrow(SRA), 6))
colnames(SRA)[-c(1:2)] <- paste0("filename", 2:6)
colnames(SRA)[2] <- "filename"

for(i in 1:length(patterns)){
  tf <- unique(files[,1][which(grepl(p1[i], files[,1]) & grepl(p2[i], files[,1]))])
  SRA[i,2:(1 + length(tf))] <- tf
}


combplates$library_ID <- paste0(combplates$plate, "_", "GG", combplates$ind, "TGCAG")
combplates$library_ID <- gsub("plate1", "SOMM270", combplates$library_ID)
combplates$sample_name <- biosample$sample_name

SRA <- merge(SRA, combplates[,c("sample_name", "library_ID")], by = "library_ID")
SRA$title <- "RAD-Seq of Danaus Plexippus: Legs, samples from Pacific Expansion"
SRA$library_strategy <- "RAD-Seq"
SRA$library_source <- "GENOMIC"
SRA$library_selection <- "Reduced Representation"
SRA$library_layout <- "paired"
SRA$platform <- "ILLUMINA"
SRA$instrument_model <- "Illumina HiSeq 4000"
SRA$design_description <- "RAD-Seq according to Ali et al. (2016)"
SRA$filetype <- "fastq"

SRA <- SRA[,c(1, 8:ncol(SRA), 2:7)]
write.table(SRA, "../data/SRA_metadata.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
