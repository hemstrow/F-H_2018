library(dplyr)
#prep
plate_info <- read.table("Raw_data/plate_info.txt", stringsAsFactors = F)
bamlist <- read.table("Raw_data/sorted.bamlist.txt", stringsAsFactors = F)
bamlist <- bamlist[,c(5,9)]

bamlist$tag <- gsub("^.+RA_.{2}", "", bamlist$V9)
bamlist$tag <- substr(bamlist$tag, 1, 8)
bamlist$plate <- gsub("_.+$", "", bamlist$V9)
bamlist$plate <- gsub("merged.", "", bamlist$plate)
bamlist$plate[bamlist$plate == "SOMM270"] <- "plate1"

colnames(bamlist)[1:2] <- c("size", "bamfile")
colnames(plate_info) <- c("tag", "plate", "well", "ind", "pop")

comb <- merge(plate_info, bamlist, by = c("plate", "tag"))
comb <- comb[order(comb$size, decreasing = T),]

#remove duplicate HAW sample (although really we should just merge them...)
comb <- comb[-which(comb$tag == "GCCAAGAC" & comb$plate == "plate1"),]

comb[comb$pop %in% c("ENA", "WNA"),]$pop <- "NAM"

#get best n of each pop
n <- nrow(comb[comb$pop == "HAW",])

out <- comb %>%
  arrange_(~ desc(size)) %>%
  group_by_(~pop) %>%
  top_n(n = n, wt = size)

out <- as.data.frame(out)
pops <- unique(out$pop)

setwd("Raw_data/")
for(i in 1:length(pops)){
  write.table(out[out$pop == pops[i],]$bamfile, paste0(pops[i], "_ess_bamlist.txt"), row.names = F, col.names = F, quote = F, fileEncoding = "UTF-8")
  write.table(comb[comb$pop == pops[i],]$bamfile, paste0(pops[i], "_bamlist.txt"), row.names = F, col.names = F, quote = F, fileEncoding = "UTF-8")
}

#NOTE these may need to be re-encoded into a unix format using something like dos2unix