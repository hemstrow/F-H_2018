# runs cross validation on filtered snps. Arguments are: 1) which sample(s) should be used for cross validation? 2) what is the outfile name?

# get arguments, define paths, library snpR
args <- commandArgs(TRUE)

.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR)

cross_samples <- as.numeric(as.character(args[1]))
outfile <- as.character(args[2])


# read in distance data
dmeta <- read.table("~/monarch/github/F-H_2018/raw_data/migration_distances.csv", sep = ",", header = T, stringsAsFactors = F)


# read in genotypes
genos <- readRDS("~/monarch/github/F-H_2018/raw_data/flt_snps.RDS")


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
rm(genos, dmeta, sampmeta)


# run the cross_validation
cv_out <- cross_validate_genomic_prediction(dat, "migration_distance",
                                            iterations = 1000000,
                                            burn_in = 3000,
                                            thin = 300, cross_samples = cross_samples, plot = F)


# save output
out <- cbind(sampleID = rownames(cv_out$comparison), samplenum = paste0(cross_samples, collapse = "_"), cv_out$comparison, h2 = cv_out$model$h2)

write.table(out, outfile, quote = F, col.names = F, row.names = F, append = T)

outRDS <- paste0(gsub(".txt", "_", outfile), cross_samples, ".RDS")

saveRDS(cv_out, outRDS)
