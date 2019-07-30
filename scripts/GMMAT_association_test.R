dmeta <- read.table("raw_data/migration_distances.csv", sep = ",", header = T, stringsAsFactors = F)
dat <- readRDS("raw_data/FHmon_geno_snpR.RDS")


# subset and add migration data
dat <- subset_snpR_data(dat, facets = "pop", subfacets = c("ENA", "WNA"))

sampmeta <- data.frame(sampleID = dat@sample.meta$samp, pop = dat@sample.meta$pop, stringsAsFactors = F)
sampmeta$sampleID <- tolower(sampmeta$sampleID)
sampmeta$sampleID <-  gsub(".+_", "", sampmeta$sampleID)
sampmeta$ord <- 1:nrow(sampmeta)

dmeta$Sample_ID <- tolower(dmeta$Sample_ID)
dmeta$Sample_ID <- gsub("ellwood", "goleta", dmeta$Sample_ID)

new.meta <- merge(sampmeta, dmeta[,-2], by.x = "sampleID", by.y = "Sample_ID", sort = F)
new.meta <- new.meta[order(new.meta$ord),]

dat <- import.snpR.data(dat[,new.meta$ord], snp.meta = dat@snp.meta, sample.meta = new.meta[,-3])
dat <- filter_snps(dat, 0.05, 0.55, .5, .5)

dat <- calc_association(dat, response = "migration_distance", formula = "migration_distance ~ pop")
plot_manhattan(dat, "gmmat_pval_migration_distance", chr = "group", log.p = T, significant = 0.00005, suggestive = 0.0005)
