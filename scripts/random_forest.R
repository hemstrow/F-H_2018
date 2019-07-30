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

# run the random forest
rf <- run_random_forest(dat, response = "migration_distance", mtry = nrow(dat)/4, num.trees = 5000, par = 11)
plot_manhattan(rf$data, "migration_distance_RF_importance", chr = "group")

source("../aus_monarchs/scripts/rf_refinement.R")

rf_refined <- refine_rf(rf, "migration_distance", 10000, trim = c(0.9, 0.1),
                        trim_cuttoffs = c(1000, 10), search_cuttoff = 2, 
                        par = 11)

saveRDS(rf_refined, "migration_distance/rf_refined.RDS")
saveRDS(rf_refined, "migration_distance/rf_raw.RDS")

plot_manhattan(rf_refined$best_model$data, "migration_distance_RF_importance", chr = "group", significant = 5e6)
pdat <- as.data.frame(rf_refined$error_delta)
ggplot(pdat, aes(x = log10(n_snps), y = prediction_error)) + geom_point()

ggplot(rf_refined$best_model$models$.base_.base$predictions, aes(x = predicted, y = pheno)) + 
  geom_point() + geom_smooth(method = "lm") +
