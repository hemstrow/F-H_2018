# run parameters
response <- "migration_distance"
mtry.scale <- 1 # percentage of snps to use in mtry
num.trees <- 100000 # number of trees in initial run. Probably doesn't need to be huge!
par <- F
trim <- .9



# define paths, library snpR
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR); library(ranger);

args <- commandArgs(TRUE)
run <- as.character(args[1])
outfile <- as.character(args[2])
cross_sample <- as.numeric(as.character(args[3]))
sub_pop <- as.character(args[4])

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
dmeta$pop[dmeta$pop == "western"] <- "WNA"
dmeta$pop[dmeta$pop == "eastern"] <- "ENA"


new.meta <- merge(sampmeta, dmeta, by.x = c("sampleID", "pop"), by.y = c("Sample_ID", "pop"))

new.meta <- new.meta[order(new.meta$ord),]
dat <- import.snpR.data(genos[,new.meta$ord], snp.meta = meta, sample.meta = new.meta, mDat = "NN")
dat <- subset_snpR_data(dat, facets = "pop", subfacets = sub_pop)
rm(genos, dmeta, sampmeta)
dat <- filter_snps(dat, 0.05, 0.55)

## remove the cross sample
rdat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-cross_sample])

# run the first random forest
rf <- run_random_forest(rdat, response = response, mtry = nrow(rdat)*mtry.scale, num.trees = num.trees, par = par)

# run the refinement
best.imp <- quantile(abs(rf$models$.base_.base$model$variable.importance), trim)
best.imp <- which(rf$models$.base_.base$model$variable.importance >= best.imp[1])
rdat <- subset_snpR_data(rdat, snps = best.imp)
rf <- run_random_forest(rdat, response = response, mtry = nrow(rdat)*mtry.scale, num.trees = num.trees, 
                        par = par, importance = "permutation", pvals = F)



# save outputs
saveRDS(rf, paste0(outfile, "_model_", run, ".RDS"))

## prediction, need to re-run the model since the impurity_corrected importance can apparently cause issues. Run with permutation.
kept.snps <- rf$data@snp.meta$.snp.id
kept.snps <- match(kept.snps, dat@snp.meta$.snp.id)

sn <- format_snps(dat, "sn")
sn <- sn[,-c(1:3)]
sn <- sn[kept.snps, cross_sample]
sn <- as.data.frame(t(sn))
colnames(sn) <- rf$models$.base_.base$model$forest$independent.variable.names
out <- predict(rf$models$.base_.base$model, sn)

write.table(data.frame(sample = cross_sample, phenotype = dat@sample.meta[cross_sample,response], prediction = out$predictions),
            paste0(outfile, "_prediction_pop_", sub_pop, "_", run, ".txt"), quote = F, col.names = F, row.names = F)
