# define paths, library snpR
.libPaths(c(.libPaths(), "/home/hemstrow/R/x86_64-pc-linux-gnu-library/3.6", "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library", "/share/apps/rmodules"))
library(snpR); library(ranger);

args <- commandArgs(TRUE)
run <- as.character(args[1])
outfile <- as.character(args[2])
cross_sample <- as.numeric(as.character(args[3]))

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
rm(genos, dmeta, sampmeta)
dat <- filter_snps(dat, 0.05, 0.55)

## remove the cross sample
rdat <- subset_snpR_data(dat, samps = (1:ncol(dat))[-cross_sample])

# run parameters
trim_cuttoffs <- c(1000, 200) # above the first level, will trim at the first percentage, below the last, will trim a single at a time
trim <- c(.9, .1) # trim at .9 before 1000 SNPs, at .1 before 200, and one per run below this.
response <- "migration_distance"
num.trees <- 50000 # number of trees
init.mtry <- nrow(rdat) # initial mtry
init.num.trees <- 50000 # number of trees in initial run. Probably doesn't need to be huge!
par <- 11

# run the first random forest
rf <- run_random_forest(rdat, response = response, mtry = init.mtry, num.trees = init.num.trees)

# run the refinement
refined_model <- refine_rf(rf = rf, 
                           response = response, 
                           trim_cuttoffs = trim_cuttoffs, 
                           num.trees = num.trees,
                           trim = trim,
                           par = par)

# save outputs
error.delta <- refined_model$error_delta
error.delta <- cbind(run = run, error.delta)
write.table(error.delta, paste0(outfile, "_delta_", run, ".txt"), quote = F, row.names = F, col.names = F)
saveRDS(refined_model, paste0(outfile, "_model_", run, ".RDS"))

# prediction, need to re-run the model since the impurity_corrected importance can apparently cause issues. Run with permutation.
kept.snps <- refined_model$best_model$data@snp.meta$.snp.id
kept.snps <- match(kept.snps, dat@snp.meta$.snp.id)

predict.rf <- run_random_forest(refined_model$best_model$data, response = response, mtry = nrow(refined_model$best_model$data), num.trees = num.trees, importance = "permutation", pvals = F)
sn <- format_snps(dat, "sn")
sn <- sn[,-c(1:3)]
sn <- sn[kept.snps, cross_sample]
sn <- as.data.frame(t(sn))
colnames(sn) <- refined_model$best_model$models$.base_.base$model$forest$independent.variable.names
out <- predict(predict.rf$models$.base_.base$model, sn)

write.table(data.frame(sample = cross_sample, phenotype = dat@sample.meta[cross_sample,response], prediction = out$predictions),
            paste0(outfile, "_prediction_", run, ".txt"), quote = F, col.names = F, row.names = F)
