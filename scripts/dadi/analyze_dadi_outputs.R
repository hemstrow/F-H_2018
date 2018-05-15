library(ggplot2); library(dplyr)
#==========import and prep data==========

res <- readLines("data/dadi_inputs/1st_pass_HGR_out.txt")

#grab the data
rdf <- data.frame(model = character(length(res)), pops = character(length(res)), theta = numeric(length(res)), ll = numeric(length(res)), AIC = numeric(length(res)), mnum = numeric(length(res)), parms = numeric(length(res)), stringsAsFactors = F)
for(i in 1:length(res)){
  tr <- strsplit(res[i], "\t")
  tr <- unlist(tr)
  rdf[i,] <- c(tr[c(2,4,6,8,10, 12)], i)
}

#define the parameters per model:
mplist <- list(
  p2_cgrowth = c("nu1B", "nu2B", "nu1f", "nu2f", "ts", "tp", "m12", "m21"), 
  p2_lgrowth_both = c("nu1B", "nu2B", "K1", "K2", "ts", "tp", "m12", "m21", "r1", "r2"),
  p2_lgrowth_1 = c("nu1B", "nu2B", "K1", "nu2f", "ts", "tp", "m12", "m21", "r1"),
  p2_lgrowth_2 = c("nu1B", "nu2B", "nu1f", "K2", "ts", "tp", "m12", "m21", "r2")
)

#what are the unique combinations of model and pop?
facets <- unique(rdf[,1:2])
wlist <- vector("list", nrow(facets))
names(wlist) <- paste0(facets[,1], "_", gsub(" ", "_", facets[,2]))

#get the results, seperated by the facets.
for(i in 1:length(wlist)){
  tdat <- which(apply(rdf, 1, function(x) all(as.character(x[1:2]) == facets[i,])))
  tdat <- as.data.frame(rdf[tdat,], stringsAsFactors = F)
  tdat[,3] <- as.numeric(tdat[,3])
  tdat[,4] <- as.numeric(tdat[,4])
  tdat[,5] <- as.numeric(tdat[,5])
  parms <- strsplit(tdat[,ncol(tdat) - 1], " ")
  po <- matrix(NA, length(parms), length(unlist(parms[1])))
  for(j in 1:length(parms)){
    po[j,] <- as.numeric(unlist(parms[[j]]))
  }
  colnames(po) <- mplist[[facets$model[i]]]
  po <- cbind(tdat[,3:5], po)
  colnames(po)[1:3] <- c("theta", "ll", "AIC")
  wlist[[i]] <- po
}


#==========analysis========
#overall between models.
rdf$AIC <- as.numeric(rdf$AIC)
rdf$ll <- as.numeric(rdf$ll)
ggplot(rdf, aes(x = AIC)) + geom_histogram() + facet_grid(model ~ pops) + theme_bw()

#best replicate per model
best.reps <- rdf %>% group_by(model) %>% group_by(pops, add = TRUE) %>% top_n(-1, AIC)
best.reps
best.reps[,-c(2,3,4,5,7)] #where are we pushing parameter bounds?
ggplot(best.reps, aes(y = AIC, x = model, color = pops)) + geom_point() + theme_bw()

#these best.reps are the starting points for run 2