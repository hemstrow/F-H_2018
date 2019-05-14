library(ggplot2); library(dplyr); source("scripts/dadi/make_dadi_parmfile_from_results.R")
#==========import and prep data==========

res <- readLines("data/dadi_inputs/cat_NH_2nd_pass_out_dportik.txt")

#grab the data
rdf <- data.frame(model = character(length(res)), pops = character(length(res)), theta = numeric(length(res)), ll = numeric(length(res)), AIC = numeric(length(res)), mnum = numeric(length(res)), parms = numeric(length(res)), stringsAsFactors = F)
for(i in 1:length(res)){
  tr <- strsplit(res[i], "\t")
  tr <- unlist(tr)
  rdf[i,] <- c(tr[c(2,4,6,8,10, 12)], i)
}

#define the parameters per model:
mplist <- list(
  vic_no_mig = c("nuA", "Ti", "s"),
  
  vic_anc_asym_mig = c("nuA", "m12", "m21", "T1", "T2", "s"),
  
  vic_sec_contact_asym_mig = c("nuA", "m12", "m21", "T1", "T2", "s"),
  
  vic_no_mig_admix_early = c("nuA", "Ti", "s", "f"),
  
  vic_no_mig_admix_late = c("nuA", "Ti", "s", "f"),
  
  vic_two_epoch_admix = c("nuA", "T1", "T2", "s", "f"),
  
  
  # exponential models
  founder_nomig_growth_pop_2 =  c("nuA", "nu2", "Ti", "s"),
  
  founder_sym_growth_pop_2 = c("nuA", "nu2", "m", "Ti", "s"),
  
  founder_asym_growth_pop_2 = c("nuA", "nu2", "m12", "m21", "Ti", "s"),
  
  founder_nomig_admix_early_growth_pop_2 = c("nuA", "nu2", "Ti", "s", "f"),
  
  founder_nomig_admix_late_growth_pop_2 = c("nuA", "nu2", "Ti", "s", "f"),
  
  founder_nomig_admix_two_epoch_growth_pop_2 = c("nuA", "nu2", "T1", "T2", "s", "f"),
  
  
  founder_nomig_growth_pop_1 = c("nuA", "nu1", "Ti", "s"),
  
  founder_sym_growth_pop_1 = c("nuA", "nu1", "m", "Ti", "s"),
  
  founder_asym_growth_pop_1 = c("nuA", "nu1", "m12", "m21", "Ti", "s"),
  
  founder_nomig_admix_early_growth_pop_1 = c("nuA", "nu1", "Ti", "s", "f"),
  
  founder_nomig_admix_late_growth_pop_1 =  c("nuA", "nu1", "Ti", "s", "f"),
  
  founder_nomig_admix_two_epoch_growth_pop_1 = c("nuA", "nu1", "T1", "T2", "s", "f"),
  
  
  founder_nomig_growth_both = c("nuA", "nu1", "nu2", "Ti", "s"),
  
  founder_sym_growth_both = c("nuA", "nu1", "nu2", "m", "Ti", "s"),
  
  founder_asym_growth_both = c("nuA", "nu1", "nu2", "m12", "m21", "Ti", "s"),
  
  founder_nomig_admix_early_growth_both = c("nuA", "nu1", "nu2", "Ti", "s", "f"),
  
  founder_nomig_admix_late_growth_both = c("nuA", "nu1", "nu2", "Ti", "s", "f"),
  
  founder_nomig_admix_two_epoch_growth_both = c("nuA", "nu1", "nu2", "T1", "T2", "s", "f"),
  
  
  # logistic models
  founder_nomig_logistic_pop_2 = c("nuA", "K2", "r2", "Ti", "s"),
  
  founder_sym_logistic_pop_2 = c("nuA", "K2", "r2", "m", "Ti", "s"),
  
  founder_asym_logistic_pop_2 = c("nuA", "K2", "r2", "m12", "m21", "Ti", "s"),
  
  founder_nomig_admix_early_logistic_pop_2 = c("nuA", "K2", "r2", "Ti", "s", "f"),
  
  founder_nomig_admix_late_logistic_pop_2 = c("nuA", "K2", "r2", "Ti", "s", "f"),
  
  founder_nomig_admix_two_epoch_logistic_pop_2 = c("nuA", "K2", "r2", "T1", "T2", "s", "f"),
  
  
  founder_nomig_logistic_pop_1 = c("nuA", "K1", "r1", "Ti", "s"),
  
  founder_sym_logistic_pop_1 = c("nuA", "K1", "r1", "m", "Ti", "s"),
  
  founder_asym_logistic_pop_1 = c("nuA", "K1", "r1", "m12", "m21", "Ti", "s"),
  
  founder_nomig_admix_early_logistic_pop_1 = c("nuA", "K1", "r1", "Ti", "s", "f"),
  
  founder_nomig_admix_late_logistic_pop_1 = c("nuA", "K1", "r1", "Ti", "s", "f"),
  
  founder_nomig_admix_two_epoch_logistic_pop_1 = c("nuA", "K1", "r1", "T1", "T2", "s", "f"),
  
  
  founder_nomig_logistic_both = c("nuA", "K1", "K2", "r1", "r2", "Ti", "s"),
  
  founder_sym_logistic_both = c("nuA", "K1", "K2", "r1", "r2", "m", "Ti", "s"),
  
  founder_asym_logistic_both = c("nuA", "K1", "K2", "r1", "r2", "m12", "m21", "Ti", "s"),
  
  founder_nomig_admix_early_logistic_both = c("nuA", "K1", "K2", "r1", "r2", "Ti", "s", "f"),
  
  founder_nomig_admix_late_logistic_both = c("nuA", "K1", "K2", "r1", "r2", "Ti", "s", "f"),
  
  founder_nomig_admix_two_epoch_logistic_both = c("nuA", "K1", "K2", "r1", "r2", "T1", "T2", "s", "f")
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
  wlist[[i]] <- cbind(pops = tdat$pops, model = tdat$model, po, stringsAsFactors = F)
}


#==========analysis========
# remove any model results that have extremely low AIC results, likely caused by integration errors
rdf$AIC <- as.numeric(rdf$AIC)
rdf$ll <- as.numeric(rdf$ll)
bads <- which(rdf$AIC < 0.1*mean(rdf$AIC))
if(length(bads) > 0){
  for(i in 1:length(bads)){
    target.wlist <- which(names(wlist) == paste0(rdf$model[bads[i]], "_", paste(unlist(strsplit(rdf$pops[bads[i]], " ")), collapse = "_")))
    target.run <- which(wlist[[target.wlist]]$AIC == rdf$AIC[bads[i]])
    wlist[[target.wlist]] <- wlist[[target.wlist]][]
  }
  rdf <- rdf[-bads,]
}





# overall between models.
rdf$mod_median <- aggregate(AIC ~ model, data = rdf, FUN = median)[match(rdf$model, sort(unique(rdf$model))),2] # add model means
(table(rdf$model)[which.min(table(rdf$model))]) # which model got through the least number of runs?
# ggplot(rdf, aes(x = AIC)) + geom_histogram() + facet_grid(model ~ pops) + theme_bw() + 
#   xlim(c(min(rdf$AIC), 10000)) + geom_vline(aes(xintercept = mod_median), color = "red")
ggplot(rdf, aes(x = AIC)) + geom_histogram() + facet_wrap(~model) + theme_bw() + 
  xlim(c(min(rdf$AIC), 10000)) + geom_vline(aes(xintercept = mod_median), color = "red")
ggplot(rdf, aes(x = model, y = log(AIC), color = pops)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(rdf, aes(x = model, y = log(AIC), color = pops)) + geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


# which model had the best AIC on average?
best.mods <- arrange(unique(rdf[,c(1,2,ncol(rdf))]), mod_median)
best.mods

# best replicate per model
best.reps <- rdf %>% group_by(model) %>% group_by(pops, add = TRUE) %>% top_n(-1, AIC)
best.reps <- arrange(best.reps, pops, model)
best.reps[,-c(3,4,5,7)] #where are we pushing parameter bounds?
ggplot(best.reps, aes(y = AIC, x = model, color = pops)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
best.reps[order(best.reps$AIC),-c(3,4,7)]



# these best.reps are the starting points for next run


# plots for the parameters
ggplot(wlist$founder_asym_growth_both_NAM_HAW, aes(x = nuA, y = nu1, color = AIC)) + geom_point() + theme_bw() + scale_color_viridis_c()
ggplot(wlist$founder_asym_growth_both_NAM_HAW, aes(x = nuA, y = Ti, color = AIC)) + geom_point() + theme_bw() + scale_color_viridis_c()
ggplot(wlist$founder_asym_growth_both_NAM_HAW, aes(x = m12, y = m21, color = AIC)) + geom_point() + theme_bw() + scale_color_viridis_c()
ggplot(wlist$founder_asym_growth_both_NAM_HAW, aes(x = nu1, y = s, color = log(theta))) + geom_point() + theme_bw() + scale_color_viridis_c()
ggplot(wlist$founder_asym_growth_both_NAM_HAW, aes(x = nu1, y = nu2, color = AIC)) + geom_point() + theme_bw() + scale_color_viridis_c()

#==========call a function to make a parameter input file using these parms=============
nuA <- c(1e-2, 100) # ancient pop size
nu1 <- c(1e-3, 100) # final pop size, pop 1, exp growth
nu2 <- c(1e-6, 100) # final pop size, pop 2, exp growth
K1 <- c(1e-2, 100) # Carrying capacity of pop 1
K2 <- c(1e-5, 100) # Carrying capacity of pop 2
r1 <- c(0, 6) # Logistic growth rate of pop 2
r2 <- c(0, 6) # Logistic growth rate of pop 2
Ti <- c(1e-5, 20) # Time to (single) split
T1 <- c(1e-5, 10)  # Time of first epoch
T2 <- c(1e-5, 10) # Time of second epoch to present
s <- c(1e-5, .5) # fraction of nuA that goes to pop 2
m12 <- c(0, 40) # Migration rate from 2 to 1
m21 <- c(0, 40) # Migration rate from 1 to 2
m <- c(0, 40) # Symmetric migration rate
f <- c(1e-5, 1) # Fraction of updated population 2 to be derived from population 1 (admixture)
bounds <- list(nuA = nuA, nu1 = nu1, nu2 = nu2, K1 = K1, K2 = K2, r1 = r1, r2 = r2, 
               Ti = Ti, T1 = T1, T2 = T2, s = s, m12 = m12, m21 = m21, m = m, f = f)

# make.parm.file.from.best(wlist, mplist, bounds, 100, 100, "False", 1, "[15,15]", "dadi/parmfiles/dadi_HGR_3RD_pass_parms.txt")
make.parm.file.from.weighted.ave(wlist, mplist, bounds, 50, 60, "False", 2, "[100,10]", "dadi/parmfiles/NH_r2_adjusted_portik.txt")
