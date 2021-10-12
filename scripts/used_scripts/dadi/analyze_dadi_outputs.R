library(ggplot2); library(dplyr); source("scripts/dadi/make_dadi_parmfile_from_results.R"); source("scripts/interpret_dadi_units.R")
#==========import and prep data==========

res <- readLines("data/dadi_inputs/cat_NH_unfolded_r3.txt")
#res <- c(res, readLines("data/dadi_inputs/cat_NH_best_r2.txt"))

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
  
  founder_nomig_admix_two_epoch_logistic_both = c("nuA", "K1", "K2", "r1", "r2", "T1", "T2", "s", "f"),
  
  # historic growth models
  founder_asym_hist_igrowth_p2 = c("nuA", "nuG", "nu2F", "m12", "m21", "Tg", "Ts", "Tg2", "s"),
  
  founder_asym_hist_3epoch_exp_growth_p1 = c("nuA", "nuG", "nu1F", "nuG2", "nu2F", "m12", "m21", "Tg", "Tg2", "Ts", "Tg3", "s")
  
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
# number of SNPs in dataset: 10000
# number of SNPs in total set: 302446
# number of sequenced sites: 1373747
mu <- 8.4e-9
g <- .3
ratio <- 10929/302446 # ratio of included snps
L <-  1373747*ratio # approx number of considered bases.
ilist <- interpret_units(wlist, mu, g, L = L)


#ggplot(ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, aes(y = log10(Tg3 + Ts), x = log10(s), color = log10(AIC))) + geom_point() + theme_bw() + scale_color_viridis_c()
#ggplot(ilist$founder_asym_hist_igrowth_p2_NAM_HAW, aes(y = log10(Tg2 + Ts), x = log10(s), color = log10(AIC))) + geom_point() + theme_bw() + scale_color_viridis_c()


# exp <- list(asym_both <- ilist$founder_asym_growth_both_NAM_HAW,
#             asym_NA <- ilist$founder_asym_growth_pop_1_NAM_HAW)
# saveRDS(exp, "dadi_run3_best_models.RDS")

#==========call a function to make a parameter input file using these parms=============
bounds <- list(nuA = c(1e-2, 100), # ancient pop size
               nu1 = c(1e-3, 100), # final pop size, pop 1, exp growth
               nu2 = c(1e-6, 100), # final pop size, pop 2, exp growth
               nuG = c(1e-3, 1000), # final pop size, historic growth period 1
               nuG2 = c(1e-3, 1000), # final pop size, historic growth period 2
               nu2F = c(1e-6, 100), # final size, pop 2, hist growth-instant growth
               nu1F = c(1e-6, 100), # final size, pop 1, hist growth
               K1 = c(1e-2, 100), # Carrying capacity of pop 1
               K2 = c(1e-5, 100), # Carrying capacity of pop 2
               r1 = c(0, 6), # Logistic growth rate of pop 2
               r2 = c(0, 6), # Logistic growth rate of pop 2
               Ti = c(1e-5, 20), # Time to (single) split
               T1 = c(1e-5, 10), # Time of first epoch
               T2 = c(1e-5, 10), # Time of second epoch to present
               Tg = c(1e-5, 20), # Time between historic growth and split/phase 2 growth
               Tg2 = c(1e-10, 20), # Time between phase two growth and split or Time between split and present
               Ts = c(1e-10, 20), # Time between split and instant p2 growth.
               Tg3 = c(1e-10, 20), # Time between instant p2 growth and present.
               s = c(1e-5, .5), # fraction of nuA that goes to pop 2
               m12 = c(0, 40), # Migration rate from 2 to 1
               m21 = c(0, 40), # Migration rate from 1 to 2
               m = c(0, 40), # Symmetric migration rate
               f = c(1e-5, 1) # Fraction of updated population 2 to be derived from population 1 (admixture)
)

#make.parm.file.from.best(wlist, mplist, bounds, 100, 100, "False", 1, "[100,10]", "dadi/parmfiles/NH_best_r4.txt")
make.parm.file.from.weighted.ave(wlist,
                                 mplist, bounds, iters = 100, reps_per_perm = 100, "True", 1, "[100,10]", "dadi/parmfiles/NH_unfolded_r4.txt")
