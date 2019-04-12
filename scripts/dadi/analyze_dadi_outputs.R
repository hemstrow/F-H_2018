library(ggplot2); library(dplyr)
#==========import and prep data==========

res <- readLines("data/dadi_inputs/cat_NAH_portik_1st_out.txt")

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
  
  founder_asym_growth_pop_2 = c("nuA", " nu2", "m12", "m21", "Ti", "s"),
  
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
  wlist[[i]] <- po
}


#==========analysis========
#overall between models.
rdf$AIC <- as.numeric(rdf$AIC)
rdf$ll <- as.numeric(rdf$ll)
ggplot(rdf, aes(x = AIC)) + geom_histogram() + facet_grid(model ~ pops) + theme_bw()

#best replicate per model
best.reps <- rdf %>% group_by(model) %>% group_by(pops, add = TRUE) %>% top_n(-1, AIC)
best.reps <- arrange(best.reps, pops, model)
best.reps[,-c(3,4,5,7)] #where are we pushing parameter bounds?
ggplot(best.reps, aes(y = AIC, x = model, color = pops)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
best.reps[order(best.reps$AIC),-c(3,4,7)]

#these best.reps are the starting points for next run


# plots for the parameters
ggplot(wlist$founder_asym_HAW_NAM, aes(x = s, y = nuA, color = AIC)) + geom_point() + theme_bw() + scale_color_viridis_c()






#==========call a function to make a parameter input file using these parms=============
nu1B <- c(0.001, 200)
nu2B <- c(0.001, 200)
nu1f <- c(0.01, 200)
nu2f <- c(0.01, 200)
K1 <- c(0.01, 20)
K2 <- c(0.01, 20)
ts <- c(.01, 10)
tp <- c(.01, 10)
m12 <- c(.01, 40)
m21 <- c(.01, 40)
r1 <- c(1, 6)
r2 <- c(1, 6)
bounds <- list(nu1B = nu1B, nu2B = nu2B, nu1f = nu1f, nu2f = nu2f,
               K1 = K1, K2 = K2, ts = ts, tp = tp, m12 = m12, m21 = m21, r1 = r1, r2 = r2)

make.parm.file.from.best(best.reps, mplist, bounds, 100, 100, "False", 1, "[15,15]", "dadi/parmfiles/dadi_HGR_3RD_pass_parms.txt")
