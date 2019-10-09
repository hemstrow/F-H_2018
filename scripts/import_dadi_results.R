import_dadi_results <- function(res, mu, g, L){
  
  res_out <- character(0)
  for(i in 1:length(res)){
    res_out <- c(res_out, readLines(res[i]))
  }
  res <- res_out
  rm(res_out)
  
  
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
  
  # interpret units
  ilist <- interpret_units(wlist, mu, g, L = L)
  
  # return
  return(list(rdf = rdf, wlist = wlist, ilist = ilist))
}