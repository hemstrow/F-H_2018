#makes a parm file for dadi where the initial values are the best given by the previous dadi pass.
#inputs:
#    best.reps: data.table or data.frame containing info on the best reps for each model/pop comparison, as given in analyze_dadi_outputs.R
#    mplist: List of parameters per model, as given in analyze_dadi_outputs.R.
#    bounds: Named list of bounds of upper and lower parameter bounds, as given in analyze_dadi_outputs.R.
#    iters: Max number of iterations for dadi to use.
#    reps_per_perm: Number of times to run each parameter set.
#    fs: is the spectra polarized? For unfolded spectra, this should be "False".
#    proj: Projection to use, format "[15, 15]".
#    ofile: Output filehandle.
#    optim: Optimization formula to use.
#    exclude: Numeric vector stipulating which rows from best.reps should be skipped when making this parmfile.

make.parm.file.from.best <- function(best.reps, mplist, bounds, iters, 
                                     reps_per_perm, fs, fp, proj, ofile, optim = "fmin",exclude = NA){
  if(!is.na(exclude)){
    best.reps <- best.reps[-exclude,]
  }
  
  #initialize vectors to store info
  rppm <- reps_per_perm*nrow(best.reps)
  mods <- rep(best.reps$model, each = reps_per_perm)
  pops <- rep(paste0("[", gsub(" ", ",", best.reps$pops), "]"), each = reps_per_perm)
  ip <- character(length(mods))
  ub <- ip
  lb <- ip
  counter <- 1
  
  #fill vectors
  for(i in 1:nrow(best.reps)){
    #add starting values
    ip[counter:(counter + reps_per_perm - 1)] <-
      paste0("[", gsub(" ", ",", best.reps$mnum[i]), "]")
    
    #grab upper and lower bounds
    tb <- unlist(mplist[names(mplist) == best.reps$model[i]])
    tb <- bounds[names(bounds) %in% tb]
    tb <- as.data.frame(tb)
    lb[counter:(counter + reps_per_perm - 1)] <- 
      paste0("[", paste0(tb[1,], collapse = ","), "]")
    ub[counter:(counter + reps_per_perm - 1)] <- 
      paste0("[", paste0(tb[2,], collapse = ","), "]")
    
    #update counter
    counter <- counter + reps_per_perm
  }
  
  #bind for printing
  parms <- cbind(mods, iters, pops, fs, fp, ip, ub, lb, proj, optim)
  
  write.table(parms, ofile, quote = F, sep = " ", row.names = F, col.names = F)
}


#makes a parm file for dadi where the initial values are AIC weighted values given by the previous dadi pass.
#inputs:
#    wlist: list containing info on reps for each model/pop comparison, as given in analyze_dadi_outputs.R (wlist)
#    mplist: List of parameters per model, as given in analyze_dadi_outputs.R.
#    bounds: Named list of bounds of upper and lower parameter bounds, as given in analyze_dadi_outputs.R.
#    iters: Max number of iterations for dadi to use.
#    reps_per_perm: Number of times to run each parameter set.
#    fs: is the spectra polarized? For unfolded spectra, this should be "False".
#    fp: factor by which to permute input data. 1 is a little, 3 is a lot.
#    proj: Projection to use, format "[15, 15]".
#    ofile: Output filehandle.
#    optim: Optimization formula to use.
#    exclude: Numeric vector stipulating which models from wlist should be skipped when making this parmfile.
#    quant: AIC quantile above which to discard data when finding weighted parameter values.
make.parm.file.from.weighted.ave <- function(wlist, mplist, bounds, iters,
                                       reps_per_perm, fs, fp, proj, ofile, optim = "fmin", exclude = NA, quant = FALSE){
  if(!is.na(exclude)){
    wlist <- wlist[-exclude,]
  }

  #initialize vectors to store info
  rppm <- reps_per_perm*length(wlist)
  mods <- rep(names(wlist), each = reps_per_perm)
  pops <- character(length(mods))
  ip <- pops
  ub <- ip
  lb <- ip
  
  
  # find the weighted average parameters for each model.
  fit.parms <- vector(mode = "list", length = length(wlist))
  names(fit.parms) <- names(wlist)
  for(i in 1:length(wlist)){
    tmod <- wlist[[i]]
    AIC <- tmod$AIC
    if(quant){
      keep <- which(AIC <= quantile(AIC, quant))
      AIC <- AIC[keep]
      tmod <- tmod[keep,]
    }
    # weights are the inverse of AIC
    wt <- 1/AIC
    tmod <- tmod[,-c(1:3)]
    fit.parms[[i]] <- round(colSums(tmod*wt)/sum(wt), 7)
  }
  
  counter <- 1
  #fill vectors
  for(i in 1:length(wlist)){
    #add starting values
    ip[counter:(counter + reps_per_perm - 1)] <-
      paste0("[", paste0(fit.parms[[i]], collapse = ","), "]")
    
    #grab upper and lower bounds
    tb <- unlist(mplist[names(mplist) == rdf$model[i]])
    tb <- bounds[names(bounds) %in% tb]
    tb <- as.data.frame(tb)
    lb[counter:(counter + reps_per_perm - 1)] <- 
      paste0("[", paste0(tb[1,], collapse = ","), "]")
    ub[counter:(counter + reps_per_perm - 1)] <- 
      paste0("[", paste0(tb[2,], collapse = ","), "]")
    
    # add pops
    pops[counter:(counter + reps_per_perm - 1)] <- paste0("[", gsub(" ", ",", wlist[[i]]$pops[1]), "]")
    
    #update counter
    counter <- counter + reps_per_perm
  }
  
  #bind for printing
  parms <- cbind(mods, iters, pops, fs, fp, ip, ub, lb, proj, optim)
  
  write.table(parms, ofile, quote = F, sep = " ", row.names = F, col.names = F)
}