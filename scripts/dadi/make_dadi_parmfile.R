############3pop-models##############
nu1B <- c(0.01, 1e-4, .5) 
nu2B <- c(0.01, 1e-4, .5)
nu1f <- c(0.1, 1e-4, 1)
nu2f <- c(0.1, 1e-4, 1)
K2 <- c(0.1, 1e-4, 1)
K3 <- c(0.1, 1e-4, 1)
ts <- c(.1, .001, 2)
tp <- c(.1, .001, 2)
m21 <- c(.1, 0, 1)
m31 <- c(.1, 0, 1)
m32 <- c(.1, 0, 1)
m23 <- c(.1, 0, 1)
r2 <- c(1, .01, 4)
r3 <- c(1, .01, 4)

#which models?
mods <- c(rep("p3_cgrowth", 20),
          rep("p3_lgrowth_both", 20),
          rep("p3_lgrowth_p3", 20),
          rep("p3_lgrowth_p2", 20),
          rep("p3_lgb_p1tb_2f", 10))

#how many iterations?
maxiters <- c(rep(10, length(mods)))

#which pops?
pops <- c(rep(c(rep("[NAM,HAW,GUA]", 10), rep("[NAM,GUA,HAW]", 10)), 4),
          rep("[NAM,HAW,GUA]", 10))

#fold the spectra?
fs <- rep("False", length(mods))

#fold the p0?
fp <- rep(3, length(mods))

#what projection?
proj <- rep("[150,15,15]", length(mods))

#what optimizer?
optimizer <- rep("fmin", length(mods))

##set ip, lb, ub
ip <- character(length(mods))
lb <- ip
ub <- ip
for(i in 1:length(ip)){
  if(mods[i] == "cgrowth"){
    tip <- paste0(c(nu1B[1], nu2B[1], nu1f[1], nu2f[1], ts[1], tp[1], m21[1], m31[1], m32[1], m23[1]), collapse = ",")
    tlb <- paste0(c(nu1B[2], nu2B[2], nu1f[2], nu2f[2], ts[2], tp[2], m21[2], m31[2], m32[2], m23[2]), collapse = ",")
    tub <- paste0(c(nu1B[3], nu2B[3], nu1f[3], nu2f[3], ts[3], tp[3], m21[3], m31[3], m32[3], m23[3]), collapse = ",")
  }
  else if (mods[i] == "lgrowth_both" | mods[i] == "lgb_p1tb_2f"){
    tip <- paste0(c(nu1B[1], nu2B[1], K2[1], K3[1], ts[1], tp[1], m21[1], m31[1], m32[1], m23[1], r2[1], r3[1]), collapse = ",")
    tlb <- paste0(c(nu1B[2], nu2B[2], K2[2], K3[2], ts[2], tp[2], m21[2], m31[2], m32[2], m23[2], r2[2], r3[2]), collapse = ",")
    tub <- paste0(c(nu1B[3], nu2B[3], K2[3], K3[3], ts[3], tp[3], m21[3], m31[3], m32[3], m23[3], r2[3], r3[3]), collapse = ",")
  }
  else if (mods[i] == "lgrowth_p2"){
    tip <- paste0(c(nu1B[1], nu2B[1], K2[1], nu2f[1], ts[1], tp[1], m21[1], m31[1], m32[1], m23[1], r2[1]), collapse = ",")
    tlb <- paste0(c(nu1B[2], nu2B[2], K2[2], nu2f[2], ts[2], tp[2], m21[2], m31[2], m32[2], m23[2], r2[2]), collapse = ",")
    tub <- paste0(c(nu1B[3], nu2B[3], K2[3], nu2f[3], ts[3], tp[3], m21[3], m31[3], m32[3], m23[3], r2[3]), collapse = ",")
  }
  else if (mods[i] == "lgrowth_p3"){
    tip <- paste0(c(nu1B[1], nu2B[1], nu1f[1], K3[1], ts[1], tp[1], m21[1], m31[1], m32[1], m23[1], r3[1]), collapse = ",")
    tlb <- paste0(c(nu1B[2], nu2B[2], nu1f[2], K3[2], ts[2], tp[2], m21[2], m31[2], m32[2], m23[2], r3[2]), collapse = ",")
    tub <- paste0(c(nu1B[3], nu2B[3], nu1f[3], K3[3], ts[3], tp[3], m21[3], m31[3], m32[3], m23[3], r3[3]), collapse = ",")
  }
  ip[i] <- paste0("[", tip, "]")
  lb[i] <- paste0("[", tlb, "]")
  ub[i] <- paste0("[", tub, "]")
}


parms <- cbind(mods, maxiters, pops, fs, fp, ip, ub, lb, proj, optimizer)

#add to rota models
parms2 <- gsub("GUA", "ROT", parms)



##############2-pop models##############

#model parameter p0, lower bounds, and upper bounds.
nuA <- c(0.01, 1e-6, 100) # ancient pop size
nu1 <- c(0.01, 1e-6, 100) # final pop size, pop 1, exp growth
nu2 <- c(0.01, 1e-6, 100) # final pop size, pop 2, exp growth
K1 <- c(0.1, 1e-5, 100) # Carrying capacity of pop 1
K2 <- c(0.1, 1e-5, 100) # Carrying capacity of pop 2
r1 <- c(1, 0, 6) # Logistic growth rate of pop 2
r2 <- c(1, 0, 6) # Logistic growth rate of pop 2
Ti <- c(.1, .001, 20) # Time to (single) split
T1 <- c(.1, .001, 10)  # Time of first epoch
T2 <- c(.1, .001, 10) # Time of second epoch to present
s <- c(0.1, 1e-5, .5) # fraction of nuA that goes to pop 2
m12 <- c(.1, 0, 40) # Migration rate from 2 to 1
m21 <- c(.1, 0, 40) # Migration rate from 1 to 2
m <- c(.1, 0, 40) # Symmetric migration rate
f <- c(0.1, 1e-5, 1) # Fraction of updated population 2 to be derived from population 1 (admixture)


#how many iters?
iters <- 50

#which models and how many reps?
reps_per_perm <- 100
mods <- c(
  # vicariance models
  vic_no_mig = T, 
  vic_anc_asym_mig = T,  
  vic_sec_contact_asym_mig = T,
  vic_no_mig_admix_early = T,  
  vic_no_mig_admix_late = T,  
  vic_two_epoch_admix = T, 
  
  # founder exponential models
  founder_nomig_growth_pop_2 = T,  
  founder_sym_growth_pop_2 = T, 
  founder_asym_growth_pop_2 = T, 
  founder_nomig_admix_early_growth_pop_2 = T, 
  founder_nomig_admix_late_growth_pop_2 = T,  
  founder_nomig_admix_two_epoch_growth_pop_2 = T, 
  
  founder_nomig_growth_pop_1 = T,  
  founder_sym_growth_pop_1 = T, 
  founder_asym_growth_pop_1 = T, 
  founder_nomig_admix_early_growth_pop_1 = T, 
  founder_nomig_admix_late_growth_pop_1 = T,  
  founder_nomig_admix_two_epoch_growth_pop_1 = T,
  
  founder_nomig_growth_both = T,  
  founder_sym_growth_both = T, 
  founder_asym_growth_both = T, 
  founder_nomig_admix_early_growth_both = T, 
  founder_nomig_admix_late_growth_both = T,  
  founder_nomig_admix_two_epoch_growth_both = T, 
  
  # founder logistic models
  founder_nomig_logistic_growth_pop_2 = T,  
  founder_sym_logistic_growth_pop_2 = T, 
  founder_asym_logistic_growth_pop_2 = T,  
  founder_nomig_admix_early_logistic_growth_pop_2 = T, 
  founder_nomig_admix_late_logistic_growth_pop_2 = T,  
  founder_nomig_admix_two_epoch_logistic_growth_pop_2 = T,
  
  founder_nomig_logistic_growth_pop_1 = T,  
  founder_sym_logistic_growth_pop_1 = T, 
  founder_asym_logistic_growth_pop_1 = T,  
  founder_nomig_admix_early_logistic_growth_pop_1 = T, 
  founder_nomig_admix_late_logistic_growth_pop_1 = T,  
  founder_nomig_admix_two_epoch_logistic_growth_pop_1 = T,
  
  founder_nomig_logistic_growth_both = T,  
  founder_sym_logistic_growth_both = T, 
  founder_asym_logistic_growth_both = T,  
  founder_nomig_admix_early_logistic_growth_both = T, 
  founder_nomig_admix_late_logistic_growth_both = T,  
  founder_nomig_admix_two_epoch_logistic_growth_both = T
)

#which pops to use?
pops <- c("[NAM,HAW]")

#is the spectra polarized?
fs <- "False"

#fold initial params?
fp <- 3

#what projection?
proj <- "[100,10]"

#which optimizer?
optim <- "fmin"

#outfile?
ofile <- "dadi/parmfiles/NH_r1_adjusted_portik.txt"
append.ofile <- F # should this be appened to an existing outfile?

##########################################
#make the parm df
#models
mods <- names(mods)[mods]
rpp <- reps_per_perm
rppm <- rpp*length(pops)
mods <- sort(rep(mods, rppm))

#how many iterations?
maxiters <- rep(iters, length(mods))

#which pops?
rppp <- rpp*length(unique(mods))
pops <- rep(pops, rppp)

#fold the spectra?
fs <- rep(fs, length(mods))

#fold initial params?
fp <- rep(fp, length(mods))

#what projection?
proj <- rep(proj, length(mods))

#what optimizer?
optim <- rep(optim, length(mods))


##set ip, lb, ub
ip <- character(length(mods))
lb <- ip
ub <- ip
for(i in 1:length(ip)){
  
  # vicariance models
  if(mods[i] == "vic_no_mig"){
    tparms <- cbind(nuA, Ti, s)
  }
  else if(mods[i] == "vic_anc_asym_mig"){
    tparms <- cbind(nuA, m12, m21, T1, T2, s)
  }
  else if(mods[i] == "vic_sec_contact_asym_mig"){
    tparms <- cbind(nuA, m12, m21, T1, T2, s)
  }
  else if(mods[i] == "vic_no_mig_admix_early"){
    tparms <- cbind(nuA, Ti, s, f)
  }
  else if(mods[i] == "vic_no_mig_admix_late"){
    tparms <- cbind(nuA, Ti, s, f)
  }
  else if(mods[i] == "vic_two_epoch_admix"){
    tparms <- cbind(nuA, T1, T2, s, f)
  }
  
  
  # exponential growth models
  else if(mods[i] == "founder_nomig_growth_pop_2"){
    tparms <- cbind(nuA, nu2, Ti, s)
  }
  else if(mods[i] == "founder_sym_growth_pop_2"){
    tparms <- cbind(nuA, nu2, m, Ti, s)
  }
  else if(mods[i] == "founder_asym_growth_pop_2"){
    tparms <- cbind(nuA,  nu2, m12, m21, Ti, s)
  }
  else if(mods[i] == "founder_nomig_admix_early_growth_pop_2"){
    tparms <- cbind(nuA, nu2, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_late_growth_pop_2"){
    tparms <- cbind(nuA, nu2, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_two_epoch_growth_pop_2"){
    tparms <- cbind(nuA, nu2, T1, T2, s, f)
  }
  
  else if(mods[i] == "founder_nomig_growth_pop_1"){
    tparms <- cbind(nuA, nu1, Ti, s)
  }
  else if(mods[i] == "founder_sym_growth_pop_1"){
    tparms <- cbind(nuA, nu1, m, Ti, s)
  }
  else if(mods[i] == "founder_asym_growth_pop_1"){
    tparms <- cbind(nuA, nu1, m12, m21, Ti, s)
  }
  else if(mods[i] == "founder_nomig_admix_early_growth_pop_1"){
    tparms <- cbind(nuA, nu1, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_late_growth_pop_1"){
    tparms <- cbind(nuA, nu1, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_two_epoch_growth_pop_1"){
    tparms <- cbind(nuA, nu1, T1, T2, s, f)
  }
  
  else if(mods[i] == "founder_nomig_growth_both"){
    tparms <- cbind(nuA, nu1, nu2, Ti, s)
  }
  else if(mods[i] == "founder_sym_growth_both"){
    tparms <- cbind(nuA, nu1, nu2, m, Ti, s)
  }
  else if(mods[i] == "founder_asym_growth_both"){
    tparms <- cbind(nuA, nu1, nu2, m12, m21, Ti, s)
  }
  else if(mods[i] == "founder_nomig_admix_early_growth_both"){
    tparms <- cbind(nuA, nu1, nu2, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_late_growth_both"){
    tparms <- cbind(nuA, nu1, nu2, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_two_epoch_growth_both"){
    tparms <- cbind(nuA, nu1, nu2, T1, T2, s, f)
  }
  
  
  # logistic growth models
  else if(mods[i] == "founder_nomig_logistic_pop_2"){
    tparms <- cbind(nuA, K2, r2, Ti, s)
  }
  else if(mods[i] == "founder_sym_logistic_pop_2"){
    tparms <- cbind(nuA, K2, r2, m, Ti, s)
  }
  else if(mods[i] == "founder_asym_logistic_pop_2"){
    tparms <- cbind(nuA, K2, r2, m12, m21, Ti, s)
  }
  else if(mods[i] == "founder_nomig_admix_early_logistic_pop_2"){
    tparms <- cbind(nuA, K2, r2, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_late_logistic_pop_2"){
    tparms <- cbind(nuA, K2, r2, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_two_epoch_logistic_pop_2"){
    tparms <- cbind(nuA, K2, r2, T1, T2, s, f)
  }
  
  else if(mods[i] == "founder_nomig_logistic_pop_1"){
    tparms <- cbind(nuA, K1, r1, Ti, s)
  }
  else if(mods[i] == "founder_sym_logistic_pop_1"){
    tparms <- cbind(nuA, K1, r1, m, Ti, s)
  }
  else if(mods[i] == "founder_asym_logistic_pop_1"){
    tparms <- cbind(nuA, K1, r1, m12, m21, Ti, s)
  }
  else if(mods[i] == "founder_nomig_admix_early_logistic_pop_1"){
    tparms <- cbind(nuA, K1, r1, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_late_logistic_pop_1"){
    tparms <- cbind(nuA, K1, r1, Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_two_epoch_logistic_pop_1"){
    tparms <- cbind(nuA, K1, r1, T1, T2, s, f)
  }
  
  else if(mods[i] == "founder_nomig_logistic_both"){
    tparms <- cbind(nuA, K1, K2, r1, r2Ti, s)
  }
  else if(mods[i] == "founder_sym_logistic_both"){
    tparms <- cbind(nuA, K1, K2, r1, r2m, Ti, s)
  }
  else if(mods[i] == "founder_asym_logistic_both"){
    tparms <- cbind(nuA, K1, K2, r1, r2m12, m21, Ti, s)
  }
  else if(mods[i] == "founder_nomig_admix_early_logistic_both"){
    tparms <- cbind(nuA, K1, K2, r1, r2Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_late_logistic_both"){
    tparms <- cbind(nuA, K1, K2, r1, r2Ti, s, f)
  }
  else if(mods[i] == "founder_nomig_admix_two_epoch_logistic_both"){
    tparms <- cbind(nuA, K1, K2, r1, r2T1, T2, s, f)
  }
  ip[i] <- paste0("[", paste0(tparms[1,], collapse = ","), "]")
  lb[i] <- paste0("[", paste0(tparms[2,], collapse = ","), "]")
  ub[i] <- paste0("[", paste0(tparms[3,], collapse = ","), "]")
}


#combine
parms <- cbind(mods, maxiters, pops, fs, fp, ip, ub, lb, proj, optim)

#write
write.table(parms, ofile, quote = F, sep = " ", row.names = F, col.names = F, append = append.ofile)


