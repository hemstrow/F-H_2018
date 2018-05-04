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
mods <- c(rep("cgrowth", 20),
          rep("lgrowth_both", 20),
          rep("lgrowth_p3", 20),
          rep("lgrowth_p2", 20),
          rep("lgb_p1tb_2f", 10))

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



##############2-pop models for GUA/ROT and HAW only!##############

#model parameter p0, lower bounds, and upper bounds.
nu1B <- c(0.01, 1e-6, .5)
nu2B <- c(0.01, 1e-6, .5)
nu1f <- c(0.1, 1e-5, 1)
nu2f <- c(0.1, 1e-5, 1)
K1 <- c(0.1, 1e-5, 1)
K2 <- c(0.1, 1e-5, 1)
ts <- c(.1, .001, 4)
tp <- c(.1, .001, 4)
m12 <- c(.1, 0, 5)
m21 <- c(.1, 0, 5)
r1 <- c(1, .01, 4)
r2 <- c(1, .01, 4)

#how many iters?
iters <- 20

#which models and how many reps?
reps_per_perm <- 50
mods <- c("p2_cgrowth", "p2_lgrowth_both", "p2_lgrowth_1", "p2_lgrowth_2")

#which pops to use?
pops <- c("[HAW,GUA]", "[HAW,ROT]", "[GUA,HAW]", "[ROT,HAW]")

#fold the spectra?
fs <- "False"

#fold initial params?
fp <- 3

#what projection?
proj <- "[15,15]"

#which optimizer?
optim <- "fmin"

#outfile?
ofile <- "dadi/parmfiles/1st_passHGR_2d_parms.txt"

##########################################
#make the parm df
#models
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
  if(mods[i] == "p2_cgrowth"){
    tparms <- cbind(nu1B, nu2B, nu1f, nu2f, ts, tp, m12, m21)
  }
  else if (mods[i] == "p2_lgrowth_both"){
    tparms <- cbind(nu1B, nu2B, K1, K2, ts, tp, m12, m21, r1, r2)
  }
  else if (mods[i] == "p2_lgrowth_1"){
    tparms <- cbind(nu1B, nu2B, K1, nu2f, ts, tp, m12, m21, r1)
  }
  else if (mods[i] == "p2_lgrowth_2"){
    tparms <- cbind(nu1B, nu2B, nu1f, K2, ts, tp, m12, m21, r2)
  }
  ip[i] <- paste0("[", paste0(tparms[1,], collapse = ","), "]")
  lb[i] <- paste0("[", paste0(tparms[2,], collapse = ","), "]")
  ub[i] <- paste0("[", paste0(tparms[3,], collapse = ","), "]")
}


#combine
parms <- cbind(mods, maxiters, pops, fs, fp, ip, ub, lb, proj, optim)

#write
write.table(parms, ofile, quote = F, sep = " ", row.names = F, col.names = F)


