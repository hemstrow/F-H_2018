mods <- c(rep("cgrowth", 20),
          rep("lgrowth_both", 20),
          rep("lgrowth_p3", 20),
          rep("lgrowth_p2", 20),
          rep("lgb_p1tb_2f", 10))
maxiters <- c(rep(10, length(mods)))
pops <- c(rep(c(rep("[NAM,HAW,GUA]", 10), rep("[NAM,GUA,HAW]", 10)), 4),
          rep("[NAM,HAW,GUA]", 10))
fs <- rep("True", length(mods))
fp <- rep(3, length(mods))
ip <- c(rep("[0.1,0.1,0.5,0.5,5,5,.1,.1,.1,.1]", 20),
        rep("[0.1,0.1,1,1,5,5,.1,.1,.1,.1,2,2]", 20),
        rep("[0.1,0.1,0.5,1,5,5,.1,.1,.1,.1,2]", 40),
        rep("[0.1,0.1,1,1,5,5,.1,.1,.1,.1,2,2]", 10))
ub <- c(rep("[1,1,5,5,20,20,10,10,10,10]", 20),
        rep("[1,1,10,10,20,20,10,10,10,10,4,4]", 20),
        rep("[1,1,5,10,20,20,10,10,10,10,4]", 40),
        rep("[1,1,10,10,20,20,10,10,10,10,4,4]", 10))
lb <- c(rep("[.001,.001,.001,.001,.01,.01,.01,.01,.01,.01]", 20),
        rep("[1e-3,1e-3,1e-3,1e-3,0.01,0.01,0.01,0.01,0.01,0.01,1,1]", 20),
        rep("[1e-3,1e-3,1e-3,1e-3,0.01,0.01,0.01,0.01,0.01,0.01,1]", 40),
        rep("[1e-3,1e-3,1e-3,1e-3,0.01,0.01,0.01,0.01,0.01,0.01,1,1]", 10))
proj <- rep("[150,15,15]", length(mods))
optimizer <- rep("fmin", length(mods))
parms <- cbind(mods, maxiters, pops, fs, fp, ip, ub, lb, proj, optimizer)
write.table(parms, "DADI/parmfiles/1st_folded_optim_snps.txt", quote = F, sep = " ", row.names = F, col.names = F)
