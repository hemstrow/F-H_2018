library(snpR); library(ggplot2); library(reshape2)
##########################
#get data
genos <- readRDS("flt_genotypes_good_inds.RDS")


##########################
#check error rates
dups <- cbind(genos[,1:3], genos[,grepl("HAW_9", colnames(genos))])
dups <- dups[dups$HAW_9 != "NN" & dups$HAW_9.1 != "NN",]

dups$eq <- ifelse(dups$HAW_9 == dups$HAW_9.1, 0, 1)
sum(dups$eq)/nrow(dups) #0.019262
##########################
#get the coancestry matrix
source("RAFM.r")

#format into RAFM format
pops <- table(substr(colnames(genos[,4:ncol(genos)]), 1, 3))
l <- list(c(names(pops)), as.numeric(pops))
Rgenos <- format_snps(genos, 3, 8, n_samp = 100, pop = l)


#for demographics, do guam rota NA Aus and HAW. Could also do just the island in HAW or just GUA/ROT.

#get the coancestries
RAFMall <- do.all(Rgenos, 100, 50, 2)

#get a matrix of the coancestry medians within and between each pop
comat <- matrix(NA, 13, 13)
for(i in 1:13){
  for(j in 1:13){
    comat[i,j] <- summary(RAFMall$theta[i,j,])[3]
  }
}

comat <- as.data.frame(comat)
colnames(comat) <- 1:13
comat$pop <- 1:13

mcomat <- melt(comat, id.vars = "pop")

#plot
ggplot(mcomat, aes(pop, variable, fill = value)) + geom_tile() + theme_bw() +
  scale_fill_gradient2(low = "white", high = "red") + geom_text(aes(label = round(value, 3)))
