# remotes::install_github("hemstrow/snpR", ref = "dev")
library(snpR)

x <- read.table("../F-H_2018/data/dadi_inputs/dadi_10kgap_snps.txt", header = T, stringsAsFactors = F)

# the make_SFS function arguments are: input data file, populations, projections (the number of
# sample gene copies (2 per individual!) to make the sfs from, lower numbers will keep more snps, maximizing snps
# kept over individuals kept tends to produce better results as long as you don't get TOO small),
# and then a logical for if the sfs should be folded or not. We can't fold the sfs if we are going to
# use the directionality index

# NAM, HAW
sf <- make_SFS(x, c("NAM", "HAW"), c(15, 15), F)
plot_sfs(sf)
NH <- calc_directionality(sf)

# HAW, GUA
sf <- make_SFS(x, c("HAW", "GUA"), c(10, 10), F)
plot_sfs(sf)
HG <- calc_directionality(sf)

# GUA, ROT
sf <- make_SFS(x, c("GUA", "ROT"), c(15, 15), F)
plot_sfs(sf)
GR <- calc_directionality(sf)

# HAW, QLD
sf <- make_SFS(x, c("HAW", "QLD"), c(10, 10), F)
plot_sfs(sf)
HQ <- calc_directionality(sf)

# HAW, NOR
sf <- make_SFS(x, c("HAW", "NOR"), c(10, 10), F)
plot_sfs(sf)
HN <- calc_directionality(sf)

sf <- make_SFS(x, c("NAM", "GUA"), c(15, 15))
plot_sfs(sf)

# NAM, QLD
sf <- make_SFS(x, c("NAM", "QLD"), c(15, 15), F)
plot_sfs(sf)
calc_directionality(sf)

# make a table of the comparisons of interest
dirs <- as.data.frame(t(combn(c("NAM", "HAW", "GUA", "QLD"), 2)), stringsAsFactors = F)
colnames(dirs)[1:2] <- c("pop1", "pop2")
dirs$di <- 0
dirs$direction <- ""
proj.sizes <- c(100, 10, 15, 20)
names(proj.sizes) <- c("NAM", "HAW", "GUA", "QLD")
for(i in 1:nrow(dirs)){
  print(dirs[i,1:2])
  sf <- make_SFS(x, c(dirs[i,1], dirs[i,2]), proj.sizes[match(dirs[i,1:2], names(proj.sizes))])
  di <- calc_directionality(sf)
  dirs$di[i] <- di
  dirs$direction[i] <- attr(di, "direction")
}


#==========plots for pres/paper================
pops <- colnames(x)[4:16]
pop.sample.sizes <- c(5, 30, 15, 100, 10, 20, 5, 5, 50, 20, 5, 10, 5)
pops <- pops[-which(pop.sample.sizes < 15)]
pop.sample.sizes <- pop.sample.sizes[-which(pop.sample.sizes < 15)]
out <- vector("list", length = (length(pops) * length(pops) - 1)/2)
iter <- 1
for(i in 1:(length(pops) - 1)){
  comps <- (i + 1):length(pops)
  for(j in comps){
    # generate SFS
    tsamps <- 10
    tsf <- make_SFS(x, c(pops[i], pops[j]), c(tsamps, tsamps))
    
    # generate a fuller SFS for the directionality
    dirtsf <- make_SFS(x, c(pops[i], pops[j]), c(pop.sample.sizes[i], pop.sample.sizes[j]))
    dir <- calc_directionality(dirtsf)
    dirN <- sum(dirtsf[-1,-1]) # no fixed sites
    
    # record data and melt
    ## scaled
    ntsf <- tsf/(sum(tsf)-tsf[1,1])
    colnames(ntsf) <- 0:(ncol(ntsf) - 1)
    ntsf <- as.data.frame(ntsf)
    ntsf[1,1] <- NA
    ntsf$p1 <- 0:(nrow(ntsf) - 1)
    ntsf <- reshape2::melt(ntsf, id.vars = "p1")
    colnames(ntsf)[2:3] <- c("p2", "Percent.Sites")
    
    ## not scaled
    colnames(tsf) <- 0:(ncol(tsf) - 1)
    tsf <- as.data.frame(tsf)
    tsf[1,1] <- NA
    N <- sum(tsf, na.rm = T)
    tsf$p1 <- 0:(nrow(tsf) - 1)
    tsf <- reshape2::melt(tsf, id.vars = "p1")
    colnames(tsf)[2:3] <- c("p2", "Num.Sites")
    
    tsf <- merge(tsf, ntsf, by = c("p1", "p2"))
    
    # save
    tsf$dir <- dir
    tsf$N <- N
    tsf$dirN <- dirN
    tsf$pop_1 <- pops[j]
    tsf$pop_2 <- pops[i]
    out[[iter]] <- tsf
    iter <- iter + 1
  }
}

out <- dplyr::bind_rows(out)
out$p1 <- as.numeric(out$p1)
out$p2 <- as.numeric(out$p2)

# make a plot
library(ggplot2)
pdf("./plots/directionality_and_sfs.pdf", width = 10, height = 10)
ggplot(out, aes(x = p1, y = p2, color = log10(Num.Sites), fill = log10(Num.Sites))) +
  facet_grid(pop_1 ~ pop_2, switch = "both") + geom_tile() +
  geom_text(aes(label = round(dir, 3), x = 1.5, y = 2.5), color = "black") + 
  ggplot2::scale_color_viridis_c(na.value = "white", option = "inferno") +
  ggplot2::scale_fill_viridis_c(na.value = "white", option = "inferno") +
  theme_bw() + 
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  theme(strip.background = element_blank())
dev.off();dev.off();
