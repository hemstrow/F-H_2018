library(RColorBrewer); library(ape); library(ggplot2); library(gridExtra); library(gridGraphics); library(data.table); 
library(snpR); library(ggtree)

# This script makes all of the R plots for the publication. Plots are as follows, 
# asterisks imply that the plot was not made with R:
# 
# Figure 1: Map of sampling locations. *




facet.order <- c("ENA", "WNA", "MAU", "OAH", "GUA", "ROT", "SAI", "SAM", "FIJ", "NCA", "NOR", "QLD", "NSW", "VIC", "NZL")
full.label <- c("Eastern North America (ENA, 45)",
                "Western North America (WNA, 40)",
                "Maui (MAU, 8)",
                "Oahu (OAH, 4)",
                "Guam (GUA, 24)",
                "Rota (ROT, 20)",
                "Saipan (SAI, 4)",
                "Samoa (31)",
                "Fiji (5)",
                "New Caledonia (NCA, 18)",
                "Norfolk Island (NOR, 16)",
                "Queensland (QLD, 44)",
                "New South Wales (NSW, 6)",
                "Victoria (VIC, 4)",
                "New Zealand (NZL, 16)")
full.label.tab <- data.frame(abrv = facet.order, full = full.label)

pal <- colorRampPalette(brewer.pal(9, "Set1"))(length(facet.order)) #palette to use
color.guide <- data.frame(pop = facet.order, color = pal, stringsAsFactors = F)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}

#=========Neighbor-joining tree========
# grab data
setwd("Raw_data/")
bamlist <- read.table("gbamlist.txt")
poorlist <- read.table("poorlist.txt")
bamlist$poor <- bamlist$V1 %in% poorlist$V1
bamlist$ord <- 1:nrow(bamlist)
plate1 <- read.csv("plate1.csv", header = T)
plate2 <- read.csv("plate2_MFMB1.csv", header = T, skip = 1)
plate3 <- read.csv("plate3_MFMB2.csv", header = T, skip = 1)

DNA274 <- "ACAGTG"
DNA275 <- "CAGATC"

plate2 <- rbind(plate2[,1:3], setNames(plate2[,4:6], colnames(plate2)[1:3]))
plate3 <- rbind(plate3[,1:3], setNames(plate3[,4:6], colnames(plate3)[1:3]))


#add index 
barcodes <- read.table("RAD barcodes.txt", header = T)
plate2 <- merge(plate2, barcodes, by = "Well")
plate3 <- merge(plate3, barcodes, by = "Well")
plate1 <- merge(plate1, barcodes, by = "Well")


#get sample info
sampinfo <- data.frame(plate = substr(bamlist[,1], 8, 13), ind = substr(bamlist[,1], 20, 27), ord = bamlist$ord, poor = bamlist$poor, stringsAsFactors = FALSE)
sampinfo$plate[sampinfo$plate == "_split"] <- "plate1"



#add sample info to plates
plate2 <- merge(sampinfo[sampinfo$plate == DNA274,], plate2, by.x = "ind", by.y = "Index")
plate3 <- merge(sampinfo[sampinfo$plate == DNA275,], plate3, by.x = "ind", by.y = "Index")
plate1 <- merge(sampinfo[sampinfo$plate == "plate1",], plate1, by.x = "ind", by.y = "Index")

#combine plates in the correct order
combplates <- rbind(plate2, plate3, plate1)

#fix up some pop IDs
pops <- substr(combplates$Pop, 1, 3)
table(pops)
pops[pops == "Gua"] <- "GUA"
pops[pops == "New"] <- "NCA"
pops[pops == "Nor"] <- "NOR"
pops[pops == "NZ2"] <- "NZL"
pops[pops == "NZR"] <- "NZL"
pops[pops == "Rot"] <- "ROT"
pops[pops == "Sai"] <- "SAI"
pops[pops == "Sam"] <- "SAM"
table(pops)
combplates$Pop <- pops

# remove duplicated HAW sample
combplates[269,]$poor <- TRUE

# split oahu and maui
combplates$Pop[combplates$Pop == "HAW"] <- ifelse(combplates$plate[combplates$Pop == "HAW"] == "plate1", "MAU", "OAH")

#get data and remove poor samples (under 1mb of sequence data)
setwd("..")
m <- read.table("data/IBS/monIBS_clean.ibsMat")
m <- as.matrix(m)
bad <- which(combplates$poor)
m <- m[-bad,-bad]
combplates <- combplates[-bad,]
combplates$Pop <- factor(combplates$Pop, levels = facet.order)
combplates$color <- color.guide$color[match(combplates$Pop, color.guide$pop)]

colnames(m) <- combplates$Pop
rownames(m) <- combplates$Pop

#make and plot the tree.
nj <- nj(m) #make tree
# 
# ##get branch colors
# tcols <- rep("black", length(nj$edge.length))
# indices <- nj$edge[nj$edge[,2] <= 281, 2]
# pcols <- combplates$color[indices]
# tcols[which(nj$edge[,2] <= 281)] <- pcols
# 

#pdf("plots/NJ_tree1.pdf", width = 11, height = 8.5)
##plot
# plot.phylo(nj, no.margin = TRUE, cex = .5, type = "unrooted", 
#            tip.color = combplates$color, lab4ut = "axial", edge.color = tcols)


nj <- tidytree::as_tibble(nj)
nj$Population <- factor(full.label.tab$full[match(nj$label, full.label.tab$abrv)], levels = full.label)
nj <- tidytree::as.treedata(nj)
fsxb <- ggtree(nj, layout = "ape") + geom_tippoint(aes(color = Population), size = 4) + 
  scale_color_manual(values = pal) + ggtitle("B") + theme(legend.position = "none")


fsxbl <- g_legend(ggplot(na.omit(tidytree::as_tibble(nj)), aes(x = parent, y = node, color = Population)) + geom_point(size = 4) + theme_bw() + scale_color_manual(values = pal))


#=========NGSrelate=============

#import metadata
setwd("Raw_data/")
#import bamlist and annotations
bamlist <- read.table("../NGSadmix/gbamlist_only_good_sorted.txt")
bamlist <- cbind(bamlist, ord = 1:nrow(bamlist))
plate1 <- read.csv("plate1.csv", header = T)
plate2 <- read.csv("plate2_MFMB1.csv", header = T, skip = 1)
plate3 <- read.csv("plate3_MFMB2.csv", header = T, skip = 1)

DNA274 <- "ACAGTG"
DNA275 <- "CAGATC"

plate2 <- rbind(plate2[,1:3], setNames(plate2[,4:6], colnames(plate2)[1:3]))
plate3 <- rbind(plate3[,1:3], setNames(plate3[,4:6], colnames(plate3)[1:3]))


#add index 
barcodes <- read.table("RAD barcodes.txt", header = T)
plate2 <- merge(plate2, barcodes, by = "Well")
plate3 <- merge(plate3, barcodes, by = "Well")
plate1 <- merge(plate1, barcodes, by = "Well")


#get sample info
sampinfo <- data.frame(plate = substr(bamlist[,1], 8, 13), ind = substr(bamlist[,1], 20, 27), ord = bamlist$ord, bam = bamlist$V1, stringsAsFactors = FALSE)
sampinfo$plate[sampinfo$plate == "_split"] <- "plate1"



#add sample info to plates
plate2 <- merge(sampinfo[sampinfo$plate == DNA274,], plate2, by.x = "ind", by.y = "Index")
plate3 <- merge(sampinfo[sampinfo$plate == DNA275,], plate3, by.x = "ind", by.y = "Index")
plate1 <- merge(sampinfo[sampinfo$plate == "plate1",], plate1, by.x = "ind", by.y = "Index")

#combine plates in the correct order
combplates <- rbind(plate2, plate3, plate1)

#fix up some pop IDs
combplates <- dplyr::arrange(combplates, ord)
pops <- substr(combplates$Pop, 1, 3)
table(pops)
pops[pops == "Gua"] <- "GUA"
pops[pops == "New"] <- "NCA"
pops[pops == "Nor"] <- "NOR"
pops[pops == "NZ2"] <- "NZL"
pops[pops == "NZR"] <- "NZL"
pops[pops == "Rot"] <- "ROT"
pops[pops == "Sai"] <- "SAI"
pops[pops == "Sam"] <- "SAM"
table(pops)
combplates$Pop <- pops

# remove duplicate hawaiian sample
combplates <- combplates[-which(combplates$bam == "SOMM270_split_RA_GGGCCAAGACTGCAG.sort.flt.bam"),]

# split oahu and maui
combplates$Pop[combplates$Pop == "HAW"] <- ifelse(combplates$plate[combplates$Pop == "HAW"] == "plate1", "MAU", "OAH")


# plot
setwd("../NGSadmix/full/pop-both/")
pop <- combplates$Pop
NGSrelate <- plot_structure("combined-merged", facet.order = facet.order, clumpp = F, facet = pop, k = 1:9,
                            alt.palette = brewer.pal(9, "Set1"))
f2b <- NGSrelate$plot + xlab("Population") + theme(strip.text = element_text(size = 10))
setwd("../../..")
#ggsave("plots/NGSadmix_plot.pdf", plot = NGSrelate, device = "pdf", width = 11, height = 8.5)

f2l1 <- g_legend(f2b)
f2b <- f2b + theme(legend.position = "none")
f2b <- f2b + ggtitle("B")




# pull in and parse the likelihoods
NGSlike <- readLines("NGSadmix/full/cat_log.txt")
NGSheaders <- NGSlike[grep("Input: ", NGSlike)]
NGSlike <- data.table(K = as.numeric(gsub("nPop=", "", stringr::str_extract(NGSheaders, "nPop=[0-9]{1,2}"))),
                      run = as.numeric(gsub("_r", "", stringr::str_extract(NGSheaders, "_r[0-9]{1,2}"))),
                      est_ln_prob = as.double(stringr::str_extract(NGSlike[grep("best like", NGSlike)], "-[0-9]+\\.[0-9]+")))

# evanno it
evanno <- NGSlike[, mean(est_ln_prob), by = K]
colnames(evanno)[2] <- "mean_est_ln_prob"
evanno$lnpK <- NA
evanno$lnppK <- NA
evanno$deltaK <- NA
evanno$sd_est_ln_prob <- NGSlike[, sqrt(var(est_ln_prob)), by = K][[2]]
evanno$lnpK[-1] <- evanno$mean_est_ln_prob[-1] - evanno$mean_est_ln_prob[-nrow(evanno)]
evanno$lnppK[-nrow(evanno)] <- abs(evanno$lnpK[-nrow(evanno)] - evanno$lnpK[-1])
evanno$deltaK[-c(1, nrow(evanno))] <- abs(evanno$lnppK)[-c(1, nrow(evanno))]/evanno$sd_est_ln_prob[-c(1, nrow(evanno))] # no reason to resolve for ln''(K)
infs <- which(is.infinite(evanno$deltaK))
if(length(infs) > 0){
  evanno$deltaK[infs] <- NA
}

evanno_m <- data.table::melt(evanno, id.vars = c("K"))
evanno_m <- evanno_m[which(evanno_m$variable %in% c("deltaK", "mean_est_ln_prob")),]
evanno_m[,variable := as.character(variable)]
evanno_m$value[evanno_m$variable == "deltaK"] <- log10(evanno_m$value[evanno_m$variable == "deltaK"])
evanno_m$variable[evanno_m$variable == "deltaK"] <- "log[10](Delta*K)"
evanno_m$variable[evanno_m$variable == "mean_est_ln_prob"] <- "bar(ln(Prob))"





pdf("plots/Figure_S2.pdf", width = 11, height = 8.5)
ggplot(evanno_m[evanno_m$K <= 9,], aes(x = K, y = value)) +
  geom_point(size = 4) + 
  facet_wrap(~variable, nrow = 2, labeller = label_parsed, scales = "free_y", strip.position = "left") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        panel.grid.minor.x  = element_blank(), 
        axis.title.y = element_blank(), 
        strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text = element_text(size = 14)) +
  scale_x_continuous(breaks = 1:9)
dev.off(); dev.off()

shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S2.jpg plots/Figure_S2.pdf")


# k = 2 or 5 are the highest deltaK


#=========PCA===================
mat <- read.table("data/IBS/monIBS_clean.covMat")
mat <- as.matrix(mat)
mat <- mat[-bad,-bad]
colnames(mat) <- colnames(m)
rownames(mat) <- rownames(m)
pca_r <- eigen(mat)
# pca_r <- prcomp(mat[-dup,][-bad.samps,])
# pca <- as.data.frame(pca_r$x) #grab the PCA vectors.
pca <- pca_r$vectors
colnames(pca) <- paste0("PC", 1:ncol(pca))
pca <- as.data.frame(pca)
pca <- cbind(pca, Population = rownames(mat))
# loadings <- (pca_r$sdev^2)/sum(pca_r$sdev^2)
# loadings <- round(loadings * 100, 2)

loadings <- pca_r$values/sum(pca_r$values)
loadings <- round(loadings*100, 2)
pca$Population <- factor(pca$Population, levels = facet.order)


# PCA <- prcomp(m)
# pplot <- as.data.frame(PCA$x)
# pplot$Population <- rownames(PCA$x)
# pplot$Population <- factor(pplot$Population, levels = facet.order)
# loadings <- summary(PCA)
# loadings$importance[2,] <- round(loadings$importance[2,], 4) * 100
#ggsave("plots/PCA.pdf", plot = PCA_plot, device = "pdf", height = 8.5, width = 11)
fsxa <- ggplot(pca, aes(PC1, PC2, fill = Population)) + 
  geom_jitter(size = 4, pch = 21, width = 0.005, height = 0.005, alpha = 0.5) + theme_bw() +
  scale_fill_manual(values = pal) + 
  xlab(label = paste0("PC1 (", loadings[1], "%)")) +
  theme(legend.position = "none") + 
  ylab(label = paste0("PC2 (", loadings[2], "%)")) +
  scale_x_reverse() + ggtitle("A") # rotated to be layed out more like the geography


#==========NGSrelate based pie charts=============
K <- 5
lat_long <- list(ENA = c(19.556050, -100.289503), WNA = c(36.625980, -121.930681),
                 OAH = c(21.3069, -158.8583), MAU = c(20.9997, -155.6581),
                 GUA = c(11.421207, 142.736584),
                 ROT = c(14.154628, 145.191535), SAI = c(17.201243, 147.750705),
                 SAM = c(-13.613179, -172.351278), FIJ = c(-17.924768, 178.081698),
                 NCA = c(-21.299579, 165.383757), NOR = c(-29.024356, 167.945279),
                 QLD = c(-27.531395, 152.919356), NSW = c(-32.751896, 151.667208),
                 VIC = c(-34.574338, 138.689131), NZL = c(-37.124496, 174.961893))
# note, GUA, ROT, and SAI are fudged so they don't overlapp as much

lat_long <- as.data.frame(lat_long)
mpd <- NGSrelate
# sstab <- matrix(c("ENA", "ENA (45)",
#                   "WNA", "WNA (40)",
#                   "HAW", "HAW (12)",
#                   "GUA", "GUA (24)",
#                   "ROT", "ROT (20)",
#                   "SAI", "SAI (4)",
#                   "SAM", "SAM (31)",
#                   "FIJ", "FIJ (5)",
#                   "NCA", "NCA (18)",
#                   "NOR", "NOR (16)",
#                   "QLD", "QLD (44)",
#                   "NSW", "NSW (6)",
#                   "VIC", "VIC (4)",
#                   "NZL", "NCL (16)"), byrow = 2, ncol = 2)
# names(lat_long) <- sstab[match(names(lat_long), sstab[,1]),2]
# mpd$plot_data$pop <- sstab[match(mpd$plot_data$pop, sstab[,1]),2]

lat_long <- as.data.frame(t(lat_long))
lat_long$pop <- rownames(lat_long)
colnames(lat_long)[1:2] <- c("lat", "long")
lat_long$long[lat_long$long < 0] <- lat_long$long[lat_long$long < 0] + 360
lat_long <- sf::st_as_sf(lat_long, coords = c("long", "lat"))
lat_long <- sf::`st_crs<-`(lat_long, "EPSG:4326")


background <- sf::st_as_sf(maps::map("world2", plot = FALSE, fill = TRUE))
lat_long <- sf::st_transform(lat_long, sf::st_crs(background))

mp <- plot_structure_map(mpd, K, "pop", lat_long, sf = list(background),
                         alt.palette = RColorBrewer::brewer.pal(8, "Set1"), 
                         scale_bar = list(dist = 1000, dist_unit = "km", transform = T, st.size = 3), crop = TRUE,
                         sf_fill_colors = "white", 
                         label_args = list(max.overlaps = 10, seed = 1212, nudge_x = -7, 
                                           nudge_y = 3, point.padding = 4.3, size = 4)) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

f2a <- mp + theme(legend.position = "none", axis.text = element_blank(),
                  axis.ticks = element_blank(), axis.title = element_blank(), 
                  axis.line = element_blank(), panel.grid = element_blank()) +
  ggtitle("A")

#ggsave("plots/NGSadmix_map_plot.pdf", plot = mp, device = "pdf")



#==========================================arranged fig 1=============
pdf("plots/Figure_2.pdf", width = 11, height = 8.5)

grid.arrange(f2a, 
             f2b,
             f2l1,
             layout_matrix = rbind(c(1, 3), c(2, 3)), widths = c(1, .2), heights = c(1, 1))
dev.off(); dev.off()

shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_2.jpg plots/Figure_2.pdf")

#==========================================arranged fig SX============
pdf("plots/Figure_S1.pdf", width = 11, height = 8.5)

grid.arrange(fsxa,
             fsxb,
             fsxbl,
             layout_matrix = rbind(c(1, 3), c(2, 3)), widths = c(1, .3), heights = c(1, 1))
dev.off(); dev.off()
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S1.jpg plots/Figure_S1.pdf")


#=========het/hom ratio=======
input <- readRDS("results/paralogs/nomaf_paralog_fix_snpR.RDS")
dat <- filter_snps(input, min_loci = 0.75)
dat <- calc_het_hom_ratio(dat)
het_hom <- get.snpR.stats(dat, stats = "het_hom_ratio")$sample
colnames(het_hom)[6] <- "Population"
het_hom$Population <- factor(het_hom$Population, c("NAM", color.guide$pop[which(color.guide$pop %in% unique(het_hom$Population))]))

hhplot <- ggplot(het_hom, aes(x = Population, y = `Het/Hom`, color = Population)) + geom_point(size = 4) +
  scale_color_manual(values = c(color.guide$color[1], color.guide$color[which(color.guide$pop %in% unique(het_hom$Population))])) +
  theme_bw() + ylab("Het Count / Hom Count") + theme(axis.text.x = element_text(size = 12),
                                                     axis.text.y = element_text(size = 12),
                                                     axis.title.x = element_text(size = 14),
                                                     axis.title.y = element_text(size = 14))
hhplot
ggsave("plots/Figure_S3.pdf", plot = hhplot, device = "pdf", height = 8.5, width = 11)
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S3.jpg plots/Figure_S3.pdf")


