library(RColorBrewer); library(ape); library(snpR); library(ggplot2); library(gridExtra); library(gridGraphics)

# This script makes all of the R plots for the publication. Plots are as follows, 
# asterisks imply that the plot was not made with R:
# 
# Figure 1: Map of sampling locations. *




facet.order <- c("ENA", "WNA", "HAW", "GUA", "ROT", "SAI", "SAM", "FIJ", "NCA", "NOR", "QLD", "NSW", "VIC", "NZL")
pal <- colorRampPalette(brewer.pal(9, "Set1"))(14) #palette to use
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


#get data and remove poor samples (under 1mb of sequence data)
setwd("..")
m <- read.table("data/IBS/monIBS_clean.ibsMat")
m <- as.matrix(m)
m <- m[-which(combplates$poor),-which(combplates$poor)]
combplates <- combplates[-which(combplates$poor),]
pal <- colorRampPalette(brewer.pal(9, "Set1"))(14) #palette to use
combplates$Pop <- factor(combplates$Pop, levels = facet.order)
combplates$color <- color.guide$color[match(combplates$Pop, color.guide$pop)]

colnames(m) <- combplates$Pop
rownames(m) <- combplates$Pop

#make and plot the tree.
nj <- nj(m) #make tree

##get branch colors
tcols <- rep("black", length(nj$edge.length))
indices <- nj$edge[nj$edge[,2] <= 281, 2]
pcols <- combplates$color[indices]
tcols[which(nj$edge[,2] <= 281)] <- pcols


#pdf("plots/NJ_tree1.pdf", width = 11, height = 8.5)
##plot
plot.phylo(nj, no.margin = TRUE, cex = .5, type = "unrooted", 
           tip.color = combplates$color, lab4ut = "axial", edge.color = tcols)
grid.echo()
f2c <- grid.grab()


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


# plot
setwd("../NGSadmix/full/pop-both/")
pop <- combplates$Pop
NGSrelate <- plot_structure("combined-merged", facet.order = facet.order, clumpp = F, facet = pop, k = 9,
                            alt.palette = brewer.pal(9, "Set1"))
f2b <- NGSrelate$plot + xlab("Population")
setwd("../../..")
#ggsave("plots/NGSadmix_plot.pdf", plot = NGSrelate, device = "pdf", width = 11, height = 8.5)

f2l2 <- g_legend(f2b)
f2b <- f2b + theme(legend.position = "none")

#=========PCA===================
PCA <- prcomp(m)
pplot <- as.data.frame(PCA$x)
pplot$Population <- rownames(PCA$x)
pplot$Population <- factor(pplot$Population, levels = facet.order)
loadings <- summary(PCA)
loadings$importance[2,] <- round(loadings$importance[2,], 4) * 100
f2a <- ggplot(pplot, aes(PC1, PC2, color = Population)) + geom_point() + theme_bw() +
  scale_color_manual(values = color.guide$color) + xlab(label = paste0("PC1 (", loadings$importance[2,1], "%)")) +
  ylab(label = paste0("PC2 (", loadings$importance[2,2], "%)"))
#ggsave("plots/PCA.pdf", plot = PCA_plot, device = "pdf", height = 8.5, width = 11)
f2l1 <- g_legend(f2a)
f2a <- ggplot(pplot, aes(PC1, PC2, color = Population)) + geom_point() + theme_bw() +
  scale_color_manual(values = color.guide$color) + xlab(label = paste0("PC1 (", loadings$importance[2,1], "%)")) +
  theme(legend.position = "none") + 
  ylab(label = paste0("PC2 (", loadings$importance[2,2], "%)")) +
  scale_y_reverse() + scale_x_reverse()


#==========NGSrelate based pie charts not used (but cool!)=============
library(scatterpie)
K <- 5
lat_long <- list(ENA = c(19.556050, -100.289503), WNA = c(36.625980, -121.930681),
                 HAW = c(19.627274, -155.493135), GUA = c(11.421207, 142.736584),
                 ROT = c(14.154628, 145.191535), SAI = c(17.201243, 147.750705),
                 SAM = c(-13.613179, -172.351278), FIJ = c(-17.924768, 178.081698),
                 NCA = c(-21.299579, 165.383757), NOR = c(-29.024356, 167.945279),
                 QLD = c(-27.531395, 152.919356), NSW = c(-32.751896, 151.667208),
                 VIC = c(-34.574338, 138.689131), NZL = c(-37.124496, 174.961893))
# note, GUA, ROT, and SAI are fudged so they don't overlapp as much

lat_long <- as.data.frame(lat_long)
mp <- plot_structure_map(NGSrelate, K, "pop", lat_long, alt.palette = RColorBrewer::brewer.pal(8, "Set1"))
f2d <- mp + theme(legend.position = "none")

#ggsave("plots/NGSadmix_map_plot.pdf", plot = mp, device = "pdf")



#==========================================arranged fig 2=============
pdf("plots/Figure_2.pdf", width = 15, height = 8.5)
grid.arrange(f2c, 
             f2b + 
               ggplot2::theme(axis.text.x = element_text(size = 9.5), 
                              axis.text.y = element_blank(), 
                              axis.ticks.y = element_blank()), 
             f2a + ggplot2::theme(axis.text = element_blank(), axis.ticks = element_blank()),
             f2d, f2l1, f2l2, 
             layout_matrix = cbind(c(3,NA,1), c(5,5,5), c(2,NA,4), c(6,6,6)), widths = c(1, .2, 1, .2), heights = c(1,.1,1))
dev.off(); dev.off()

#=========het/hom ratio=======
input <- readRDS("results/paralogs/nomaf_paralog_fix_snpR.RDS")
dat <- filter_snps(input, min_loci = 0.75)
dat <- calc_het_hom_ratio(dat)
het_hom <- dat@sample.stats
colnames(het_hom)[2] <- "Population"
colnames(het_hom)[5] <- "Het_Hom"
het_hom$Population <- factor(het_hom$Population, c("NAM", color.guide$pop[which(color.guide$pop %in% unique(het_hom$Population))]))

hhplot <- ggplot(het_hom, aes(x = Population, y = Het_Hom, color = Population)) + geom_point() +
  scale_color_manual(values = c(color.guide$color[1], color.guide$color[which(color.guide$pop %in% unique(het_hom$Population))])) +
  theme_bw() + ylab("Het Count / Hom Count")
ggsave("plots/het_hom.pdf", plot = hhplot, device = "pdf", height = 8.5, width = 11)


