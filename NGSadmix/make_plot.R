library(gridExtra); library(RColorBrewer); library(pophelper)

#import metadata
setwd("../Raw_data/")
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

#################

#make plots

library(pophelper)

setwd("../NGSadmix/full/")

#prepare CLUMPP outputfiles
qfiles <- list.files(full.names = T, pattern = "qopt")
qlist <- readQ(qfiles)
#clumppExport(qlist, useexe = T)
collectClumppOutput(filetype = "both")

#prepare palette
cbp <- brewer.pal(9, "Set1")

#metadata
spmeta <- data.frame(ind = combplates$ID,
                     pop = combplates$Pop,
                     stringsAsFactors = F)

spmeta <- spmeta[order(pops),]

#get the clummp data
mq <- readQ(list.files("pop-both/", full.names = T, pattern = "merged"))

#re-sort individuals by population
for(i in 1:length(mq)){
  mq[[i]] <- mq[[i]][order(pops),]
}


#fix things so that the cluster order is the same in all sets
##function to do this... include in snpR somewhere? It's handy.
fix_clust <- function(x){
  
  #loop through each q object
  for (i in 2:length(x)){
    #see which columns in the previous run are the most similar to each column
    
    #initialize mapping df
    mdf <- data.frame(tcol = 1:ncol(x[[i]]), pcol = numeric(ncol(x[[i]])),
                      ed = numeric(ncol(x[[i]])))
    
    #loop through each column and find where to map it.
    for (j in 1:ncol(x[[i]])){
      
      #intialize euc distance vector
      elist <- numeric(ncol(x[[i - 1]]))
      
      #compare to each other col.
      for(k in 1:ncol(x[[i-1]])){
        #save euclidian dist
        elist[k] <- sum((x[[i]][,j] - x[[i-1]][,k])^2)
      }
      
      #save results
      mdf[j,2] <- which.min(elist)
      mdf[j,3] <- min(elist)
    }
    
    #reassign clusters in this qdf
    ##which is the new cluster? Probably that with the most distance to any original clusters.
    dups <- duplicated(mdf[,2]) | duplicated(mdf[,2], fromLast = T)
    nc <- which.max(mdf[dups,3])
    mdf[dups,2][nc] <- nrow(mdf)
    mdf <- mdf[order(mdf[,2]),]
    
    ##reasign clusters
    tdf <- x[[i]]
    tdf <- tdf[,mdf[,1]]
    
    ##replace object in x with the re-arranged qfile.
    colnames(tdf) <- colnames(x[[i]])
    x[[i]] <- tdf 
  }
  
  return(x)
}
##run it
mqf <- fix_clust(mq)


#pops
pop <- as.data.frame(spmeta[,2], stringsAsFactors = F)
colnames(pop) <- "pop"


#sort the individuals within each population based on the qvals in the last 
##function to do this
sort_inds <- function(x, pop, cluster = "first", q = "last"){
  
  #get which pop to use
  if(q == "last"){
    q <- length(x)
  }
  
  #get order to stick individual in
  lx <- x[[q]]
  upops <- unique(pop)
  lx$s <- 1:nrow(lx)
  
  #get the sorting cluster priority:
  if(cluster == "first"){
    cseq <- (ncol(lx)-1):1
  }
  else if (cluster == "last"){
    cseq <- 1:(ncol(lx)-1)
  }
  else if (is.numeric(cluster)){
    if(length(cluster) == ncol(lx)){
      cseq <- cluster
    }
    else{
      if(length(cluster) < ncol(lx)){
        cseq <- c((1:ncol(lx))[-which(1:ncol(lx) %in% cluster)], rev(cluster))
      }
      else{
        stop("Cluster length is longer than number of clusters in x element q.\n")
      }
    }
  }
  
  for(i in 1:nrow(upops)){
    tx <- lx[which(pop == upops[i,]),]
    for(j in cseq){
      tx <- tx[order(tx[,j]),]
    }
    lx[which(pop == upops[i,]),] <- tx
  }
  
  #order all the datasets like this.
  ord <- lx$s
  for(i in 1:(length(x))){
    x[[i]] <- x[[i]][ord,]
  }
  
  return(x)
}

mqs <- sort_inds(mqf, pop, cluster = c(6,7,8), q = 7) #using the k = 6 plot, which is the last one where anything new shows up, sort by clusters 4(UPD/ASP), 5(CLF), then 6 (OPL)


arr.pops <- c("ENA", "WNA", "HAW", "GUA", "ROT", "SAI", "SAM", "FIJ", "NCA", "NOR", "QLD", "NSW", "VIC", "NZL") #specify order for sorting populations in figure

#prepare plot
p <- plotQ(mqs,
           returnplot = T, exportplot = F, imgoutput = "join", clustercol = cbp,
           grplab = pop,
           grplabsize = 3, grplabcol = "black", splabcol = "black", splabsize = 10, 
           grplabangle = 90, grplabpos = 1, grplabheight = 1,
           splab = paste0("K=", 2:9), pointsize = 8, divsize = 1, 
           ordergrp = T, subsetgrp = arr.pops)



#plotQ(readQ("pop_K4/pop_K4-combined-merged.txt"),
#      returnplot = T, exportplot = F, basesize = 11)


setwd("../../")