find_non_paralogous_sections <- function(paralog_dir, paralog_file_pattern, unique_chromsome_file, base_pair_buffer){
  
  flist <- list.files(paralog_dir, paralog_file_pattern)
  flist <- paste0(paralog_dir, flist)
  dat <- vector("list", length(flist))
  for(i in 1:length(dat)){
    dat[[i]] <- readr::read_table2(flist[i], col_names = F)
  }
  dat <- dplyr::bind_rows(dat)
  
  
  #====================prepare paralog list=====================
  # filter to likely paralogs
  paralog_list <- dat[dat[,5] > 10, 1:2]
  rm(dat)
  # unique entries
  paralog_list <- unique(paralog_list)
  colnames(paralog_list) <- c("chr", "position")
  paralog_list <- arrange(paralog_list, chr, position) # sort
  
  # now need to highlight the regions to run. These exclude each site within 1kb of a paralog (basically the tag)
  # to do so, loop through each possible chr and find the acceptable regions
  chr_opts <- read.table(unique_chromsome_file, header = F, stringsAsFactors = F)
  chr_lengths <- chr_opts[,2]
  chr_opts <- chr_opts[,1]
  chr_opts <- gsub(">", "", chr_opts)
  
  # loop through each possible chromosome, building a series of arguments listing "ok" sites
  out <- vector("list", length(chr_opts))
  names(out) <- chr_opts
  
  
  good.sections <- character(0)
  for(i in 1:length(chr_opts)){
    # if this chr is all good, note and skip
    if(!chr_opts[i] %in% paralog_list$chr){
      out[[i]] <- paste0(chr_opts[i], ":1-")
      next
    }
    
    # initialize
    t.chr.bads <- unlist(paralog_list[which(paralog_list$chr == chr_opts[i]), 2])
    t.chr.bads <- sort(t.chr.bads)
    b.end <- 0
    g.start <- 1
    
    # for each bad entry, need to set a 1kb restriction area
    for(j in 1:length(t.chr.bads)){
      
      # if this is the last one, everything after this bad window is good
      if(j == length(t.chr.bads)){
        # as long as there are actually more bps in the chr, note that everything else is good
        if(chr_lengths[j] <= t.chr.bads[j] + base_pair_buffer + 1){
          good.sections <- c(good.sections, paste0(chr_opts[i], ":", t.chr.bads[j] + base_pair_buffer + 1, "-"))
        }
        next
      }
      
      # if this bad window start is greater than the previous bad window end, time to print a new good section (unless there is a perfect overlap!)
      if(t.chr.bads[j] - base_pair_buffer > b.end){
        if(g.start < t.chr.bads[j] - base_pair_buffer - 1){
          good.sections <- c(good.sections, paste0(chr_opts[i], ":", g.start, "-", t.chr.bads[j] - base_pair_buffer - 1))
        }
        g.start <- t.chr.bads[j] + base_pair_buffer + 1
        b.end <- t.chr.bads[j] + base_pair_buffer
      }
      
      # otherwise we are still in the same window
      else{
        b.end <- t.chr.bads[j] + base_pair_buffer
        g.start <- t.chr.bads[j] + base_pair_buffer + 1
        next()
      }
    }
  }
  
  return(good.sections)
}