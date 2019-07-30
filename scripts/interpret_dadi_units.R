interpret_units <- function(wlist, mu, g, L){
  mu <- mu*L
  
  for(i in 1:length(wlist)){
    wlist[[i]]$theta  <- wlist[[i]]$theta/(4*mu)
    
    start.col <- which(colnames(wlist[[i]]) == "AIC") + 1
    for(j in start.col:ncol(wlist[[i]])){
      t.col.name <- colnames(wlist[[i]])[j]
      
      if(grepl("^nu", t.col.name) | t.col.name %in% c("s", "f") | grepl("^K", t.col.name)){
        wlist[[i]][j] <- wlist[[i]]$theta * wlist[[i]][j]
      }
      else if(grepl("^m", t.col.name)){
        wlist[[i]][j] <- wlist[[i]][j]/(wlist[[i]]$theta*2)
      }
      else if(grepl("^T", t.col.name)){
        wlist[[i]][j] <- 2 * wlist[[i]]$theta * wlist[[i]][j] * g
      }
    }
  }
  
  return(wlist)
}