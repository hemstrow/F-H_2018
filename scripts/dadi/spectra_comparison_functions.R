import_spectra <- function(x){
  x <- readLines(x)
  
  # get sf
  sf <- unlist(strsplit(x[2], " "))
  sf <- as.numeric(sf)
  samps <- as.numeric(unlist(strsplit(x[1], " "))[1:2])
  sf <- matrix(sf, samps[2], samps[1])
  
  # get mask
  mask <- as.numeric(unlist(strsplit(x[3], " ")))
  mask <- matrix(mask, samps[2], samps[1])
  
  # mask
  sf[which(mask == 1)] <- NA
  
  # add pops attribute
  pops <- unlist(strsplit(x[1], " "))[4:5]
  if(!is.na(pops[1])){
    pops <- gsub("\"", "", pops)
    attr(sf, "pops") <- rev(pops)
  }
  return(sf)
}
plot_comp_spectra <- function(obs, mod, resid, viridis.option = "inferno",
                              resid.viridis.option = "magma"){

  # reformat using snpR, but need to combine plots on constant scale so can't
  # just take these
  p1 <- snpR::plot_sfs(obs)
  attr(mod, "pop") <- attr(obs, "pop")
  p2 <- snpR::plot_sfs(mod)
  
  # combine reformatted data
  p1d <- p1$data
  p2d <- p2$data
  p1d$source <- "data"
  p2d$source <- "model"
  pdc <- rbind(p1d, p2d)
  pops <- attr(obs, "pop")
  
  # plot the two spectra
  pc <- ggplot2::ggplot(pdc, ggplot2::aes(x = p1, y = p2, fill = log10(N))) + 
    ggplot2::geom_tile() +
    ggplot2::facet_grid(.~source) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_viridis_c(na.value = "white", option = viridis.option) +
    ggplot2::scale_fill_viridis_c(na.value = "white", option = viridis.option) +
    ggplot2::xlab(pops[1]) + ggplot2::ylab(pops[2]) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(strip.background = element_blank())
  
  pc_leg <- ggpubr::get_legend(pc)
  
  pc <- pc + ggplot2::theme(legend.position = "none")
  
  # get and plot the residuals
  attr(resid, "pop") <- attr(obs, "pop")
  pr <- snpR::plot_sfs(resid, log = F)
  pr <- suppressMessages(pr + 
                           ggplot2::scale_color_viridis_c(na.value = "white", 
                                                          direction = 1,
                                                          option = resid.viridis.option) +
                           ggplot2::scale_fill_viridis_c(na.value = "white", 
                                                         direction = 1,
                                                         option = resid.viridis.option) +
                           ggplot2::labs(fill = "Residual", 
                                         color = "Residual"))
  
  pr_leg <- ggpubr::get_legend(pr)
  
  pr <- pr + ggplot2::theme(legend.position = "none")
  
  # plot residual histogram
  hist <- ggplot2::ggplot(data.frame(resid = as.numeric(resid)), aes(x = resid)) +
    geom_histogram(color = "steelblue", fill = "steelblue") + 
    theme_bw() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .1))) +
    ggplot2::xlab("Residuals")
  
  return(list(comp = pc, pc_legend = pc_leg, resid = pr, resid_legend = pr_leg, resid_hist = hist))
}
arrange_comp_spectra <- function(out_plots){
  gridExtra::grid.arrange(out_plots$comp, out_plots$pc_legend, 
                          out_plots$resid, out_plots$resid_legend,
                          out_plots$resid_hist,
                          layout_matrix = rbind(c(1, 1, 2), c(5, 3, 4)),
                          nrow = 2, ncol = 3, widths = c(1, 1, .2))
}

pick_and_plot_comp <- function(selected_ilist, rdf, quant_vars, quant_cuts = c(.5, .5),
                               log_vars = c(T, T),
                               projection,
                               polarized,
                               model = NULL,
                               viridis.option = "inferno",
                               resid.viridis.option = "magma"){

  plot_run_IDs <- vector("list", length = 4)
  thresholds <- vector("list", length = 2)
  
  # note the vars in each quantile (I know this is repetitive, but I'm too lazy right now to make it cleaner...)
  if(log_vars[1]){
    thresholds[[1]] <- quantile(log10(selected_ilist[,quant_vars[1]]), quant_cuts[1])
    cat("Var 1 threshold: ", thresholds[[1]], "\n")
    thresholds[[1]] <- which(log10(selected_ilist[, quant_vars[1]]) > thresholds[[1]])
  }
  else{
    thresholds[[1]] <- quantile(selected_ilist[,quant_vars[1]], quant_cuts[1])
    cat("Var 1 threshold: ", thresholds[[1]], "\n")
    thresholds[[1]] <- which(selected_ilist[, quant_vars[1]] > thresholds[[1]])
  }
  if(log_vars[2]){
    thresholds[[2]] <- quantile(log10(selected_ilist[,quant_vars[2]]), quant_cuts[2])
    cat("Var 2 threshold: ", thresholds[[2]], "\n")
    thresholds[[2]] <- which(log10(selected_ilist[, quant_vars[2]]) > thresholds[[2]])
  }
  else{
    thresholds[[2]] <- quantile(selected_ilist[,quant_vars[2]], quant_cuts[2])
    cat("Var 2 threshold: ", thresholds[[2]], "\n")
    thresholds[[2]] <- which(selected_ilist[, quant_vars[2]] > thresholds[[2]])
  }
  
  
  # fetch run IDs for minimum AIC scores in each quadrant
  plot_run_IDs[[1]] <- selected_ilist[-sort(union(thresholds[[1]], thresholds[[2]])),c("run_ID", "pass")][which.min(selected_ilist[-sort(union(thresholds[[1]], thresholds[[2]])),]$AIC)[1],] # bottom left
  plot_run_IDs[[2]] <- selected_ilist[thresholds[[1]][-which(thresholds[[1]] %in% thresholds[[2]])],c("run_ID", "pass")][which.min(selected_ilist[thresholds[[1]][-which(thresholds[[1]] %in% thresholds[[2]])],]$AIC)[1],] # bottom right
  plot_run_IDs[[3]] <- selected_ilist[thresholds[[2]][-which(thresholds[[2]] %in% thresholds[[1]])],c("run_ID", "pass")][which.min(selected_ilist[thresholds[[2]][-which(thresholds[[2]] %in% thresholds[[1]])],]$AIC)[1],] # top left
  plot_run_IDs[[4]] <- selected_ilist[sort(intersect
                                           (thresholds[[1]], thresholds[[2]])),c("run_ID", "pass")][which.min(selected_ilist[sort(intersect(thresholds[[1]], thresholds[[2]])),]$AIC)[1],] # top right
  plot_run_IDs <- unlist(lapply(plot_run_IDs, paste, collapse = " "))
  
  # fetch matching AICs from data to get the runs to plot
  rdf <- data.table::as.data.table(rdf)
  rdf2 <- rdf[,unique_ids := do.call(paste, .SD), .SDcols = c("parms", "pass")]
  plot_parms <- rdf2[which(rdf2$unique_ids %in% plot_run_IDs),]
  plot_parms <- plot_parms[match(plot_parms$unique_ids, plot_run_IDs),] # order to match
  
  # call the python script and plot each spectra
  ## prep other parms
  projection <- paste0("[", paste0(projection, collapse = ","), "] ")
  if(is.null(model)){
    model <- selected_ilist[1, "model"]
  }
  model <- paste0(model, " ")
  polarized <- ifelse(polarized, "True", "False")
  
  out_plots <- vector("list", 4)
  raw_plots <- out_plots
  for(i in 1:length(out_plots)){
    cmd <- paste0("python scripts/dadi/save_2d_spectra.py [",
                  paste0(unlist(strsplit(plot_parms[i,]$mnum, " ")), collapse = ","), "] ",
                  model,
                  "data/dadi_inputs/dadi_10kgap_snps.txt ",
                  "[", paste0(unlist(strsplit(plot_parms[i,]$pops, " ")), collapse = ","), "] ",
                  projection,
                  polarized
    )
    shell(cmd)
    raw_plots[[i]] <- plot_comp_spectra(import_spectra("obs_spectra.txt"), 
                                        import_spectra("modeled_spectra.txt"),
                                        import_spectra("residual.txt"))
    
    out_plots[[i]] <- ggpubr::as_ggplot(arrange_comp_spectra(raw_plots[[i]]))
    Sys.sleep(2) # tends to crash without this -- maybe rendering lag or something?
  }
  return(list(plots = out_plots, plot_point_ids = plot_run_IDs, raw_plots = raw_plots))
}