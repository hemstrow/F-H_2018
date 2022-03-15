#======================set parameters==================
source("scripts/dadi/dadi_reading_functions.R")
library(dplyr); library(ggplot2)
mu <- 8.4e-9
g <- .3
ratio <- 9370/302446 # ratio of included snps
L <-  1373747*ratio # approx number of considered bases.
wb <- openxlsx::createWorkbook()

#======================read in and summarize results===========================
pass1 <- import_dadi_results("data/dadi_inputs/cat_NH_unfolded_r1.txt",
                             mu = mu, g = g, L = L)
pass2 <- import_dadi_results("data/dadi_inputs/cat_NH_unfolded_r2.txt",
                             mu = mu, g = g, L = L)
pass3 <- import_dadi_results("data/dadi_inputs/cat_NH_unfolded_r3.txt",
                             mu = mu, g = g, L = L)
pass4 <- import_dadi_results("data/dadi_inputs/cat_NH_unfolded_r4.txt", 
                             mu = mu, g = g, L = L)  


# AIC dist per model
rdf <- rbind(cbind(pass1$rdf, pass = 1),
             cbind(pass2$rdf, pass = 2),
             cbind(pass3$rdf, pass = 3),
             cbind(pass4$rdf, pass = 4))
rdf$pass <- as.factor(rdf$pass)

rdf$model[rdf$model == "founder_asym_growth_both"] <- "Found and Grow"
rdf$model[rdf$model == "founder_asym_hist_3epoch_exp_growth_p1"] <- "Three Epoch"
rdf$model[rdf$model == "founder_asym_hist_igrowth_p2"] <- "Zhan"
rdf$model[rdf$model == "founder_nomig_admix_two_epoch_growth_both"] <- "Two Epoch"

# minimum AIC per model
mins <- tapply(rdf$AIC, rdf$model, min)
sort(mins)

# plot AIC scores in each model in each pass
pdf("./plots/Figure_S5.pdf")
ggplot(rdf,aes(x = model, y = log(AIC), color = pass)) + 
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Set1"))
dev.off(); dev.off()
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S5.jpg plots/Figure_S5.pdf")


# best overall model
best.reps <- tapply(rdf$AIC, rdf[,c("model", "pass")], min)
best.reps[order(best.reps[,4]),] # best model in pass 4
best.reps[which(rowSums(best.reps == min(best.reps, na.rm = T)) == 1),,drop = F] # best overall


# model averages and summaries
aic_means <- tapply(rdf$AIC, rdf[,c("model", "pass")], mean, na.rm = T)
aic_mins <- tapply(rdf$AIC, rdf[,c("model", "pass")], min, na.rm = T)
aic_sd <- tapply(rdf$AIC, rdf[,c("model", "pass")], sd)
ll_means <- tapply(rdf$ll, rdf[,c("model", "pass")], mean)
ll_max <- tapply(rdf$ll, rdf[,c("model", "pass")], max)
ll_sd <- tapply(rdf$ll, rdf[,c("model", "pass")], sd)
n_runs <- table(rdf[,c("model", "pass")])

comb <- list(aic_mean = aic_means, aic_min = aic_mins, aic_sd = aic_sd, 
             ll_mean = ll_means, ll_max = ll_max, ll_sd = ll_sd, 
             n = n_runs)
for(i in 1:length(comb)){
  comb[[i]] <- reshape2::melt(comb[[i]])
  comb[[i]]$var <- names(comb)[i]
}
comb <- dplyr::bind_rows(comb)
dcomb <- reshape2::dcast(comb, model ~ var + pass)

openxlsx::addWorksheet(wb, "model_statistics")
openxlsx::writeDataTable(wb, "model_statistics", dcomb, tableStyle = "none")
openxlsx::saveWorkbook(wb, "plots/Table_S4.xlsx", overwrite = T)

#=======================grab and prepare data for three best models and Zahn model=============================
# Three Epoch
ilist <- rbind(cbind(pass1$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 1),
               cbind(pass2$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 2),
               cbind(pass3$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 3),
               cbind(pass4$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 4))
ilist$pass <- as.factor(ilist$pass)

# Found and grow
ilist2 <- rbind(cbind(pass1$ilist$founder_asym_growth_both_NAM_HAW, pass = 1),
                cbind(pass2$ilist$founder_asym_growth_both_NAM_HAW, pass = 2),
                cbind(pass3$ilist$founder_asym_growth_both_NAM_HAW, pass = 3),
                cbind(pass4$ilist$founder_asym_growth_both_NAM_HAW, pass = 4))
ilist2$pass <- as.factor(ilist2$pass)

# Two Epoch
ilist3 <- rbind(cbind(pass1$ilist$founder_nomig_admix_two_epoch_growth_both_NAM_HAW, pass = 1),
                cbind(pass2$ilist$founder_nomig_admix_two_epoch_growth_both_NAM_HAW, pass = 2),
                cbind(pass3$ilist$founder_nomig_admix_two_epoch_growth_both_NAM_HAW, pass = 3),
                cbind(pass4$ilist$founder_nomig_admix_two_epoch_growth_both_NAM_HAW, pass = 4))
ilist3$pass <- as.factor(ilist3$pass)

# Zahn
ilist4 <- rbind(cbind(pass1$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 1),
                cbind(pass2$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 2),
                cbind(pass3$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 3),
                cbind(pass4$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 4))
ilist4$pass <- as.factor(ilist4$pass)

# combine data
## merge: 1 and 2
ilist1.1 <- ilist[,c("model", "AIC", "nuA", "Tg", "Ts", "Tg2", "Tg3", "s", "nu1F", "nu2F", "nuG2", "m12", "m21", "pass", "run_ID", "pops")]
ilist1.1$AncientTime <- ilist1.1$Tg + ilist1.1$Ts + ilist1.1$Tg2 + ilist1.1$Tg3
ilist1.1$SplitTime <- ilist1.1$Ts + ilist1.1$Tg3
colnames(ilist1.1)[c(9,10,11)] <- c("Ne_NA", "Ne_Ha", "Ne_Split")
ilist2.1 <- ilist2[,c("model", "AIC", "nuA", "nu1", "nu2", "m12", "m21", "Ti", "s", "pass", "run_ID", "pops")]
colnames(ilist2.1)[c(4,5,8)] <- c("Ne_NA", "Ne_Ha", "SplitTime")
ilist2.1$Ne_Split <- ilist2.1$nuA
ilistm <- merge(ilist1.1, ilist2.1, all = T)
ilistm$model <- ifelse(ilistm$model == "founder_asym_growth_both", "Found and Grow", "Three Epoch")

## merge in 3
ilist3.1 <- ilist3
ilist3.1$SplitTime <- ilist3.1$T1 + ilist3.1$T2
colnames(ilist3.1)[7:8] <- c("Ne_NA", "Ne_Ha")
ilistm <- merge(ilistm, ilist3.1, all = T)
ilistm$model[ilistm$model == "founder_nomig_admix_two_epoch_growth_both"] <- "Two Epoch"

## merge in 4
ilist4.1 <- ilist4
ilist4.1$SplitTime <- ilist4.1$Ts + ilist4.1$Tg2
ilist4.1$AncientTime <- ilist4.1$Tg + ilist4.1$Tg2 + ilist4.1$Ts
ilist4.1$Ne_NA <- ilist4.1$nuG - ilist4.1$s
ilist4.1$Ne_Split <- ilist4.1$nuG
ilist4.1$Ne_Ha <- ilist4.1$nu2F
ilistm <- merge(ilist4.1, ilistm, all = T)
ilistm[ilistm$model == "founder_asym_hist_igrowth_p2",]$model <- "Zhan"




#======================plot some spectra=====================
source("scripts/dadi/spectra_comparison_functions.R")
# plot spectra in each quadrant for each model
exit <- T
reticulate::repl_python() # need to do this to add python to the path, apparently
exit

quant_cuts <- c(.5, .5)

## 1
spectra1 <- pick_and_plot_comp(selected_ilist = ilistm[ilistm$model == "Three Epoch",], 
                               rdf = rdf, 
                               quant_vars = c("s", "SplitTime"), 
                               quant_cuts = quant_cuts, 
                               log_vars = c(T, T), 
                               projection = c(100, 10), 
                               polarized = T, 
                               model = "founder_asym_hist_3epoch_exp_growth_p1")

## 2
spectra2 <- pick_and_plot_comp(selected_ilist = ilistm[ilistm$model == "Found and Grow",], 
                               rdf = rdf, 
                               quant_vars = c("s", "SplitTime"), 
                               quant_cuts = quant_cuts, 
                               log_vars = c(T, T), 
                               projection = c(100, 10), 
                               polarized = T, 
                               model = "founder_asym_growth_both")

## 3
spectra3 <- pick_and_plot_comp(selected_ilist = ilistm[ilistm$model == "Two Epoch",], 
                               rdf = rdf, 
                               quant_vars = c("s", "SplitTime"), 
                               quant_cuts = quant_cuts, 
                               log_vars = c(T, T), 
                               projection = c(100, 10), 
                               polarized = T, 
                               model = "founder_nomig_admix_two_epoch_growth_both")


##4
spectra4 <- pick_and_plot_comp(selected_ilist = ilistm[ilistm$model == "Zhan",], 
                               rdf = rdf, 
                               quant_vars = c("s", "SplitTime"), 
                               quant_cuts = quant_cuts, 
                               log_vars = c(T, T), 
                               projection = c(100, 10), 
                               polarized = T, 
                               model = "founder_asym_hist_igrowth_p2")

#===============prepare plots: model comparisons===============
## set the correct points to highlight (those that were used for spectra comparisons)
ilistm$point_ids <- paste(ilistm$run_ID, ilistm$pass)
ilistm$highlight <- ifelse(ilistm$point_ids %in% c(spectra1$plot_point_ids, 
                                                   spectra2$plot_point_ids,
                                                   spectra3$plot_point_ids,
                                                   spectra4$plot_point_ids), 1, 0)

ilistm$model <- factor(ilistm$model, c("Three Epoch", "Found and Grow", "Two Epoch", "Zhan"))


## grab legend
dpb.time <- ggplot(ilistm, aes(y = log10(SplitTime), x = log10(s), color = log10(AIC), shape = pass)) + 
  geom_point(size = 4) + geom_jitter() + theme_bw() + scale_color_viridis_c() + facet_wrap(~model) +
  ylab("log10(Time since Hawaii establishment)") + xlab("log10(Ne Hawaii Establishment)") +
  labs(color = bquote(~log[10]~ '(AIC)')) + scale_shape_manual(values = 15:18)

legend <- ggpubr::get_legend(dpb.time)

## make plots
dpb.time <- ggplot(ilistm, aes(y = log10(SplitTime), x = log10(s), color = log10(AIC), shape = pass)) + 
  geom_point() + theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_blank()) + 
  scale_color_viridis_c() + facet_wrap(~model, ncol = 4) +
  geom_point(data = ilistm[which(ilistm$highlight == 1),], 
             aes(y = log10(SplitTime), x = log10(s)), color = "red", size = 3) +
  ylab(bquote(~log[10]~ '(Years Since Split)')) + xlab(bquote('Establishment' ~log[10]~ '(' ~N[e]~ '), Hawaii')) + 
  scale_y_continuous(labels = function(x) sprintf("%.5f", x)) + scale_shape_manual(values = 15:18)

dpb.end.size <- ggplot(ilistm, aes(x = log10(Ne_NA), y = log10(Ne_Ha), color = log10(AIC), shape = pass)) + 
  geom_point() + theme_bw()+ 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_color_viridis_c() + facet_wrap(~model, ncol = 4) +
  geom_point(data = ilistm[which(ilistm$highlight == 1),], aes(x = log10(Ne_NA), y = log10(Ne_Ha)), color = "red", size = 3) +
  xlab(bquote('Current' ~log[10]~ '(' ~N[e]~ '), North America')) + ylab(bquote('Current' ~log[10]~ '(' ~N[e]~ '), Hawaii')) + 
  scale_y_continuous(labels = function(x) sprintf("%.5f", x)) + scale_shape_manual(values = 15:18)

dpb.mig <- ggplot(ilistm, aes(x = m12, y = m21, color = log10(AIC), shape = pass)) + 
  geom_point() + theme_bw()+ 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +scale_color_viridis_c() + facet_wrap(~model, ncol = 4) +
  geom_point(data = ilistm[which(ilistm$highlight == 1),], aes(x = m12, y = m21), color = "red", size = 3) +
  xlab("Migration: Hawaii -> NA") + ylab("Migration: NA -> Hawaii") +
  scale_y_continuous(labels = function(x) sprintf("%.5f", x)) + scale_shape_manual(values = 15:18)

## combine plots
pdf("plots/Figure_3.pdf", width = 11, height = 8.5)
gridExtra::grid.arrange(dpb.time, dpb.mig, dpb.end.size, legend,
                        layout_matrix = matrix(c(1,2,3,4,4,4), nrow = 3, ncol = 2),
                        widths = c(1,.1))
dev.off();dev.off();
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_3.jpg plots/Figure_3.pdf")

#=================prepare plots: spectra=================
# want four graphs: one with the real, one with the models, one with the resid hists, and one with the resid plot
## get heatmap data, residuals, etc all combined
all_spectra <- list(spectra1, spectra2, spectra3, spectra4)
all_spectra_dat <- list("vector", 4*4)
residuals <- all_spectra_dat
residual_heatmap <- all_spectra_dat

model_key <- data.frame(i = 1:4, model = c("Three Epoch",
                                           "Found and Grow",
                                           "Two Epoch",
                                           "Zhan"))

quad_key <- data.frame(j = 1:4, c("BL", "BR", "TL", "TR"))
tracker <- 1
for(i in 1:length(all_spectra)){
  for(j in 1:length(all_spectra[[i]]$raw_plots)){
    all_spectra_dat[[tracker]] <- all_spectra[[i]]$raw_plots[[j]]$comp$data
    all_spectra_dat[[tracker]]$model <- model_key[i, 2]
    all_spectra_dat[[tracker]]$quadrant <- quad_key[j, 2]
    
    residuals[[tracker]] <- all_spectra[[i]]$raw_plots[[j]]$resid_hist$data
    residuals[[tracker]]$model <- model_key[i, 2]
    residuals[[tracker]]$quadrant <- quad_key[j, 2]
    
    residual_heatmap[[tracker]] <- all_spectra[[i]]$raw_plots[[j]]$resid$data
    residual_heatmap[[tracker]]$model <- model_key[i, 2]
    residual_heatmap[[tracker]]$quadrant <- quad_key[j, 2]
    tracker <- tracker + 1
  }
}

all_spectra_dat <- dplyr::bind_rows(all_spectra_dat)
residuals <- dplyr::bind_rows(residuals)
residual_heatmap <- dplyr::bind_rows(residual_heatmap)



# function to make plots for each model
make_spectra_plot <- function(spectra_dat, residuals_dat, residual_heatmap_dat){
  rd <- spectra1$raw_plots[[1]]$comp$data[which(spectra1$raw_plots[[1]]$comp$data$source == "data"),]
  
  
  max_spec <- log10(max(c(spectra_dat$N, rd$N), na.rm = T))
  min_spec <- log10(min(c(spectra_dat$N, rd$N), na.rm = T))
  
  spectra_dat$quadrant <- factor(spectra_dat$quadrant, levels = c("TL", "TR", "BL", "BR"))
  residual_heatmap_dat$quadrant <- factor(residual_heatmap_dat$quadrant, levels = c("TL", "TR", "BL", "BR"))
  residuals_dat$quadrant <- factor(residuals_dat$quadrant, levels = c("TL", "TR", "BL", "BR"))
  residuals_dat$quadrant <- factor(residuals_dat$quadrant, levels = c("TL", "TR", "BL", "BR"))

  
  # mod spectra
  ms <- ggplot(spectra_dat[spectra_dat$source == "model",], 
               aes(x = p1, y = p2, fill = log10(N))) + 
    geom_tile() +
    facet_wrap(~quadrant, nrow = 2) +
    theme_bw() +
    scale_color_viridis_c(na.value = "white", option = "inferno", limits = c(min_spec, max_spec)) +
    scale_fill_viridis_c(na.value = "white", option = "inferno", limits = c(min_spec, max_spec)) +
    xlab("Derived Allele Count (HAW)") + ylab("Derived Allele Count (NAM)") +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0,10, by = 2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0,100, by = 20)) +
    ggtitle("Model") +
    theme(strip.background = element_blank())
  
  ms_leg <- cowplot::get_legend(ms)
  
  ms <- ms + theme(legend.position = "none")
  
  # real spectra
  rs <- ggplot(rd, 
               aes(x = p1, y = p2, fill = log10(N))) + 
    geom_tile() +
    theme_bw() +
    scale_color_viridis_c(na.value = "white", option = "inferno", limits = c(min_spec, max_spec)) +
    scale_fill_viridis_c(na.value = "white", option = "inferno", limits = c(min_spec, max_spec)) +
    xlab("Derived Allele Count (HAW)") + ylab("Derived Allele Count (NAM)") +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0,10, by = 2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0,100, by = 20)) +
    theme(strip.background = element_blank(), legend.position = "none") +
    ggtitle("Data")
  
  # residual hist
  his <- ggplot(residuals_dat, aes(x = resid)) +
    geom_histogram(color = "steelblue", fill = "steelblue", bins = 50) + 
    theme_bw() +
    facet_grid(quadrant~model) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    # scale_y_log10(expand = expansion(mult = c(0, .1))) +
    xlab("Residuals") + theme(strip.background = element_blank()) 
  
  # residual plot
  resid_max <- max(abs(residual_heatmap_dat$N), na.rm = T)
  resid <- ggplot(residual_heatmap_dat, 
                  aes(x = p1, y = p2, fill = N)) + 
    geom_tile() +
    facet_wrap(~quadrant, nrow = 2) +
    theme_bw() +
    # scale_color_viridis_c(na.value = "white", option = "magma", direction = -1) +
    # scale_fill_viridis_c(na.value = "white", option = "magma", direction = -1) +
    scico::scale_fill_scico(limits = c(-1*resid_max, resid_max), palette = "vik") +
    scico::scale_color_scico(limits = c(-1*resid_max, resid_max), palette = "vik") +
    xlab("Derived Allele Count (HAW)") + ylab("Derived Allele Count (NAM)") +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0,10, by = 2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0,100, by = 20)) +
    theme(strip.background = element_blank())
  
  
  # arrange
  gridExtra::grid.arrange(rs, ms, ms_leg,
                          layout_matrix = matrix(c(1,1,2,2,3,3), nrow = 2, ncol = 3),
                          widths = c(1,1,.1))
  
  return(list(rs = rs, ms = ms, ms_leg = ms_leg, resid = resid, histogram = his))
}


# call for each model



pdf("plots/Figure_4.pdf", width = 15, height = 8.5)
f4 <- make_spectra_plot(all_spectra_dat[which(all_spectra_dat$model == "Three Epoch"),],
                        residuals[which(residuals$model == "Three Epoch"),],
                        residual_heatmap[which(residual_heatmap$model == "Three Epoch"),])
dev.off();dev.off();
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_4.jpg plots/Figure_4.pdf")


pdf("plots/Figure_S6.pdf", width = 15, height = 8.5)
fS6 <- make_spectra_plot(all_spectra_dat[which(all_spectra_dat$model == "Found and Grow"),],
                        residuals[which(residuals$model == "Found and Grow"),],
                        residual_heatmap[which(residual_heatmap$model == "Found and Grow"),])
dev.off();dev.off()
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S6.jpg plots/Figure_S6.pdf")


pdf("plots/Figure_S7.pdf", width = 15, height = 8.5)
fS7<- make_spectra_plot(all_spectra_dat[which(all_spectra_dat$model == "Two Epoch"),],
                         residuals[which(residuals$model == "Two Epoch"),],
                         residual_heatmap[which(residual_heatmap$model == "Two Epoch"),])
dev.off();dev.off()
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S7.jpg plots/Figure_S7.pdf")


pdf("plots/Figure_S8.pdf", width = 15, height = 8.5)
fS8 <- make_spectra_plot(all_spectra_dat[which(all_spectra_dat$model == "Zhan"),],
                         residuals[which(residuals$model == "Zhan"),],
                         residual_heatmap[which(residual_heatmap$model == "Zhan"),])
dev.off();dev.off()
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S8.jpg plots/Figure_S8.pdf")


pdf("plots/Figure_S9.pdf", width = 11, height = 8.5)
gridExtra::grid.arrange(f4$resid + ggtitle("Three Epoch"), 
                        fS6$resid + ggtitle("Found and Grow"),
                        fS7$resid + ggtitle("Two Epoch"),
                        fS8$resid + ggtitle("Zhan"),
                        ncol = 2)
dev.off();dev.off()
shell("C://usr/bin/gswin64c.exe -sDEVICE=jpeg -r288 -o plots/Figure_S9.jpg plots/Figure_S9.pdf")
