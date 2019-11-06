source("scripts/import_dadi_results.R"); source("scripts/interpret_dadi_units.R")
library(dplyr); library(ggplot2)
mu <- 8.4e-9
g <- .3
ratio <- 9370/302446 # ratio of included snps
L <-  1373747*ratio # approx number of considered bases.


# pass1 <- import_dadi_results(c("data/dadi_inputs/cat_NH_1st_pass_dportik.txt", "data/dadi_inputs/cat_hg_r1_fix_combined.txt"),
#                              mu = mu, g = g, L = L)
# pass2 <- import_dadi_results(c("data/dadi_inputs/cat_NH_best_r2.txt", "data/dadi_inputs/cat_NH_3e_re_best_r2.txt"),
#                              mu = mu, g = g, L = L)
# pass3 <- import_dadi_results("data/dadi_inputs/cat_NH_best_r3.txt",
#                              mu = mu, g = g, L = L)
# pass4 <- import_dadi_results("data/dadi_inputs/cat_NH_best_r4.txt", mu = mu, g = g, L = L)


pass1 <- import_dadi_results(c("data/dadi_inputs/cat_NH_1st_pass_dportik.txt", "data/dadi_inputs/cat_NH_hg_r1.txt"),
                             mu = mu, g = g, L = L)
pass2 <- import_dadi_results(c("data/dadi_inputs/cat_NH_2nd_pass_out_dportik.txt", "data/dadi_inputs/cat_NH_hg_r2_fix_weighted.txt"),
                             mu = mu, g = g, L = L)
pass3 <- import_dadi_results(c("data/dadi_inputs/cat_NH_r3.out", "data/dadi_inputs/cat_NH_hg_r3.txt"),
                             mu = mu, g = g, L = L)
pass4 <- import_dadi_results("data/dadi_inputs/cat_NH_all_r4.txt", mu = mu, g = g, L = L)


# AIC dist per model
rdf <- rbind(cbind(pass1$rdf, pass = 1),
             cbind(pass2$rdf, pass = 2),
             cbind(pass3$rdf, pass = 3),
             cbind(pass4$rdf, pass = 4))
rdf$pass <- as.factor(rdf$pass)

rdf$model[rdf$model == "founder_asym_growth_both"] <- "Found and Grow"
rdf$model[rdf$model == "founder_asym_hist_3epoch_exp_growth_p1"] <- "Three Epoch Found and Grow"
rdf$model[rdf$model == "founder_asym_hist_igrowth_p2"] <- "Zhan Model"

## only the "good" models
mins <- tapply(rdf$AIC, rdf$model, min)
goods <- names(mins)[which(mins <= 2500)]

gmodels <- rdf[rdf$model %in% goods,]

pdf("./plots/dadi_model_comparison.pdf")
ggplot(gmodels,aes(x = model, y = log(AIC), color = pass)) + 
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Set1"))
dev.off(); dev.off()

# best overall model
best.reps <- rdf %>% group_by(model) %>% group_by(pops, add = TRUE) %>% top_n(-1, AIC)
best.reps <- arrange(best.reps, pops, model)
View(pass1$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW)

# plots for the best model
# best model across all passes
ilist <- rbind(cbind(pass1$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 1),
               cbind(pass2$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 2),
               cbind(pass3$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 3),
               cbind(pass4$ilist$founder_asym_hist_3epoch_exp_growth_p1_NAM_HAW, pass = 4))
ilist$pass <- as.factor(ilist$pass)

# best model in pass 4
ilist2 <- rbind(cbind(pass1$ilist$founder_asym_growth_both_NAM_HAW, pass = 1),
                cbind(pass2$ilist$founder_asym_growth_both_NAM_HAW, pass = 2),
                cbind(pass3$ilist$founder_asym_growth_both_NAM_HAW, pass = 3),
                cbind(pass4$ilist$founder_asym_growth_both_NAM_HAW, pass = 4))

ilist2$pass <- as.factor(ilist2$pass)

# combine data
## merge
ilist3.1 <- ilist[,c("model", "AIC", "nuA", "Tg", "Ts", "Tg2", "Tg3", "s", "nu1F", "nu2F", "nuG2", "m12", "m21", "pass")]
ilist3.1$AncientTime <- ilist3.1$Tg + ilist3.1$Ts + ilist3.1$Tg2 + ilist3.1$Tg3
ilist3.1$SplitTime <- ilist3.1$Ts + ilist3.1$Tg3
colnames(ilist3.1)[c(9,10,11)] <- c("Ne_NA", "Ne_Ha", "Ne_Split")
ilist3.2 <- ilist2[,c("model", "AIC", "nuA", "nu1", "nu2", "m12", "m21", "Ti", "s", "pass")]
colnames(ilist3.2)[c(4,5,8)] <- c("Ne_NA", "Ne_Ha", "SplitTime")
ilist3.2$Ne_Split <- ilist3.2$nuA
ilist3 <- merge(ilist3.1, ilist3.2, all = T)
ilist3$model <- ifelse(ilist3$model == "founder_asym_growth_both", "Found and Grow", "Three Epoch Found and Grow")


## function to grab legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}

## grab legend
dpb.time <- ggplot(ilist3, aes(y = log10(SplitTime), x = log10(s), color = AIC, shape = pass)) + 
  geom_point() + geom_jitter() + theme_bw() + scale_color_viridis_c() + facet_wrap(~model) +
  ylab("log10(Time since Hawaii establishment)") + xlab("log10(Ne Hawaii Establishment)")

legend <- g_legend(dpb.time)

## make plots
dpb.time <- ggplot(ilist3, aes(y = log10(SplitTime), x = log10(s), color = AIC, shape = pass)) + 
  geom_point() + theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "white")) + 
  scale_color_viridis_c() + facet_wrap(~model) +
  # geom_point(data = ilist3[ilist3$best == 1,], aes(y = log10(SplitTime), x = log10(Ne_NA)), color = "red") +
  ylab(bquote(~log[10]~ '(Years Since Split)')) + xlab(bquote('Establishment' ~log[10]~ '(' ~N[e]~ '), Hawaii')) + 
  scale_y_continuous(labels = function(x) sprintf("%.5f", x))


dpb.end.size <- ggplot(ilist3, aes(x = log10(Ne_NA), y = log10(Ne_Ha), color = AIC, shape = pass)) + 
  geom_point() + theme_bw()+ 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_color_viridis_c() + facet_wrap(~model) +
  # geom_point(data = ilist3[ilist3$best == 1,], aes(x = log10(s), y = log10(Ne_Ha)), color = "red") +
  xlab(bquote('Current' ~log[10]~ '(' ~N[e]~ '), North America')) + ylab(bquote('Current' ~log[10]~ '(' ~N[e]~ '), Hawaii')) + 
  scale_y_continuous(labels = function(x) sprintf("%.5f", x))

dpb.mig <- ggplot(ilist3, aes(x = m12, y = m21, color = AIC, shape = pass)) + 
  geom_point() + theme_bw()+ 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +scale_color_viridis_c() + facet_wrap(~model) +
  # geom_point(data = ilist3[ilist3$best == 1,], aes(x = m12, y = m21), color = "red") +
  xlab("Migration: Hawaii -> NA") + ylab("Migration: NA -> Hawaii") +
  scale_y_continuous(labels = function(x) sprintf("%.5f", x))

## combine plots
pdf("plots/dadi_summary_plot.pdf")
gridExtra::grid.arrange(dpb.time, dpb.end.size, dpb.mig, legend,
                        layout_matrix = matrix(c(1,2,3,5,5,5), nrow = 3, ncol = 2),
                        widths = c(1,.1))
dev.off();dev.off();

#============try comparing the 3 epoch to the zahn model==========
ilist4 <- rbind(cbind(pass1$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 1),
               cbind(pass2$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 2),
               cbind(pass3$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 3),
               cbind(pass4$ilist$founder_asym_hist_igrowth_p2_NAM_HAW, pass = 4))
ilist4$pass <- as.factor(ilist4$pass)

ilist4.1 <- ilist4
ilist4.1$SplitTime <- ilist4.1$Ts + ilist4.1$Tg2
ilist4.1$AncientTime <- ilist4.1$Tg + ilist4.1$Tg2 + ilist4.1$Ts
ilist4.1$Ne_NA <- ilist4.1$nuG - ilist4.1$s
ilist4.1$Ne_Split <- ilist4.1$nuG
ilist4.1$Ne_Ha <- ilist4.1$nu2F
ilist4.1$AncientGrowthDelta <- ilist4.1$nuG - ilist4.1$nuA
ilist1.1 <- ilist
ilist1.1$AncientGrowthDelta <- ilist1.1$nuG - ilist1.1$nuA
ilist1.1 <- cbind(ilist1.1, ilist3.1[,c("AncientTime", "SplitTime", "Ne_NA", "Ne_Split", "Ne_Ha")])
ilist5 <- merge(ilist4.1, ilist1.1, all = T)
ilist5[ilist5$model == "founder_asym_hist_3epoch_exp_growth_p1",]$model <- "Three Epoch Found and Grow"
ilist5[ilist5$model == "founder_asym_hist_igrowth_p2",]$model <- "Zhan"


zc.time <- ggplot(ilist5, aes(y = log10(SplitTime), x = log10(s), color = AIC, shape = pass)) + 
  geom_point() + theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "white")) + 
  scale_color_viridis_c() + facet_wrap(~model) +
  # geom_point(data = ilist3[ilist3$best == 1,], aes(y = log10(SplitTime), x = log10(Ne_NA)), color = "red") +
  ylab(bquote(~log[10]~ '(Time Since Split)')) + xlab(bquote('Establishment' ~log[10]~ '(' ~N[e]~ '), Hawaii')) + 
  scale_y_continuous(labels = function(x) sprintf("%.5f", x))

zc.end.size <- ggplot(ilist5, aes(x = log10(Ne_NA), y = log10(Ne_Ha), color = AIC, shape = pass)) + 
  geom_point() + theme_bw()+ 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_color_viridis_c() + facet_wrap(~model) +
  xlab("log10(North America Ne)") + ylab("log10(Hawaii Ne)") +
  scale_y_continuous(labels = function(x) sprintf("%.5f", x))

zc.mig <- ggplot(ilist5, aes(x = m12, y = m21, color = AIC, shape = pass)) + 
  geom_point() + theme_bw()+ 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +scale_color_viridis_c() + facet_wrap(~model) +
  xlab("Migration: Hawaii -> NA") + ylab("Migration: NA -> Hawaii") +
  scale_y_continuous(labels = function(x) sprintf("%.5f", x))

pdf("plots/zhan_compare_supplment.pdf")
gridExtra::grid.arrange(zc.time, zc.end.size, zc.mig, legend,
                        layout_matrix = matrix(c(1,2,3,5,5,5), nrow = 3, ncol = 2),
                        widths = c(1,.1))
dev.off();dev.off()
