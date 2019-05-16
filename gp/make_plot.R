library(ggplot2)

# read in and add column names to data
dat <- read.table("gp/BayesB_leave_one_out_cv.txt", header = F)
colnames(dat) <- c("sampleID", "samplenum", "observed", "predicted", "h2")


# make the plot and a basic linear model
tlm <- lm(observed~predicted, dat)
lm_summary <- summary(tlm)
lm_anova <- anova(tlm)
r2 <- lm_summary$r.squared

plot1 <-  ggplot(dat, aes(x = observed, y = predicted, color = h2)) + geom_point() + geom_smooth(method = "lm") +
  theme_bw() + scale_color_viridis_c() + 
  labs(color = expression("h"^{2}), x = "observed BV", y = "predicted BV")

pdf(file = "gp/gp_leave_one_out.pdf", plot1)
plot1
dev.off(); dev.off()
