library(ggplot2)
library(ggfortify)
library(lme4)
library(lmerTest)
library(tidyverse)

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

setwd('~/Box/skinner/projects_analyses/Project SAPAP3/ca')
load('MiceOperantsData.rdata')

df <- as.tibble(MiceOperantsData$mice_bin_df) %>% mutate(mouseID = as.factor(mouseID), mag_entry = magzineentry, time_sc = scale(timepoint), 
                                                         time_sc_neg_inv = scale(-1/timepoint), treatment = as.factor(treatment))
# learning curves by response type
# setwd('~/Box/skinner/projects_analyses/Project SAPAP3/ca/plots')
setwd('~/OneDrive/papers/sapap3_ca/plots/')

knots = 1

r <- ggplot(df, aes(timepoint, action, color = type, lty = treatment)) + binomial_smooth(formula = y ~ splines::ns(x, knots)) 
m <- ggplot(df, aes(timepoint, mag_entry, lty = treatment)) + binomial_smooth(formula = y ~ splines::ns(x, knots))
pdf('ca_beh.pdf', width = 12, height = 6)
ggarrange(r,m, ncol = 2)
dev.off()

ri <- ggplot(df, aes(timepoint, action, color = type)) + binomial_smooth(formula = y ~ splines::ns(x, knots)) + facet_grid(cols = vars(mouseID), rows = NULL)
mi <- ggplot(df, aes(timepoint, mag_entry, color = type)) + binomial_smooth(formula = y ~ splines::ns(x, knots))+ facet_grid(cols = vars(mouseID), rows = NULL)
pdf('ca_ind_beh.pdf', width = 20, height = 10)
ggarrange(ri,mi, ncol = 1, nrow = 2)
dev.off()

p1 <- ggplot(df, aes(timepoint, action, color = type, lty = treatment,  group = mouseID)) + 
  binomial_smooth(formula = y ~ splines::ns(x, knots)) + facet_grid(~type)
p2 <- ggplot(df, aes(timepoint, mag_entry, lty = treatment,  group = mouseID)) + 
  binomial_smooth(formula = y ~ splines::ns(x, knots))

# by mouse
pdf('ca_ind_learning_curves_mag_entry.pdf', width = 10, height = 10)
ggarrange(p1,p2, ncol = 1, nrow = 2)
dev.off()



# by mouse
pdf('ca_ind_magazine.pdf', width = 20, height = 20)
ggplot(df, aes(timepoint, mag_entry, color = mouseID)) + binomial_smooth(formula = y ~ splines::ns(x, 3))
dev.off()


