library(ggplot2)
library(ggfortify)
library(lme4)
library(lmerTest)
library(tidyverse)
library(car)
library(emmeans)
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

setwd('~/Box/skinner/projects_analyses/Project SAPAP3/ca')
load('MiceOperantsData.rdata')

#######
# preliminary behavioral analyses of calcium reversal SAPAP3 data
# timepoint -- time within run
# 

df <- as.tibble(MiceOperantsData$mice_bin_df) %>% mutate(mouseID = as.factor(mouseID), mag_entry = magzineentry, time_sc = scale(timepoint), 
                                                         time_sc_neg_inv = scale(-1/timepoint), treatment = as.factor(treatment))
# learning curves by response type
# setwd('~/Box/skinner/projects_analyses/Project SAPAP3/ca/plots')
setwd('~/OneDrive/papers/sapap3_ca/plots/')
m1 <- glmer(action ~ type*time_sc + (1|mouseID),family = binomial(link = "logit"),data = df)
summary(m1)
Anova(m1, '3')
em1 <- as_tibble(emmeans(m1, specs = c("type", "time_sc"), at = list(time_sc = c(-2,0,2))))
em1$`Likelihood of responding` <- em1$emmean
pdf('ca_rev_m1.pdf', width = 8, height = 6)
ggplot(em1, aes(time_sc, `Likelihood of responding`, group = type, color = type)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) #+ theme_dark() + scale_color_viridis_d(option = "inferno")
dev.off()

# hyperbolic time inferior to linear time -- can treat as approximately linear
# m1a <- glmer(action ~ type*time_sc_neg_inv + (1|mouseID),family = binomial(link = "logit"),data = df)
# summary(m1a)
# Anova(m1a, '3')
# anova(m1, m1a)



m2 <- glmer(action ~ type*time_sc*treatment + (1|mouseID),family = binomial(link = "logit"),data = df)
summary(m2)
Anova(m2, '3')
vif(m2)
anova(m1,m2)
em2 <- as_tibble(emmeans(m2, specs = c("type", "time_sc", "treatment"), at = list(time_sc = c(-2,0,2))))
em2$`Likelihood of responding` <- em2$emmean
pdf('ca_rev_m2.pdf', width = 8, height = 6)
ggplot(em2, aes(time_sc, `Likelihood of responding`, color = treatment, group = treatment)) + geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
  facet_wrap(~type)#+ theme_dark() + scale_color_viridis_d(option = "inferno")
dev.off()
# cache results
save('ca_models.RData')
