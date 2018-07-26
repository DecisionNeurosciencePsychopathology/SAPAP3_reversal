## GLME on binned response data from Manning et al (2018) 
#prepare data as follows:
#timestamp data as excel spreadsheet (e.g. "SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx") with three columns, A) subject ID, B) Correct lever press timestamp and C) Incorrect lever press timestamp. where there are more values in B than C just have zeros
#Blind details as excel spreadsheet (e.g. "Reversal cFos cohort blind.xlsx") with two columns, A) subject ID ("Phsyical Tag") and B)Genotype
#cFos cell counts as csv file (e.g. "cfos_040618_finalmice.xlsx") with columns including subject ID and cFos measurments arranged/labelled by region of interest

library(readr)
library(lme4)
library(ggplot2)
#library(dplyr)
library(tidyr)
library(psych)
library(gdata)
library(R.matlab)
library(xtable)
library(Hmisc)
library(nnet)
library(reshape2)
library(ggbiplot)
library(corrplot)
library(lsmeans)
library(factoextra)
library(ggfortify)
library(readxl)
library(linbin)
library(MASS)
library(lattice)
library(fastICA)
library(plotly)
library(stargazer)
library(rio)

# read in the data
setwd("~/code/SAPAP3_reversal/")
d <- read_excel("~/code/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
View(d)


# read in genotype
g <- read_excel("~/code/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
g$id <- g$`Physical Tag`
g <- g[,2:3]


# read in cfos cell counts

cfos <- read_csv("~/code/SAPAP3_reversal/cFos_040618_finalmice.csv")
names(cfos)[names(cfos)=="ID"] <- "id"
names(cfos)[names(cfos)=="correct"] <- "tot_correct"
names(cfos)[names(cfos)=="incorrrect"] <- "tot_incorrect"
names(cfos)[names(cfos)=="Genotype"] <- "genotype01"
names(cfos)[names(cfos)=="NAcS"] <- "NAcS"
names(cfos)[names(cfos)=="NAcC"] <- "NAcC"


View(cfos)
# check PCA on cfos

just_rois <- cfos[,7:16]
#head(just_rois)
feedROIcor <- cor(just_rois)
pdf("cfos correlations by regions_no.pdf", width=12, height=12)
corrplot(feedROIcor, cl.lim=c(0,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-feedROIcor, sig.level=0.5, insig = "blank")
dev.off()


cfos.pca = prcomp((just_rois),scale = TRUE, center = TRUE)

# run PCA and write component scores
cfos_pcas <- get_pca_ind(cfos.pca)
cfos$val1 <- cfos_pcas$coord[,1]
# hdf$val2 <- cfos_pcas$coord[,2]
# hdf$val3 <- cfos_pcas$coord[,3]
# hdf$val4 <- cfos_pcas$coord[,4]

# look
summary(cfos.pca)

# check the variance explained by various factors
plot(cfos.pca,type = 'l')

# arrange data for analysis
# make tall cfos
cshort <- cfos[,c(1:2,7:16)]
cshort <- as.data.frame(cshort)
t_cfos <- melt(cshort,id.vars = c("id","genotype01"))
View(t_cfos)
names(t_cfos)[names(t_cfos)=="value"] <- "cfos"
names(t_cfos)[names(t_cfos)=="variable"] <- "region"

#test effect of genotype on regional cfos and plot
summary(cm1 <- lm(cfos ~ region*genotype01, data = t_cfos))
anova(cm1)

pdf("cfos by region and genotype.pdf", width=12, height=6)
boxplot(cfos ~ genotype01 + region, data = t_cfos, main = "cfos activity by SAPAP3 genotype (blue = WT, purple = KO) and region",
        xlab = "Region", ylab = "cfos", varwidth = TRUE, col =  cm.colors(2))
dev.off()

# arrange timestamps in bins
# binsize is in 0.1 seconds (from MEDPC)
d$id <- d$`subject ID`
d$corr <- d$`Correct lever press`
d$inc <- d$`Incorrect lever press`
d <- d[,4:6]
d[d==0] <- NA
View(d)
ids <- unique(d$id)
describe(d$corr)
describe(d$inc)
# end <- 18000
start <- 0
binsize <- 40
maxbin <- 16800

# censor data for mice that finished in <30 minutes
c <- data.frame()
for (id in ids)
## insert the if condition that sets a different end for some ids, or just set end to equal the time of last response
  {end <- max(c(max(d$inc[d$id==id]),max(na.omit(d$corr[d$id==id]))))
  # make sure we are not removing idle time at the end of session for mice with full 30 minutes testing
  if (end>maxbin) {end <- 18000}
  bins <- seq(start,end,binsize)
  t=cut(d$corr[d$id==id],c(bins), labels = FALSE)
b <- table( factor(t, levels = 1:length(bins)-1))
b <- as.data.frame(b)
b$id <- id
names(b)[names(b)=="Var1"] <- "time"
names(b)[names(b)=="Freq"] <- "corr"
tinc=cut(d$inc[d$id==id],c(bins), labels = FALSE)
i <- table( factor(tinc, levels = 1:length(bins)-1))
i <- as.data.frame(i)
b$inc <- i$Freq
c <- rbind(c,b)
}

c <- c[c$time!="0",]

c$t.ord <- as.ordered(c$time)
c$t.num <- as.numeric(c$time)

# look at data
plot(c$time,c$inc)
plot(c$time,c$corr)

hist(c$inc)
hist(c$corr)
length(c$inc)

ggplot(aes(x=t.num, y=corr),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Correct responses") +
  stat_smooth(data=c, aes(x=t.num, y=corr), size=0.4, alpha=0.1, se=FALSE, method="loess")

ggplot( aes(x=t.num, y=inc),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Incorrect responses") +
  stat_smooth(data=c, aes(x=t.num, y=inc), size=0.4, alpha=0.1, se=FALSE, method="loess")

# set up negative binomial models
# estimate overdispersion parameter theta
theta.corr <- theta.ml(c$corr, mean(c$corr), 513, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)
theta.inc <- theta.ml(c$inc, mean(c$inc), 513, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)


# for correct presses
#summary(mcorr1 <- glm(corr ~ t.ord + (1:id), family = negative.binomial(theta = theta.corr), data = c))
summary(mcorr2 <- glm(corr ~ t.num + (1:id), family = negative.binomial(theta = theta.corr), data = c))

#anova(mcorr1,mcorr2)

# for incorrect presses
summary(minc1 <- glm(inc ~ t.num + (1:id), family = negative.binomial(theta = theta.corr), data = c))

# combine all responses into one variable with correct/incorrect as factor
c.all <- melt(c,id.vars = c("id","time","t.num","t.ord"),value.name = "response", variable.name = "type")
names(c.all)[names(c.all)=="value"] <- "response"
names(c.all)[names(c.all)=="variable"] <- "type"

View(c.all)

nbfit <- suppressWarnings(fitdistr(c.all$response, "negative binomial"))
nbfit$estimate
logLik(nbfit)

c.all$response <- as.numeric(c.all$response)

# plot log(responses) over time
pdf("Response curves gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = c.all) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=c.all, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter")
dev.off()
ggplot(aes(x=t.num, y=response, color = type),data = c.all) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses") +
  geom_count(position = "jitter")

# merge with genotype
# corr/inc separated
hdf <- merge(c,g) # horizontal data frame
#View(hdf)

# plot responses (individual points)
ggplot(aes(x=t.num, y=corr, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess", size = 1) + theme_bw(base_size=20) +
  xlab("Time, 4s bins") + ylab("Correct responses") +
  geom_count(position = "jitter")

ggplot( aes(x=t.num, y=inc, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess",  size = 1) + theme_bw(base_size=20) +
  xlab("Time, 4s bins") + ylab("Incorrect responses") +
  geom_count(position = "jitter")


# all together
bdf <- merge(c.all,g)
#View(bdf)

names(cfos)[names(cfos)=="grooming time"] <- "grooming.time"

# merge with cell counts
hdfc <- merge(hdf,cfos)
bdfc <- merge(bdf,cfos)

# differnt ways to graph model of responses across the session
pdf("Log response curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = bdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter") +
  facet_grid(Genotype ~ .)
dev.off()


pdf("Log response curves by genotype gam combined responses.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = Genotype),data = bdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = Genotype), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()

pdf("Log response curves by genotype gam incorrect responses.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(inc + .1), color = Genotype),data = hdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=hdf, aes(x=t.num, y=log(inc + .1), color = Genotype), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()


pdf("Log response curves by genotype gam correct responses.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(corr + .1), color = Genotype),data = hdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=hdf, aes(x=t.num, y=log(corr + .1), color = Genotype), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()


#  quick look at genotype
# correct
summary(mcorrg1 <- glm(corr ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.corr), data = hdf))
#incorrect
summary(mincg1 <- glm(inc ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.inc), data = hdf))

# both
# calculate theta for both responses
theta.resp <- theta.ml(na.omit(c.all$response), mean(na.omit(c.all$response)), 972, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)
#summary(mrespg1 <- glm(response ~ 1 + t.num*type + type*Genotype +  (1:id), family = negative.binomial(theta = theta.resp), data = bdf))
#summary(mrespg1 <- glm(response ~ 1 + t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.resp), data = bdf))
summary(mrespg1 <- glm(response ~ t.num*type*Genotype +  (1:id), family = negative.binomial(theta = theta.resp), data = bdf))
car::Anova(mrespg1)

#test cFos ROI (alone) model e.g. PrL
#swap in other regions of interest to check them in the model
#consider correcting p value if running multiple models for cfos in multiple regions
summary(mrespg2 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + PrL*t.num*Genotype + PrL*type*Genotype + PrL*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg2)

#quick graphs to summarize interactions
#for "list" use 3 values of cFos measurement ~mean, mean + 1 SD, mean - 1 SD
hist(cfos$PrL)
mean(cfos$PrL)
sd(cfos$PrL)
lsmip(mrespg2, PrL ~ type | Genotype, at = list(PrL = c(200,400,600)), ylab = "log(response rate)", xlab = "type ", type = "predicted" )
lsmip(mrespg2, PrL ~ Genotype | type, at = list(PrL = c(200,400,600)), ylab = "log(response rate)", xlab = "type ", type = "predicted" )
lsmip(mrespg2, PrL ~ t.num | type, at = list(PrL = c(200,400,600),t.num  = c(1,length(bins)/2, length(bins))), ylab = "log(response rate)", xlab = "time ", type = "predicted" )

#PrL GRAPH 
summary(mrespg2 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + PrL*t.num*Genotype + PrL*type*Genotype + PrL*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg2)

leastsquare = lsmeans::lsmeans(mrespg2, pairwise ~ PrL:Genotype:type, at = list(PrL = c(200,400,600)), adjust='tukey')
library(dplyr)

CLD = cld(leastsquare, alpha=0.05, Letters=letters,
          adjust='tukey')
CLD$type <- recode_factor(CLD$type,corr="Correct",inc="Incorrect")
pdf(file = "PrL plot publish.pdf", width = 10, height = 6)
pd = position_dodge(15)    ### How much to jitter the points on the plot
ggplot(CLD, aes( x = PrL, y = lsmean, color = Genotype, label = .group)) +
  facet_wrap( ~ type) +
  
  geom_point(shape  = 16,
             size   = 4,
             position = pd) +
  
  geom_errorbar(
    aes(ymin  =  asymp.LCL,
        ymax  =  asymp.UCL),
    width =  0.2,
    size  =  0.7,
    position = pd
  ) +  theme_bw() +  theme(
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  ) +  ylab("Log-probability of response") + xlab("Prelimbic cortex cfos") +
  ggtitle ("Behavior in reversal phase by genotype, regional cfos level, and response type") +
  labs( caption  = paste0(
    "\n",
    "Boxes indicate the LS mean.\n",
    "Error bars indicate the 95% ",
    "confidence interval of the LS mean, Sidak method for 12 estimates. \n",
    "Means sharing a letter are ",
    "not significantly different ",
    "(Tukey-adjusted comparisons for 12 estimates)."),
    hjust = 0.5 ) +
  geom_text(nudge_x = c(20, 10, 20, 10, 10, 10,10, 10, 10, 10, 10, 10),
            nudge_y = c(0,  0, 0,  0, 0 , 0,0,0, .05,  -.05, 0, 0), color   = "black") +
  # nudge_y = 0, color   = "black") +
  
  scale_color_manual(values = c("magenta3","black"))
dev.off()
