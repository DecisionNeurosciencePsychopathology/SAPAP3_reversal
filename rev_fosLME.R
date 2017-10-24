## LME on binned response data from Lizzie Manning and Suzanne Ahmari's free operant deterministic reversal experiment

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
# read in the data
# setwd("/Users/manninge/Documents/GitHub/SAPAP3_reversal/")
# d <- read_excel("/Users/manninge/Documents/GitHub/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
setwd("~/code/SAPAP3_reversal/")
d <- read_excel("~/code/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
View(d)


# read in genotype
# g <- read_excel("/Users/manninge/Documents/GitHub/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
g <- read_excel("~/code/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
g$id <- g$`Physical Tag`
g <- g[,2:3]


# read in cfos cell counts
# cfos <- read_csv("/Users/manninge/Documents/GitHub/SAPAP3_reversal/SAPAP3 reversal cFos density_Final mice.csv")
cfos <- read_csv("~/code/SAPAP3_reversal/SAPAP3 reversal cFos density_Final mice.csv")
names(cfos)[names(cfos)=="ID"] <- "id"
names(cfos)[names(cfos)=="correct"] <- "tot_correct"
names(cfos)[names(cfos)=="incorrrect"] <- "tot_incorrect"
names(cfos)[names(cfos)=="Genotype"] <- "genotype01"
names(cfos)[names(cfos)=="NAc S"] <- "NAccS"
names(cfos)[names(cfos)=="NAc C"] <- "NAccC"


View(cfos)
# check PCA on cfos

just_rois <- cfos[,7:18]
#head(just_rois)
feedROIcor <- cor(just_rois)
pdf("cfos correlations by regions.pdf", width=12, height=12)
corrplot(feedROIcor, cl.lim=c(0,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-feedROIcor, sig.level=0.5, insig = "blank")
dev.off()

cfos.pca = prcomp((just_rois),scale = TRUE, center = TRUE)

# PCA for Dorsal striatum
DS_rois <- cfos[,13:15]
DS_cfos.pca = prcomp((DS_rois),scale = TRUE, center = TRUE)
DS_cfos_pcas <- get_pca_ind(DS_cfos.pca)
cfos$val1 <- DS_cfos_pcas$coord[,1]
summary(DS_cfos.pca)
plot(DS_cfos.pca,type = 'l')

#PCA for DLS/CMS
DCS_rois <- cfos[,c(13,15)]
DCS_cfos.pca = prcomp((DCS_rois),scale = TRUE, center = TRUE)
DCS_cfos_pcas <- get_pca_ind(DCS_cfos.pca)
cfos$val7 <- DCS_cfos_pcas$coord[,1]

# PCA for Ventral striatum
VS_rois <- cfos[,16:18]
VS_cfos.pca = prcomp((VS_rois),scale = TRUE, center = TRUE)
VS_cfos_pcas <- get_pca_ind(VS_cfos.pca)
cfos$val2 <- VS_cfos_pcas$coord[,1]

# PCA for Nucleus accumbens
NAc_rois <- cfos[,17:18]
NAc_cfos.pca = prcomp((NAc_rois),scale = TRUE, center = TRUE)
NAc_cfos_pcas <- get_pca_ind(NAc_cfos.pca)
cfos$val3 <- NAc_cfos_pcas$coord[,1]

# PCA for mPFC
mPFC_rois <- cfos[,c(7,11)]
mPFC_cfos.pca = prcomp((mPFC_rois),scale = TRUE, center = TRUE)
mPFC_cfos_pcas <- get_pca_ind(mPFC_cfos.pca)
cfos$valm6 <- mPFC_cfos_pcas$coord[,1]

# PCA for OFC
OFC_rois <- cfos[,c(7:8)]
OFC_cfos.pca = prcomp((OFC_rois),scale = TRUE, center = TRUE)
OFC_cfos_pcas <- get_pca_ind(OFC_cfos.pca)
cfos$val5 <- OFC_cfos_pcas$coord[,1]

# PCA for PFC
PFC_rois <- cfos[,c(7:8,11:12)]
PFC_cfos.pca = prcomp((PFC_rois),scale = TRUE, center = TRUE)
PFC_cfos_pcas <- get_pca_ind(PFC_cfos.pca)
cfos$valm8 <- PFC_cfos_pcas$coord[,1]

# PCA for PrL-NACs CIRCUIT
PrLNAcS_rois <- cfos[,c(5,12)]
PrLNAcS_cfos.pca = prcomp((PrLNAcS_rois),scale = TRUE, center = TRUE)
PrLNAcS_cfos_pcas <- get_pca_ind(PrLNAcS_cfos.pca)
cfos$valm9 <- PrLNAcS_cfos_pcas$coord[,1]

#Re-merge cell counts with other data

# run PCA (code from Alex 081717)
#ds.pca = prcomp((just_rois),scale = TRUE, center = TRUE)
# and write component scores
#ds_pcas <- get_pca_ind(ds.pca)
#cfos$val1 <- cfos_pcas$coord[,1]


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


## Quick attempt at ICA -- probably over-analyzing these data
# a <- fastICA(just_rois, 5,  alg.typ = "parallel", fun = "logcosh", alpha = 1,
#              method = "C", row.norm = FALSE, maxit = 200,
#              tol = 0.0001, verbose = TRUE)
# par(mfrow = c(1, 3))
# plot(a$X, main = "Pre-processed data")
# plot(a$X %*% a$K, main = "PCA components")
# plot(a$S, main = "ICA components")
# icaloads <- as.data.frame(a$A)
# names(icaloads) <- names(just_rois)
# heatmap(a$A, labCol = names(just_rois))

cfos$genotype01 <- as.factor(cfos$genotype01)

ggplot(aes(x=genotype01, y=IL),data = cfos) +
  geom_count()
ggplot(aes(x=genotype01, y=DMS),data = cfos) +
  geom_count()
ggplot(aes(x=genotype01, y=NAccS),data = cfos) +
  geom_count()


# boxplot(IL ~ genotype01, data = cfos, main = "cfos activity by genotype",
#         xlab = "SAPAP3 genotype", ylab = "cfos, IL", varwidth = TRUE, col =  cm.colors(3))
#
# boxplot(mOFC ~ genotype01, data = cfos, main = "cfos activity by genotype",
#         xlab = "SAPAP3 genotype", ylab = "cfos, mOFC", varwidth = TRUE, col =  cm.colors(3))
#
# boxplot(lOFC ~ genotype01, data = cfos, main = "cfos activity by genotype",
#         xlab = "SAPAP3 genotype", ylab = "cfos, lOFC", varwidth = TRUE, col =  cm.colors(3))

# make tall cfos
cshort <- cfos[,c(1:2,7:18)]
cshort <- as.data.frame(cshort)
t_cfos <- melt(cshort,id.vars = c("id","genotype01"))
names(t_cfos)[names(t_cfos)=="value"] <- "cfos"
names(t_cfos)[names(t_cfos)=="variable"] <- "region"

summary(cm1 <- lm(cfos ~ region*genotype01, data = t_cfos))
anova(cm1)

pdf("cfos by region and genotype.pdf", width=12, height=6)
boxplot(cfos ~ genotype01 + region, data = t_cfos, main = "cfos activity by SAPAP3 genotype (blue = WT, purple = KO) and region",
        xlab = "Region", ylab = "cfos", varwidth = TRUE, col =  cm.colors(2))
dev.off()
# only striatal
# str_rois <- cfos[,13:18]
# head(str_rois)
# feedROIcor <- cor(str_rois)
# corrplot(feedROIcor, cl.lim=c(0,1),
#          method = "circle", tl.cex = 1, tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = 1-feedROIcor, sig.level=0.6, insig = "blank")





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
# binsize <- 1000
# try smaller bins:
binsize <- 40
maxbin <- 16800


c <- data.frame()
for (id in ids)
## Lizzie, here you can insert the if condition that sets a different end for some ids, or just set end to equal the time of last response
  {end <- max(c(max(d$inc[d$id==id]),max(na.omit(d$corr[d$id==id]))))
  # make sure we are not removing idle time at the end of sessions for those 5-6 mice
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

plot(c$time,c$inc)
plot(c$time,c$corr)

hist(c$inc)
hist(c$corr)
length(c$inc)

ggplot(aes(x=t.num, y=corr),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Correct responses") +
  stat_smooth(data=c, aes(x=t.num, y=corr), size=0.4, alpha=0.1, se=FALSE, method="loess")

ggplot( aes(x=t.num, y=inc),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Incorrect responses") +
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

# AD: we should not need manual censoring

#remove bins with no data
# censored <- c(54:55,124:126,156:162,230:234,249:252,335:342,359:360,371:378,411:414,449:450,540:541,610:612,642:648,716:720,735:738,821:828,845:846,857:864,897:900,935:936)
# censoredtemp <- censored*1000/binsize;
# c.all[c(censoredtemp),c('response')] <- "NA"
# c.all[c(54,124:126,156:162,230:234,249:252,335:342,359:360,371:378,411:414,449:450,540,610:612,642:648,716:720,735:738,821:828,845:846,857:864,897:900,935,936),c('response')] <- "NA"
c.all$response <- as.numeric(c.all$response)

pdf("Response curves gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = c.all) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses (log)") +
  stat_smooth(data=c.all, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter")
dev.off()
ggplot(aes(x=t.num, y=response, color = type),data = c.all) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses") +
  geom_count(position = "jitter")

# merge with genotype
# corr/inc separated
hdf <- merge(c,g) # horizontal data frame
#View(hdf)

ggplot(aes(x=t.num, y=corr, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess", size = 1) + theme_bw(base_size=20) +
  xlab("Time, 200s bins") + ylab("Correct responses") +
  geom_count(position = "jitter")

ggplot( aes(x=t.num, y=inc, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess",  size = 1) + theme_bw(base_size=20) +
  xlab("Time, 200s bins") + ylab("Incorrect responses") +
  geom_count(position = "jitter")


# all together
bdf <- merge(c.all,g)
#View(bdf)

names(cfos)[names(cfos)=="grooming time"] <- "grooming.time"

# merge with cell counts
hdfc <- merge(hdf,cfos)
bdfc <- merge(bdf,cfos)


pdf("Log response curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = bdf) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter") +
  facet_grid(Genotype ~ .)
dev.off()

pdf("Linear response curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=response, color = type),data = bdf) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses") +
  geom_count(position = "jitter") +
  facet_grid(Genotype ~ .)
dev.off()


#  quick look at genotype
# correct
summary(mcorrg1 <- glm(corr ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.corr), data = hdf))
#incorrect
summary(mincg1 <- glm(inc ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.inc), data = hdf))

summary(lmc1 <- lm(cfos$tot_correct ~ cfos$genotype01))
summary(lmi1 <- lm(cfos$tot_incorrect ~ cfos$genotype01))
summary(lmg1 <- lm(log(cfos$grooming.time) ~ cfos$genotype01))

# both
# calculate theta for both responses
# AD: this was previously not working, and all the models were run on the old theta estimate, which was incorrect
# LM 082717: yes I see that now the Genotype x type x time is highly significant
theta.resp <- theta.ml(na.omit(c.all$response), mean(na.omit(c.all$response)), 972, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)
summary(mrespg1 <- glm(response ~ 1 + t.num*type + type*Genotype +  (1:id), family = negative.binomial(theta = theta.resp), data = bdf))
# LM 082717 CAR::ANOVA output did not provide analysis related to genotype x type x time ineteraction?
car::Anova(mrespg1)

# plot over time by genotype
#lsmip(mrespg1, Genotype ~ 1 + t.num | type, at = list(t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
ls.respg1 <- lsmeans(mrespg1,"type", by = "Genotype")
plot(ls.respg1, type ~ response, horiz=F,ylab = "Response levels", xlab = "type")



# grooming time -- I understand this to be the index of OCD-like behavior
summary(mrespg2 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + grooming.time*type + grooming.time*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg2)


#DLS (alone) model
#LM 082717: did not show DLS* time *type * genotype test, and reason why the model was set up to not test this?
summary(mrespg3 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + DLS*t.num*Genotype + DLS*type*Genotype + DLS*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg3)
summary(mrespg3full <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + DLS*t.num*Genotype*type +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
anova(mrespg3full,mrespg3, test = "LR")

lsmip(mrespg3, DLS ~ type | Genotype, at = list(DLS = c(10,100,200)), ylab = "log(response rate)", xlab = "type ", type = "predicted" )
lsmip(mrespg3, DLS ~ Genotype | type, at = list(DLS = c(10,100,200)), ylab = "log(response rate)", xlab = "type ", type = "predicted" )
lsmip(mrespg3, DLS ~ t.num | type, at = list(DLS = c(10,100,200),t.num  = c(1,length(bins)/2, length(bins))), ylab = "log(response rate)", xlab = "time ", type = "predicted" )

# AD: DLS activity associated with learning of correct and has minimal impact on incorrect
lsmip(mrespg3, DLS*type ~ t.num |  Genotype , at = list(DLS = c(10,100,200),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "type", type = "predicted" )


# stronger extinction at stronger activity, more perseverative at lower activity

# AD: this is a NS interaction, do not plot
# ls.respg3 <- lsmeans(mrespg3,"t.num", by = "DLS", at = list(DLS = c(10,100,200),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])))
# plot(ls.respg3, type ~ response, horiz=F,ylab = "Response levels", xlab = "time, 200s bins")

# better reversal at higher DLS cfos levels
ls.respg3a <- lsmeans(mrespg3,"type", by = "DLS", at = list(DLS = c(10,100,200)))
plot(ls.respg3a, type ~ response, horiz=F,ylab = "Response levels", xlab = "Response type")


# AD: this a NS interaction, do not plot
# stronger extinction in WT as a Fx(DLS) than in KO
# lsmip(mrespg3, DLS ~ t.num | Genotype, at = list(DLS = c(10,100,200),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )


#NAccS
summary(mrespg4 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + NAccS*t.num*type*Genotype + (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg4)
# AD: please see this 4-way plot.  Even though the 4-way interaction is NS, it helps you better understand the data
lsmip(mrespg4, NAccS*type ~ t.num | Genotype , at = list(NAccS = c(50,150,250),t.num = c(1,length(bins)/2, length(bins))), ylab = "log(response rate)", xlab = "Time ", type = "predicted" )
lsmip(mrespg4, NAccS ~ type | Genotype , at = list(NAccS = c(50,150,250)), ylab = "log(response rate)", xlab = "response type ", type = "predicted" )
lsmip(mrespg4, NAccS ~ Genotype | type , at = list(NAccS = c(50,150,250)), ylab = "log(response rate)", xlab = "response type ", type = "predicted" )

summary(mrespg4_restr <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + NAccS*t.num*type + NAccS*t.num*Genotype + NAccS*type*Genotype + (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg4_restr)
lsmip(mrespg4_restr, NAccS ~ t.num | Genotype, at = list(NAccS = c(50,150,250),t.num  = c(1,length(bins)/2, length(bins))), ylab = "log(response rate)", xlab = "Genotype", type = "predicted" )
lsmip(mrespg4_restr, NAccS ~ Genotype | type , at = list(NAccS = c(50,150,250)), ylab = "log(response rate)", xlab = "response type ", type = "predicted" )


summary(mrespg4c <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + NAccC*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg4c)
# Genotype * NAccC X T.NUM intercation
summary(mrespg4D <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + NAccC*t.num*Genotype + NAccC*type*Genotype + NAccC*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg4D)
lsmip(mrespg4D, NAccC ~ t.num | Genotype, at = list(NAccC = c(50,150,300),t.num  = c(1,length(bins)/2, length(bins))), ylab = "log(response rate)", xlab = "Genotype", type = "predicted" )

lsmip(mrespg4c, NAccC ~ Genotype , at = list(NAccC = c(50,150,300)), ylab = "log(response rate)", xlab = "Genotype ", type = "predicted")
lsmip(mrespg4c, NAccC ~ t.num | type , at = list(NAccC = c(50,150,300),t.num  = c(1,length(bins)/2, length(bins))), ylab = "log(response rate)", xlab = "Time ", type = "predicted" )
lsmip(mrespg4c, NAccC ~ type | Genotype , at = list(NAccC = c(50,150,300)), ylab = "log(response rate)", xlab = "response type ", type = "predicted" )
lsmip(mrespg4c, NAccC ~ Genotype | type , at = list(NAccC = c(50,150,300)), ylab = "log(response rate)", xlab = "response type ", type = "predicted" )

hist(cfos$NAccC)

# summary(mrespg4a <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + NAccC*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
# car::Anova(mrespg4a)


# anova(mrespg3,mrespg4, test = "Rao")

# looking at PFC ROIS

#lOFC: no interactions found
#LM 082717: using "DLS-style" (see 5a) model showed lOFC*genoytpe*type interaction
summary(mrespg5 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + lOFC*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg5, lOFC ~ type | Genotype , at = list(lOFC = c(400,800,1200)), ylab = "log(response rate)", xlab = "Type ", type = "predicted" )
lsmip(mrespg5, lOFC ~ Genotype | type , at = list(lOFC = c(400,800,1200)), ylab = "log(response rate)", xlab = "Type ", type = "predicted" )


car::Anova(mrespg5, type = "III")
summary(mrespg5a <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + lOFC*t.num*Genotype + lOFC*type*Genotype + lOFC*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg5a)
hist(cfos$lOFC)
lsmip(mrespg5a, lOFC ~ Genotype | type , at = list(lOFC = c(400,800,1200)), ylab = "log(response rate)", xlab = "Type ", type = "predicted" )


#PrL influences response rate in KO mice, not WT mice. Low PrL actiivty associated with more perseverative behaviour, high PrL associated with better extinction
#IL-NAcS typically mediates extinction
## AD: you cannot really interpret this model because it aggregates over all mice
# LM 082717: yep was just mucking around with lm, can't remember exactly why....
lm_PrL <- lm(response ~ PrL*t.num*type*Genotype, data=bdfc); summary(lm_PrL)
lm <- lm(response ~t.num*type*Genotype, data=bdfc); summary(lm)
anova(lm_PrL, lm)

#LM 082717 again, using DLS-style model (6a) find  type x gene x PrL and type x time X PrL interactions
summary(mrespg6 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + PrL*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg6)

summary(mrespg6a <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + PrL*t.num*Genotype + PrL*type*Genotype + PrL*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg6a)

# AD: careful, this is qualified by higher-order interactions including type
#lsmip(mrespg6, PrL ~ t.num | Genotype , at = list(PrL = c(200,400,600),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
lsmip(mrespg6a, PrL ~ type | Genotype , at = list(PrL = c(200,400,600)), ylab = "log(response rate)", xlab = "Type ", type = "predicted" )

## save regression stats into a table
# quick and dirty example, default is coefficient/SE
stargazer(mrespg6a,type="html", out="PrL.htm")

# pretty example: note the very low p cutoffs, you can do (0.05,0.01, 0.001) or so; Please double-check variable labels
stargazer(mrespg6a,  type="html", digits = 2 ,single.row=TRUE,  star.cutoffs = c(0.05, 0.005, 0.0001), report = 'vtcs*',
          dep.var.labels=c("Responses"), covariate.labels=c("Time","Type: correct/incorrect", "Genotype: WT vs. SAPAP3 KO",
                                                            "ROI",   "Time*type", "Time*genotype", "Type*genotype",  "Time*ROI",
                                                            "Genotype*ROI", "Type*ROI", "Time*genotype*ROI", "Type*genotype*ROI",
                                                            "Time*type*ROI"), out="prl_pretty.htm")
# it can also combine multiple models, but unfortunately not in the way you want
stargazer(mrespg5a,mrespg6a, type="html", out="lOFC_PrL.htm")


#IL: no interaction found
summary(mrespg7 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + IL*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))

#mOFC: Trend for genotype * mOFC * response time or type interactions: pattern (response time) looks like NAcS
# LM 082717: again using similar model structure to DLS (8a) I found gene x type x mOFC interaction
summary(mrespg8 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + mOFC*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg8)

summary(mrespg8a <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + mOFC*t.num*Genotype + mOFC*type*Genotype + mOFC*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg8a)
# AD: careful, this is a NS interaction
# lsmip(mrespg8, mOFC ~ t.num | Genotype , at = list(mOFC = c(400,700,1000),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
lsmip(mrespg8, mOFC ~ type | Genotype , at = list(mOFC = c(400,700,1000)), ylab = "log(response rate)", xlab = "Type", type = "predicted" )

lsmip(mrespg8, mOFC ~ Genotype | type , at = list(mOFC = c(400,700,1000)), ylab = "log(response rate)", xlab = "Genotype", type = "predicted" )



#M2: no interaction found
# AD: comes through if you cut the time at 12000 s
summary(mrespg9 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + M2*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))

#averaged region analyses
#Dorsal striatum (DLS/DMS/CMS): similar to DLS
summary(mrespg11 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val1*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg11)
lsmip(mrespg11, val1 ~ t.num | Genotype , at = list(val1 = c(-2,0,2),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time", type = "predicted" )
lsmip(mrespg11, val1 ~ type | Genotype , at = list(val1 = c(-2,0,2)), ylab = "log(response rate)", xlab = "type", type = "predicted" )
lsmip(mrespg11, val1 ~ t.num | type , at = list(val1 = c(-2,0,2),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time", type = "predicted" )

#better reversal with DS Cfos
ls.respg11a <- lsmeans(mrespg11,"type", by = "val1", at = list(val1 = c(-2,0,2)),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)]))
plot(ls.respg11a, type ~ response, horiz=F,ylab = "Response levels", xlab = "Response type")
plot(ls.respg11a, type ~ response | t.num , horiz=F,ylab = "Response levels", xlab = "Response type")
#DCS
summary(mrespg11C <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val7*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg11C)

ls.respg11b <- lsmeans(mrespg11,"Genotype", by = "val1", at = list(val1 = c(-2,0,2)))
plot(ls.respg11b, type ~ response , horiz=F,ylab = "Response levels", xlab = "Genotype")

#Ventral striatum (NAcC/NAcS/VMS): Not similar to NAcS
summary(mrespg12 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val2*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))

#Nucleus accumbens (NAcC/NAcS): Similar to NAcS
summary(mrespg13 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val3*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg13, val3 ~ t.num | Genotype , at = list(val3 = c(-1.5,0,1.5),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time", type = "predicted" )
car::Anova(mrespg12)

#mPFC (mOFC/PrL): Similar to NAcS (not similar to PrL if IL or IL/lOFC are added to PrL)
summary(mrespg16 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + valm6*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg16, valm6 ~ t.num | Genotype , at = list(valm6 = c(-2,0,2),t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time", type = "predicted" )
lsmip(mrespg16, valm6 ~ type | Genotype , at = list(valm6 = c(-2,0,2)), ylab = "log(response rate)", xlab = "Type", type = "predicted" )

#oFC (mOFC/lOFC): trend genotype x cfos interaction (no interaction with behaviour)
summary(mrespg15 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val5*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
#PFC(OFC and mPFC)
summary(mrespg17 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val7*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg17)

#PrL-NAcS circuit
summary(mrespg17 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + valm9*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg17)

#stuff for deciding how to graph interactions
hist(just_rois$M1)
hist(just_rois$NAccS)
hist(cfos$CMS)
hist(bdfc$grooming.time)
sd(cfos$NAccS)
median(cfos$NAccS)
mean(cfos$NAccS)
mean(cfos$val3)

#checking regions quickly
summary(mrespg10 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + DLS*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
summary(mrespg10 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + NAccC*t.num*Genotype + NAccC*type*Genotype + NAccC*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg10)
lsmip(mrespg10, CMS ~ type | Genotype , at = list(grooming.time = c(25,75,150)), ylab = "log(response rate)", xlab = "Type", type = "predicted" )
lsmip(mrespg10, CMS ~ type | Genotype , at = list(CMS = c(25,75,150)), ylab = "log(response rate)", xlab = "Type", type = "predicted" )


# conclusion: strong effects of cfos, weak effect of genotype, + interactions

# pretty graph
library(multcompView)

leastsquare = lsmeans(mrespg3, pairwise ~ DLS:Genotype:type, at = list(DLS = c(10,100,200)), adjust='tukey')
library(dplyr)

CLD = cld(leastsquare, alpha=0.05, Letters=letters,
          adjust='tukey')
CLD$type <- recode_factor(CLD$type,corr="Correct",inc="Incorrect")
pdf(file = "DLS plot pretty PRETTY.pdf", width = 10, height = 6)
pd = position_dodge(15)    ### How much to jitter the points on the plot
ggplot(CLD, aes( x = DLS, y = lsmean, color = Genotype, label = .group)) +
  facet_wrap( ~ type) +
  
  geom_point(shape  = 15,
             size   = 4,
             position = pd) +
  scale_y_reverse() +
  
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
  ) +  ylab("Log-probability of response") + xlab("Dorsolateral striatum cfos") +
  ggtitle ("Behavior in reversal phase by genotype, regional cfos level, and response type") +
  labs( caption  = paste0(
      "\n",
      "Boxes indicate the LS mean.\n",
      "Error bars indicate the 95% ",
      "confidence interval of the LS mean, Sidak method for 8 estimates. \n",
      "Means sharing a letter are ",
      "not significantly different ",
      "(Tukey-adjusted comparisons for 8 estimates)."),
    hjust = 0.5 ) +
  geom_text(nudge_x = c(20, 10, 20, 10, 10, 10,10, 10, 10, 10, 10, 10),
              nudge_y = c(0,  0, 0,  0, 0 , 0,0,0, .05,  -.05, 0, 0), color   = "black") +
  # nudge_y = 0, color   = "black") +
  
    scale_color_manual(values = c("blue", "red"))
dev.off()
