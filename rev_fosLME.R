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


# read in the data
setwd("/Users/manninge/Documents/GitHub/SAPAP3_reversal/")
d <- read_excel("/Users/manninge/Documents/GitHub/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
#setwd("~/code/SAPAP3_reversal/")
# d <- read_excel("~/code/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
#View(d)


# read in genotype
g <- read_excel("/Users/manninge/Documents/GitHub/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
#g <- read_excel("~/code/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
g$id <- g$`Physical Tag`
g <- g[,2:3]


# read in cfos cell counts
cfos <- read_csv("/Users/manninge/Documents/GitHub/SAPAP3_reversal/SAPAP3 reversal cFos density_Final mice.csv")
#cfos <- read_csv("~/code/SAPAP3_reversal/SAPAP3 reversal cFos density_Final mice.csv")
names(cfos)[names(cfos)=="ID"] <- "id"
names(cfos)[names(cfos)=="correct"] <- "tot_correct"
names(cfos)[names(cfos)=="incorrrect"] <- "tot_incorrect"
names(cfos)[names(cfos)=="Genotype"] <- "genotype01"
names(cfos)[names(cfos)=="NAc S"] <- "NAccS"
names(cfos)[names(cfos)=="NAc C"] <- "NAccC"


#View(cfos)
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

#Re-merge cell counts with other data

# run PCA (code from Alex 081717)
#ds.pca = prcomp((just_rois),scale = TRUE, center = TRUE)
# and write component scores
#ds_pcas <- get_pca_ind(ds.pca)
#cfos$val1 <- cfos_pcas$coord[,1]


# run PCA and write component scores
cfos_pcas <- get_pca_ind(cfos.pca)
# hdf$val1 <- cfos_pcas$coord[,1]
# hdf$val2 <- cfos_pcas$coord[,2]
# hdf$val3 <- cfos_pcas$coord[,3]
# hdf$val4 <- cfos_pcas$coord[,4]

# look
summary(cfos.pca)

# check the variance explained by various factors
plot(cfos.pca,type = 'l')

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
end <- 18000
start <- 0
binsize <- 1000
bins <- seq(start,end,binsize)

c <- data.frame()
for (id in ids)
{t=cut(d$corr[d$id==id],c(bins), labels = FALSE)
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

ggplot(aes(x=t.num, y=corr),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Correct responses") +
  stat_smooth(data=c, aes(x=t.num, y=corr), size=0.4, alpha=0.1, se=FALSE, method="loess")

ggplot( aes(x=t.num, y=inc),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Incorrect responses") +
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
pdf("Response curves gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = c.all) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses (log)") +
  stat_smooth(data=c.all, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter")
dev.off()
ggplot(aes(x=t.num, y=response, color = type),data = c.all) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses") +
  geom_count(position = "jitter")

# merge with genotype
# corr/inc separated
hdf <- merge(c,g) # horizontal data frame
#View(hdf)

ggplot(aes(x=t.num, y=corr, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess", size = 1) + theme_bw(base_size=20) +
  xlab("Time, 100s bins") + ylab("Correct responses") +
  geom_count(position = "jitter")

ggplot( aes(x=t.num, y=inc, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess",  size = 1) + theme_bw(base_size=20) +
  xlab("Time, 100s bins") + ylab("Incorrect responses") +
  geom_count(position = "jitter")


# all together
bdf <- merge(c.all,g)
#View(bdf)

names(cfos)[names(cfos)=="grooming time"] <- "grooming.time"

# merge with cell counts
hdfc <- merge(hdf,cfos)
bdfc <- merge(bdf,cfos)


pdf("Log esponse curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = bdf) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter") +
  facet_grid(. ~ Genotype)
dev.off()

pdf("Linear esponse curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=response, color = type),data = bdf) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses") +
  geom_count(position = "jitter") +
  facet_grid(. ~ Genotype)
dev.off()


#  quick look at genotype
# correct
summary(mcorrg1 <- glm(corr ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.corr), data = hdf))
#incorrect
summary(mincg1 <- glm(inc ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.inc), data = hdf))

summary(lmc1 <- lm(cfos$tot_correct ~ cfos$genotype01))
summary(lmi1 <- lm(cfos$tot_incorrect ~ cfos$genotype01))
summary(lmg1 <- lm(cfos$grooming.time ~ cfos$genotype01))

# both
# calculate theta for both responses
theta.resp <- theta.ml(c.all$response, mean(c.all$response), 1026, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)
summary(mrespg1 <- glm(response ~ 1 + t.num*type+ type*t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.resp), data = bdf))
car::Anova(mrespg1)

# plot over time by genotype
#lsmip(mrespg1, Genotype ~ 1 + t.num | type, at = list(t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
ls.respg1 <- lsmeans(mrespg1,"type", by = "Genotype")
plot(ls.respg1, type ~ response, horiz=F,ylab = "Response levels", xlab = "type")



# grooming time -- I understand this to be the index of OCD-like behavior
summary(mrespg2 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + grooming.time*t.num*type +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg2)


#DLS (alone) model
summary(mrespg3 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + DLS*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg3)

lsmip(mrespg3, DLS ~ t.num | type, at = list(DLS = c(10,100,200),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )

# stronger extinction at stronger activity, more perseverative at lower activity
ls.respg3 <- lsmeans(mrespg3,"t.num", by = "DLS", at = list(DLS = c(10,100,200),t.num = c(1000, 9000, 18000)))
plot(ls.respg3, type ~ response, horiz=F,ylab = "Response levels", xlab = "time, 100s bins")

# better reversal at higher DLS cfos levels
ls.respg3a <- lsmeans(mrespg3,"type", by = "DLS", at = list(DLS = c(10,100,200)))
plot(ls.respg3a, type ~ response, horiz=F,ylab = "Response levels", xlab = "Response type")


# stronger extinction in WT as a Fx(DLS) than in KO
lsmip(mrespg3, DLS ~ t.num | Genotype, at = list(DLS = c(10,100,200),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )



summary(mrespg4 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + NAccS*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg4)
lsmip(mrespg4, NAccS ~ t.num | Genotype , at = list(NAccS = c(50,150,300),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
lsmip(mrespg4, NAccS ~ type | Genotype , at = list(NAccS = c(50,150,300)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )

# anova(mrespg3,mrespg4, test = "Rao")

# looking at PFC ROIS

#lOFC: no interactions found
summary(mrespg5 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + lOFC*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg5)

#PrL influences response rate in KO mice, not WT mice. Low PrL actiivty associated with more perseverative behaviour, high PrL associated with better extinction
#IL-NAcS typically mediates extinction
summary(mrespg6 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + PrL*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg6, PrL ~ t.num | Genotype , at = list(PrL = c(200,400,600),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )

#IL: no interaction found
summary(mrespg7 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + IL*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))

#mOFC: Trend for genotype * mOFC * response time or type interactions: pattern (response time) looks like NAcS
summary(mrespg8 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + mOFC*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg8, mOFC ~ t.num | Genotype , at = list(mOFC = c(400,700,1000),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
lsmip(mrespg8, mOFC ~ type | Genotype , at = list(mOFC = c(400,700,1000), ylab = "log(response rate)", xlab = "Type", type = "predicted" )
lsmip(mrespg8, mOFC ~ type | Genotype , at = list(mOFC = c(400,700,1000)), ylab = "log(response rate)", xlab = "Type", type = "predicted" )
lsmip(mrespg8, mOFC ~ Genotype , at = list(mOFC = c(400,700,1000)), ylab = "log(response rate)", xlab = "Genotype", type = "predicted" )
      


#M2: no interaction found
summary(mrespg9 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + M2*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))

#averaged region analyses
#Dorsal striatum (DLS/DMS/CMS): similar to DLS
summary(mrespg11 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val1*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg11, val1 ~ t.num | Genotype , at = list(val1 = c(-2,0,2),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
lsmip(mrespg11, val1 ~ t.num | type , at = list(val1 = c(-2,0,2),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )

#better reversal with DS Cfos
ls.respg11a <- lsmeans(mrespg11,"type", by = "val1", at = list(val1 = c(-2,0,2)),t.num = c(1000, 9000, 18000))
plot(ls.respg11a, type ~ response, horiz=F,ylab = "Response levels", xlab = "Response type")
plot(ls.respg11a, type ~ response | t.num , horiz=F,ylab = "Response levels", xlab = "Response type")

ls.respg11b <- lsmeans(mrespg11,"Genotype", by = "val1", at = list(val1 = c(-2,0,2)))
plot(ls.respg11b, type ~ response , horiz=F,ylab = "Response levels", xlab = "Genotype")

#Ventral striatum (NAcC/NAcS/VMS): Not similar to NAcS
summary(mrespg12 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val2*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))

#Nucleus accumbens (NAcC/NAcS): Similar to NAcS
summary(mrespg13 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val3*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg13, val3 ~ t.num | Genotype , at = list(val3 = c(-1.5,0,1.5),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )


#mPFC (mOFC/PrL): Similar to NAcS (not similar to PrL if IL or IL/lOFC are added to PrL)
summary(mrespg16 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + valm6*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg16, valm6 ~ t.num | Genotype , at = list(valm6 = c(-2,0,2),t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
lsmip(mrespg16, valm6 ~ type | Genotype , at = list(valm6 = c(-2,0,2)), ylab = "log(response rate)", xlab = "Type", type = "predicted" )

#oFC (mOFC/lOFC): trend genotype x cfos interaction (no interaction with behaviour)
summary(mrespg15 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + val5*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))

#stuff for deciding how to graph interactions
hist(just_rois$M1)
hist(just_rois$NAccS)
hist(cfos$mOFC)
hist(bdfc$t.num)
sd(cfos$val1)
median(cfos$val1)
mean(cfos$val1)
mean(cfos$val3)

#checking regions quickly
summary(mrespg10 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + M1*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
lsmip(mrespg10, M1 ~ type | Genotype , at = list(M1 = c(10,60,110)), ylab = "log(response rate)", xlab = "Type", type = "predicted" )




# conclusion: strong effects of cfos, weak effect of genotype, + interactions
