# some very simple simulations to make sure that decreasing bin size to 20s does not result in mis-specifying the distributions
# did not have the time to make this work really smoothly, but you can see that the best bin size can be as small as 3-5s
# first, run this script, then call the function from command line, e.g. fbins(20), and it will plot histograms of simulated vs. actual data.
fbins <- function(binsize)
{
library(readr)
library(lme4)
library(ggplot2)
#library(dplyr)
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
# setwd("/Users/manninge/Documents/GitHub/SAPAP3_reversal/")
# d <- read_excel("/Users/manninge/Documents/GitHub/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
("~/code/SAPAP3_reversal/")
d <- read_excel("~/code/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
# View(d)


d$id <- d$`subject ID`
d$corr <- d$`Correct lever press`
d$inc <- d$`Incorrect lever press`
d <- d[,4:6]
d[d==0] <- NA
ids <- unique(d$id)
describe(d$corr)
describe(d$inc)
# end <- 18000
start <- 0
# binsize <- 1000
# try smaller bins:
# binsize <- 20
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

c.all <- melt(c,id.vars = c("id","time","t.num","t.ord"),value.name = "response", variable.name = "type")
names(c.all)[names(c.all)=="value"] <- "response"
names(c.all)[names(c.all)=="variable"] <- "type"

nbfit <- suppressWarnings(fitdistr(c.all$response, "negative binomial"))
print(nbfit$estimate)
# logLik(nbfit)
# AIC(nbfit, k = 2)
simulated <- rnegbin(nbfit$n,nbfit$estimate[1],nbfit$estimate[2])
actual <- c.all$response

# just simple visual diagnostics
histogram(~ simulated + actual)
}

