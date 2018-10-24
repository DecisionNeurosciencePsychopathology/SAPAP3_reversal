## LME on binned response data from Lizzie Manning and Suzanne Ahmari's free operant deterministic reversal experiment

library(readr)
library(lme4)
library(ggplot2)
#library(dplyr)
library(tidyverse)
library(psych)
library(gdata)
library(R.matlab)
library(xtable)
library(Hmisc)
library(nnet)
library(reshape2)
# library(ggbiplot)
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
library(ggraph)
library(corrr)
library(igraph)
library(ggpubr)
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

# WT
just_rois <- cfos[cfos$genotype01==0,7:18]
#head(just_rois)
ROIcorWT <- just_rois %>%
  correlate(method = "spearman") %>%
  stretch()
graph_cors <- ROIcorWT %>%
   filter(abs(r) > .618) %>%
  graph_from_data_frame(directed = F)
g1 <- ggraph(graph_cors, layout = "igraph", algorithm = "circle") +
  geom_edge_link(aes(edge_alpha = r, edge_width = r, color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(0.618, 1), colors = c("red", "yellow")) +
  geom_node_point(color = "light grey", size = 14) +
  geom_node_text(aes(label = name), repel = F, size = 3) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = "Connectivity in WT")
# KO
just_rois <- cfos[cfos$genotype01==1,7:18]
#head(just_rois)
ROIcorKO <- just_rois %>%
  correlate(method = "spearman") %>%
  stretch()
graph_cors <- ROIcorKO %>%
  filter(abs(r) > .618) %>%
  graph_from_data_frame(directed = F)
g2 <- ggraph(graph_cors, layout = "igraph", algorithm = "circle") +
  geom_edge_link(aes(edge_alpha = r, edge_width = r, color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(0.618, 1), colors = c("red", "yellow")) +
  geom_node_point(color = "light grey", size = 14) +
  geom_node_text(aes(label = name), repel = F, size = 3) +
  theme_graph(base_family = 'Helvetica') +
  labs(title = "Connectivity in KO")
gall <- ggarrange(g1,g2,ncol = 2, common.legend = T)
ggsave("differential_connectivity.pdf", width = 12, height = 6.5)

h <- t.test(ROIcorKO$r,ROIcorWT$r, paired = TRUE)

ggsave("connectivity in WT.pdf", width = 8, height = 8)