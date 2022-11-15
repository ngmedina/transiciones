####################################################################@
###                     CLIMATE ANCESTORS ESTIMATION ####
####################################################################@

########################################@
#               LIBRARIES
########################################@
library(phytools)
require(data.table)
library(ggplot2)
library(ape)

########################################@
#               PATHS
########################################@

# ###PORTATIL VV
# path.to.data <- "file:///C:/Users//VV.5039878//Desktop//In prep_1_Coca_AsPG_Climate//Data//"
# path.to.results <- "file:///C:/Users//VV.5039878//Desktop//In prep_1_Coca_AsPG_Climate//Results//"
main.path <- "/Users/virginiavalcarcel/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In\ prep_Coca_Tropical\ vs.\ temperate\ Araliacae_Statistics"
setwd(main.path)
path.to.data <- "/Users/virginiavalcarcel/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In\ prep_Coca_Tropical\ vs.\ temperate\ Araliacae_Statistics/Data/"
path.to.results <- "/Users/virginiavalcarcel/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In\ prep_Coca_Tropical\ vs.\ temperate\ Araliacae_Statistics/Results/"


####MAC DESPACHO VV
main.path <- "/Users/vvalnun/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/"
setwd(main.path)
path.to.data <- "/Users/vvalnun/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Data/"
path.to.results <- "/Users/vvalnun/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Results/"
path.to.environment <- "/Users/vvalnun/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Env/"

## PORTATIL WINDOWS NGM
main.path <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/"
setwd(main.path)
path.to.data <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Data/"
path.to.results <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Results/"

########################################@
#                 DATA              #####
########################################@
# AsPG random climate values
aral <- read.csv("Data/Araliaceae_clima_PCA_2021_01_10.csv")
aral$Genus<-as.character(aral$Genus)
aral <- aral[!(aral$Genus =="Didymopanax"), ]
aral <- aral[!(aral$Genus =="Cephalopanax"), ]
aral <- aral[!(aral$Genus =="Frodinia"), ]
aral <- aral[!(aral$Genus =="Crepinella"), ]
aral$Genus<-as.factor(aral$Genus)

# load trees nucleo
nucleo_ultra <- ape::read.tree('Data/nucleo_ultra.newick')

# load trees plasto
plasto_ultra <- ape::read.tree('Data/plasto_ultra.newick')

# save date for version control
date <- "2022_07_11"

########################################################'
## Run twice one with nucleo and the other with plasto #####
########################################################'
# First run
run <- "nucleo"
ultra <- nucleo_ultra

# Second run
run <- "plasto"
ultra <- plasto_ultra

###################################################

####################################################################@
###           ANCESTORS CLIMATE ESTIMATION - PROBABILISTIC      ####
####################################################################@

# split aral in matrices by genus
spl_aral <- split(aral, f = aral$Genus)

# check the selected genus
names <- names(spl_aral)

# calculate density kernels
pdf_of_data <- lapply(spl_aral, function(mat) density(mat$PC1, from= 0, to=1, bw=0.1))

# IMPORTANT!! ######'
# RUN ONLY ONCE: Save the result and use the saved result (see code below)

# # Number of iterations to run
# N <- 10000
# 
# # This part of the code selects N points in the distribution of each taxon proportionally to their frequency
# # for debugging:
#  nm <- names(pdf_of_data)[1]
#  x.new <- lapply(names(pdf_of_data), function(nm) {
#    dens_tmp <- pdf_of_data[[nm]]
#    aral_tmp <- spl_aral[[nm]]
#    rnorm(N, sample(aral_tmp$PC1, size = N, replace = TRUE), dens_tmp$bw)
#  }
#  )
#  
# # # transform result from list to matrix
#  x.new <- do.call(cbind, x.new)
# # 
# # # names of species
# colnames(x.new) <- names
# write.csv(x.new, paste0(path.to.data, "/rnd1.csv"))

### END OF RUN ONLY ONCE #####'

# we save the selection of points in case we need it later on
x.new <- read.csv(paste0(path.to.data, "rnd1.csv"), as.is = T)
x.new <- x.new[,-1]

#select the correct nuclear tree
data.order <-  match(colnames(x.new),ultra$tip.label)
x.new <-  x.new[,data.order]

# we reconstruct the ancestral niches for the N combinations of points 
# for debugging
# rnd <- x.new[1,]
anc_clim_prob <- apply(x.new, 1, function(rnd){
  fastAnc(ultra,rnd,vars=TRUE,CI=TRUE)#estimate ancestral states with 95CI
}
)

# reformat the results to a matrix
# x <- anc_clim[[1]]
anc_clim_mt_prob <- lapply(anc_clim_prob, function(x) x$ace)
anc_clim_mt_prob <- do.call(rbind, anc_clim_mt_prob)

anc_clim_prob <- as.data.frame(anc_clim_mt_prob)

# save results
write.csv(anc_clim_mt_prob, paste0(path.to.results, "ancestors_climate_prob_", run, date, ".csv"), row.names=FALSE)

# # nodes potentially interesting to plot: N7 --> 25 (Den,(Che,Gam)); N14 --> 34 (Hed,Mer); N12 --> 30 (Kal,(Met,Mac))

###########################################################################@
####                    CLIMATE ANCESTORS ESTIMATION - SEMIQUANTITATIVE ####
###########################################################################@

########################################@
#               LIBRARIES
########################################@

## first triy just with this. If you get problem 1, then go to lines 170-173
library(phytools)
library(dplyr)
library(tidyverse)#to use column 1 as rownames
library(ggstance)
library(viridis)
library(wesanderson)


# load data to skip previous steps
anc_clim_prob <- read.csv(paste0(path.to.results, "ancestors_climate_prob_" ,run, date, ".csv"))

##########################################################@
#                 MEAN, MEDIAN AND RANGE VALUES       ####
##########################################################@
summary_aral <- aral %>%
  group_by(Genus)%>%
  summarize(meanPC1=mean(PC1, na.rm = TRUE),meanPC2=mean(PC2, na.rm = TRUE),
            medianPC1=median(PC1, na.rm = TRUE), medianPC2= median(PC2, na.rm = TRUE),
            maxPC1=max(PC1), minPC1 = min(PC1),
            maxPC2=max(PC2), minPC2= min(PC2)) 

# match tree and clim names
data.order <- match(ultra$tip.label, summary_aral$Genus)
summary_aral <- summary_aral[data.order, ]

# reconstruct ancestral values
var <- summary_aral[,2][[1]]
anc_clim_list <- apply(summary_aral[c("meanPC1","medianPC1","maxPC1","minPC1")], 2, function(var) 
  fastAnc(ultra, var, vars = TRUE, CI = TRUE)
  )

# extract clim values of ancestors
anc_clim <- lapply(anc_clim_list, '[[', 1)
anc_clim <- as.data.frame(t(do.call(rbind, anc_clim)))
anc_clim$node <- rownames(anc_clim)

# write.csv(anc_clim, paste0(path.to.results, "ancestors_climate_", run, date, ".csv"), row.names=FALSE)

# Plot ancestral reconstructions
# colores igual que NGM
color <- wes_palette("Zissou1", 100, type = "continuous")

pdf(paste0("Climate_Anc_ML_MeanPC1_", run, date, "_color.pdf"))
x <- summary_aral$meanPC1
names(x) <- summary_aral$Genus
objPC1 <- contMap(ultra,x, method = "anc.ML", plot=FALSE)
objPC1 <- setMap(objPC1,colors=color)#to change colours 
plot(objPC1,legend=0.7*max(nodeHeights(ultra)),
     fsize=c(0.7,0.9))
dev.off()


pdf(paste0("Climate_Anc_ML_MeanPC2_", run, date, "_color.pdf"))
x <- summary_aral$meanPC2
names(x) <- summary_aral$Genus

objPC2 <- contMap(ultra,x ,method = "anc.ML",plot=FALSE)
objPC2 <- setMap(objPC2,colors=color)#to change colours 
plot(objPC2,legend=0.7*max(nodeHeights(ultra)),
     fsize=c(0.7,0.9))
dev.off()


####################################################################@
###                     CLIMATE ANCESTORS ESTIMATION - DISCRETE #####
####################################################################@

########################################@
#               LIBRARIES
########################################@


## first try just with this. If you get problem 1, then go to lines below
library(phytools)

## only to solve problem 1 
# library(devtools)
# install_github("liamrevell/phytools",quiet=TRUE) # only if the tree has polytomies
# library(phytools)
library(ape)

########################################@
#                 DATA              
########################################@
# AsPG climatic characterization of genera under expert criterion (Temperate versus tropical, column 2)
data <- read.table ("data_discrete.txt", dec=",", row.names=1, header=T, sep="\t")

run <- "plasto"
ultra <- nucleo_ultra
run <- "nucleo"
ultra <- plasto_ultra

# name changes
# Asian Schefflera: Heptapleurum
# Neotropical Schefflera: Sciadophyllum
rownames(data)[rownames(data) =="Schefflera_Asian"] <- "Heptapleurum"
rownames(data)[rownames(data) == "Schefflera_Neotropical"] <- "Sciadophyllum"

## Problem 1: The tree has polytomies and ace does not work and retrieve Error in ace(data.ordered[, 2], phy, type = "discrete", model = "ER") : 
#"phy" is not rooted AND fully dichotomous. And when using phy<-multi2di(phy) a new error appears when runnning ace: Error in ace(data.ordered[, 2], phy, type = "discrete", model = "ER") : 
#some branches have length zero or negative
## Solution: Liam Revell has this solution (http://blog.phytools.org/2015/09/imperceptible-update-to-rerootingmethod.html):
phy<-multi2di(ultra) 
phy$edge.length[ultra$edge.length==0]<-1e-8 # only if the tree has polytomies

##Problem 2: taxa order in table does not match tip order in the tree
## Solution: 
data.order = match(ultra$tip.label, rownames(data))
data.ordered = data[data.order, ]
#vec<-c("Tetrapanax","Heteropanax","Schefflera_Asian","Schefflera_Neotropical","Gamblea","Chengiopanax","Dendropanax","Brassaiopsis","Eleutherococcus","Macropanax","Metapanax","Kalopanax","Fatsia","Oreopanax","Merrilliopanax","Hedera","Oplopanax")
#rename_labels <- function(tree, vec) {
  #phy$tip.label <- vec
  #return(tree)
#}
#phy2 <- rename_labels(tree = phy, vec = vec)

###############################################################################@
#           ANCESTORS CLIMATE ESTIMATION USING ML CHARACTER RECONSTRUCTION #####
###############################################################################@

#ER model (you can repeat, SYM, ARD, etc...). 

ERreconstructionER <- ace(data.ordered[,1],phy, type="discrete", model="ER")## column 1 is as rownames
ERreconstructionER$loglik
ERreconstructionER$rates
ERreconstructionER$lik.anc #da las probabilidades de cada estado del car?cter en cada nodo

pdf(paste0("Climate_Anc_ML_ER_discrete_", run, date, ".pdf"))
plot(phy)
tiplabels(pch=22, bg= c("red","red","blue","red","red","red","red","blue","blue","red","red","red","blue","blue","red","red","blue","blue"), frame="n", cex=1.5)
nodelabels(pie=ERreconstructionER$lik.anc,piecol=setNames(c("blue","red"),
                                                         c("a","b")),cex=0.6)
dev.off()

#SYM model (you can repeat, SYM, ARD, etc...). 

ERreconstructionSYM <- ace(data.ordered[,1],phy, type="discrete", model="SYM")## column 1 is as rownames
ERreconstructionSYM$loglik
ERreconstructionSYM$rates
ERreconstructionSYM$lik.anc #da las probabilidades de cada estado del car?cter en cada nodo

pdf(paste0("Climate_Anc_ML_SYM_discrete_", run, date, ".pdf"))
plot(phy)
tiplabels(pch=22, bg= c("red","red","blue","red","red","red","red","blue","blue","red","red","red","blue","blue","red","red","blue","blue"), frame="n", cex=1.5)
nodelabels(pie=ERreconstructionSYM$lik.anc,piecol=setNames(c("blue","red"),
                                                          c("a","b")),cex=0.6)
dev.off()


#ARD model (you can repeat, SYM, ARD, etc...). 

ERreconstructionARD <- ace(data.ordered[,1],phy, type="discrete", model="ARD")## column 1 is as rownames
ERreconstructionARD$loglik
ERreconstructionARD$rates
ERreconstructionARD$lik.anc #da las probabilidades de cada estado del car?cter en cada nodo

pdf(paste0("Climate_Anc_ML_ARD_discrete_", run, date, ".pdf"))
plot(phy)
tiplabels(pch=22, bg= c("red","red","blue","red","red","red","red","blue","blue","red","red","red","blue","blue","red","red","blue","blue"), frame="n", cex=1.5)
nodelabels(pie=ERreconstructionARD$lik.anc,piecol=setNames(c("blue","red"),
                                                           c("a","b")),cex=0.6)
dev.off()


#######Arboles planos para base figura cambiando tiplabel

pdf(paste0("tree", run, ".pdf"))
plot(ultra)
dev.off()

pdf(paste0("tree_", run, "_con node label.pdf"))
plot(ultra)
nodelabels()
dev.off()

