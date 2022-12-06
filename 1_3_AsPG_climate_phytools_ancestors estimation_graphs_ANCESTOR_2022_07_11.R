####################################################################
###                     CLIMATE ANCESTORS ESTIMATION
####################################################################

########################################
#               LIBRARIES
########################################
library(ggplot2)
library(ggstance)
library(tidyverse)
library(wesanderson)
library(scales)

# check this!! https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

########################################
#               PATHS
########################################

# ###PORTATIL VV
# path.to.data <- "file:///C:/Users//VV.5039878//Desktop//In prep_1_Coca_AsPG_Climate//Data//"
# path.to.results <- "file:///C:/Users//VV.5039878//Desktop//In prep_1_Coca_AsPG_Climate//Results//"

####MAC DESPACHO VV
# path.to.data <- "/Users/vvalnun/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Data/"
# path.to.results <- "/Users/vvalnun/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Results/"
# path.to.environment <- "/Users/vvalnun/Dropbox/HEDERA/MAC-DESPACHO/HEDERA/ESTADISTICA/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Env/"

### PORTATIL WINDOWS NGM
main.path <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/"
setwd(main.path)

path.to.data <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Data/"
path.to.results <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Results/"

# save date for version control
date <- "2022_07_11"

########################################
#                 DATA
########################################

# climatic data of the tips
aral <- read.csv("Data/Araliaceae_clima_PCA_2021_01_10.csv")

########################################################'
## Run twice one with nucleo and the other with plasto #####
########################################################'
# First run
run <- "nucleo"
anc_clim <- read.csv(paste0(path.to.results, "ancestors_climate_nucleo", date, ".csv"))
anc_clim_prob <- read.csv(paste0(path.to.results, "ancestors_climate_prob_nucleo", date, ".csv"))

# Second run
run <- "plasto"
anc_clim <- read.csv(paste0(path.to.results, "ancestors_climate_plasto", date, ".csv"))
anc_clim_prob <- read.csv(paste0(path.to.results, "ancestors_climate_prob_plasto", date, ".csv"))

# extract limits for each transitional region
quantiles_lat <- aral %>%
  group_by(C1_Latitude) %>%
  summarise(.,minPC1 = min(PC1), Q10 = quantile(PC1, prob = 0.10), Q90 = quantile(PC1, prob = 0.90), maxPC1 = max(PC1), na.rm=TRUE)

limits_trans_lat <- c(t_sub_temp_low = quantiles_lat$Q10[1], # temp
                      t_sub_temp_up = quantiles_lat$Q90[2], # subt-temp
                      t_sub_trop_low = quantiles_lat$Q10[3], # subt-temp
                      t_sub_trop_up = quantiles_lat$Q90[1] # subt-temp
                      )

# summarize data
summary_aral <- aral %>%
  group_by(Genus)%>%
  summarize(meanPC1=mean(PC1, na.rm = TRUE),meanPC2=mean(PC2, na.rm = TRUE),
            medianPC1=median(PC1, na.rm = TRUE), medianPC2= median(PC2, na.rm = TRUE),
            maxPC1=max(PC1), minPC1 = min(PC1),
            maxPC2=max(PC2), minPC2= min(PC2)) 

# change name of node in anc_clim to match name in anc_clim_prob
anc_clim$node <- paste0("X", anc_clim$node)
anc_clim_prob <- tidyr::gather(anc_clim_prob, node, value, factor_key=TRUE)

# colors for the nodes by mean PC value
# PC mean value by node
anc_clim_prob <- anc_clim_prob %>%
  group_by(node) %>%
  mutate(mean_PC=mean(value))

# extract limits for each transitional region
qlat <- aral %>%
  group_by(C1_Latitude) %>%
  summarise(.,minPC1 = min(PC1), Q10 = quantile(PC1, prob = 0.10), Q90 = quantile(PC1, prob = 0.90), maxPC1 = max(PC1), na.rm=TRUE)

limits <- c(t_sub_temp_low = qlat$Q10[1], # temp
            t_sub_temp_up = qlat$Q90[2], # subt-temp
            t_sub_trop_low = qlat$Q10[3], # subt-temp
            t_sub_trop_up = qlat$Q90[1] # subt-temp
            )

# select colors 
# color palette
# to establish the number of colors in each ramp
dlimits <- abs(round(100*(c(-6, limits)-c(limits,4))))

# to generate colors use this link https://colordesigner.io/gradient-generator
# create one ramp per category
coltemp <- colorRampPalette(c("#5d478b", "#5864ae", "#4b82cd", "#33a1e8",
                              "#00bffe", "#00bffe", "#00caff", "#00d4ff",
                              "#00deff", "#00e7ef","#00eed6", "#00f5b5"))(dlimits[1])
coltemp_subt <- colorRampPalette(c("#00fa8e", "#00fd5e", 
                                   "#00ff05", "#00ff05", "#5bfb00", "#7ff700",
                                   "#9af300", "#b0ef00"))(dlimits[2])
colsubt <- colorRampPalette(c("#c4ea00", "#d5e600", "#e4e100", "#f2dc00",
                              "#fed700", "#fed700", "#ffcf02", "#ffc706",
                              "#ffc00c", "#ffb812"))(dlimits[3])
colsubt_trop <- colorRampPalette(c("#ffb118", "#ffa91d", "#ffa222", "#fe9a26", 
                                   "#fc932a", "#fc932a", "#fc8a21", "#fd8018",
                                   "#fd760f", "#fd6b06"))(dlimits[4])
coltrop <- colorRampPalette(c("#fe5f00", "#fe5100", "#fe4100", "#fe2d00",
                              "#fe0000", "#fe0000", "#ea0000", "#d60000",
                              "#c30000", "#b00000"))(dlimits[5])

# paste together de colors of all the ramps
color <- c(coltemp, coltemp_subt,colsubt,colsubt_trop,coltrop)

# standardize meanPC1 values to rage between 1 and 100, extremes = limits of the PCA
color01 <- 1000* round(rescale(c(-6,unique(anc_clim_prob$mean_PC),4)), digits = 2)

# select color proportional to value using the standardized values
color_anc <- color[color01]

ggplot() +
  # horizontal boxplots & density plots
  geom_density(data = anc_clim_prob, aes(x= value, fill = node)) +
  scale_fill_manual(values = color_anc) +
  # results of the other reconstruction methods
  # geom_linerangeh(data = anc_clim,   aes(y = -0.3, xmin = minPC1,
  #                 xmax = maxPC1), colour = "lightgrey", size = 2) +
  geom_linerangeh(data = anc_clim,   aes(y = -0.3, xmin = medianPC1-0.05,
                                       xmax = medianPC1+0.05), colour = "red", size = 2) +
  geom_linerangeh(data = anc_clim,   aes(y = -0.3, xmin = meanPC1-0.05,
                                          xmax = meanPC1+0.05), colour = "orange", size = 2) +
  # horizontal lines with limits between regions
  geom_vline(xintercept = limits_trans_lat, linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 1)) +
  xlim(c(-6, 4)) + # same range as in TIPS FIGURE
  ylab("") +
  xlab("First axis of the PCA") +
  facet_grid(node ~ .) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA),
        strip.text.y = element_text(angle=360)
        )
# ggsave(paste0(path.to.results, "Climate_Anc_M_", run, "_probabilistic_with_ranges_newramp", date,".svg"), width = 8, height = 10) 


# # each genera in one plot
# # i <- (1:length(unique(anc_clim_prob$node)))[[1]]
# n <- length(unique(anc_clim_prob$node))
# 
# lapply(1:n, function(i){
#   col <- color_anc[[i]]
#   node <- unique(anc_clim_prob$node)[[i]]
#   anc_clim_sel <- anc_clim[anc_clim$node == node,]
#   anc_clim_prob_sel <- anc_clim_prob[anc_clim_prob$node == node,]
#   
#   ggplot() +
#     # horizontal boxplots & density plots
#     geom_density(data = anc_clim_prob_sel, aes(x= value),fill = col) +
#     # geom_linerangeh(data = anc_clim,   aes(y = -0.3, xmin = minPC1,
#     #                xmax = maxPC1), colour = "lightgrey", size = 2) +
#     geom_linerangeh(data = anc_clim_sel,   aes(y = 0, xmin = medianPC1-0.1,
#                                                xmax = medianPC1+0.1), colour = "red", size = 6) +
#     geom_linerangeh(data = anc_clim_sel,   aes(y = 0, xmin = meanPC1-0.1,
#                                                xmax = meanPC1+0.1), colour = "orange", size = 6) +
#     # geom_vline(xintercept = limits_trans_lat) +
#     scale_y_continuous(breaks = c(0,0.5, 1)) +
#     ylim(c(0,1.6)) +
#     xlim(c(-6,4)) +
#     ylab("") +
#     xlab("") +
#     theme_bw() +
#     theme(
#       panel.background = element_rect(fill='transparent'), #transparent panel bg
#       plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
#       panel.grid.major = element_blank(), #remove major gridlines
#       panel.grid.minor = element_blank(), #remove minor gridlines
#       legend.background = element_rect(fill='transparent'), #transparent legend bg
#       legend.box.background = element_rect(fill='transparent'), #transparent legend panel
#       panel.border = element_blank(), # remove panel border
#       axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
#       axis.line.y = element_line(size = 1, linetype = "solid", colour = "white"),
#       axis.ticks.x = element_blank(),
#       axis.ticks.y = element_blank(),
#       axis.text.y = element_blank()
#     )
#   ggsave(paste0("Results/pdf_nodes_",run , node, ".png"), width = 12, height = 4)    
# })
# 
# 
# # probar esto
# # https://www.datanovia.com/en/blog/elegant-visualization-of-density-distribution-in-r-using-ridgeline/
# 
