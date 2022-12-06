####################################################################
###                     CLIMATE TIPS ESTIMATION
####################################################################

########################################'
#               LIBRARIES          #####
########################################'

# library(wesanderson)
library(dplyr)
library(ggplot2)
library(ggstance)
library(ggtree)
library(cowplot)
library(scales)
library(tidyr)

########################################'
#               PATHS             #####
########################################'

### PORTATIL WINDOWS NGM

main.path <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/"
setwd(main.path)

path.to.data <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Data/"
path.to.results <- "C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/2020_1_In prep_Coca_Tropical vs. temperate Araliacae_Statistics/Results/"

########################################'
#               LOAD DATA #####
########################################'

aral <- read.csv("Data/Araliaceae_clima_PCA_2021_01_10.csv")

# AsPG tip probabilistic climate values
aral <- read.csv("Data/Araliaceae_clima_PCA_2021_01_10.csv")
aral$Genus<-as.character(aral$Genus)
aral <- aral[!(aral$Genus =="Didymopanax"), ]
aral <- aral[!(aral$Genus =="Cephalopanax"), ]
aral <- aral[!(aral$Genus =="Frodinia"), ]
aral <- aral[!(aral$Genus =="Crepinella"), ]
aral$Genus<-as.factor(aral$Genus)

# save date for version control
date <- "2022_07_11"


########################################################'
## Run twice one with nucleo and the other with plasto #####
########################################################'

#############"
# run nucleo 
#############"

run <- "nucleo"

# load tree
nucleo_ultra <- ape::read.tree('Data/nucleo_ultra.newick')
tree_ultra <- nucleo_ultra

# match Genus in tree and tip climate data
aral$Genus <- factor(aral$Genus, levels = tree_ultra$tip.label)
aral <- aral[!is.na(aral$Genus),]

# load internal nodes probabilistic climate data
anc_clim <- read.csv(paste0(path.to.results, "ancestors_climate_nucleo", date, ".csv"))
anc_clim_prob <- read.csv(paste0(path.to.results, "ancestors_climate_prob_nucleo", date, ".csv"))

# change name of node in anc_clim to match name in anc_clim_prob
colnames(anc_clim_prob) <- as.character(anc_clim$node)
# from short to long format
anc_clim_prob <- tidyr::gather(anc_clim_prob, node, value, factor_key=TRUE)

#############"
# run plasto 
#############"

run <- "plasto"

# load tree
plasto_ultra <- ape::read.tree('Data/plasto_ultra.newick')
tree_ultra <- plasto_ultra

# match Genus in tree and tip climate data
aral$Genus <- factor(aral$Genus, levels = tree_ultra$tip.label)
aral <- aral[!is.na(aral$Genus),]

# load internal nodes probabilistic climate data
anc_clim <- read.csv(paste0(path.to.results, "ancestors_climate_nucleo", date, ".csv"))
anc_clim_prob <- read.csv(paste0(path.to.results, "ancestors_climate_prob_nucleo", date, ".csv"))

# change name of node in anc_clim to match name in anc_clim_prob
colnames(anc_clim_prob) <- as.character(anc_clim$node)
# from short to long format
anc_clim_prob <- tidyr::gather(anc_clim_prob, node, value, factor_key=TRUE)



########################################'
#               CREATE PLOTS         #####
########################################'


##### Phylotenetic tree

# generate labels for tips
d <- data.frame(label = tree_ultra$tip.label)

tree <- ggtree(tree_ultra) %<+% d +  xlim(NA,1.2) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab(aes(label = label))


tree

# reorder the labels of the database to match tip labels
aral$Genus <- factor(aral$Genus, levels = get_taxa_name(tree))

# summarize climatic data
summary_aral <- aral %>%
  group_by(Genus)%>%
  summarize(meanPC1=mean(PC1, na.rm = TRUE),meanPC2=mean(PC2, na.rm = TRUE),
            medianPC1=median(PC1, na.rm = TRUE), medianPC2= median(PC2, na.rm = TRUE),
            maxPC1=max(PC1), minPC1 = min(PC1),
            maxPC2=max(PC2), minPC2= min(PC2)) 

#### Generate color gradient 

# extract limits for each transitional region
qlat <- aral %>%
  group_by(C1_Latitude) %>%
  summarise(.,minPC1 = min(PC1), Q10 = quantile(PC1, prob = 0.10), Q90 = quantile(PC1, prob = 0.90), maxPC1 = max(PC1), na.rm=TRUE)

limits <- c(t_sub_temp_low = qlat$Q10[1], # temp
                      t_sub_temp_up = qlat$Q90[2], # subt-temp
                      t_sub_trop_low = qlat$Q10[3], # subt-temp
                      t_sub_trop_up = qlat$Q90[1] # subt-temp
            )


# color palette
# to establish the number of colors in each ramp
dlimits <- abs(round(100*(c(-6, limits)-c(limits,4))))

# to generate colors we used this webpage https://colordesigner.io/gradient-generator
# we selected the maximum and minimum values so that they matched the colors in the
# maps of the paper and then some intermediate colors to build a ramp per
# category

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

df <- data.frame(y= rep(0, length = length(color)), x = seq(from =-6, to= 4, length= length(color)))

# reference bar
barplot <- ggplot(df,aes(x= x, y = y)) +
  geom_point(color = color, size=8, shape=0) +
  geom_vline(xintercept = limits, 
             linetype = "dashed") +
  xlim(c(-6, 4)) +
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none")

barplot


#### Density plot for TIPS

# select colors
# standardize meanPC1 values to rage between 1 and 100, extremes = limits of the PCA
color01 <- 1000* round(rescale(c(-6,summary_aral$meanPC1,4)), digits = 2)

# select color proportional to value using the standardized values
color_tip <- color[color01[c(-1,-length(color01))]]

# plot 
density_plot_tips <- ggplot() +
  # horizontal boxplots & density plots
  geom_density(data = aral, aes(x= PC1, fill = Genus)) +
  scale_fill_manual(values =  color_tip) +
  # uncomment to show range of values as a bar in the lower part of the density plot
  # geom_linerangeh(data = anc_clim,   aes(y = -0.3, xmin = minPC1,
  #                xmax = maxPC1), colour = "lightgrey", size = 2) +
  geom_linerangeh(data = summary_aral,   aes(y = -0.3, xmin = meanPC1-0.05,
  xmax = meanPC1+0.05), colour = "red", size = 2) +
  geom_linerangeh(data = summary_aral,   aes(y = -0.3, xmin = medianPC1-0.05,
  xmax = medianPC1+0.05), colour = "orange", size = 2) +
  geom_vline(xintercept = limits, 
             linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 1)) +
  xlim(c(-6, 4)) + 
  ylab("") +
  xlab("First axis of the PCA") +
  facet_grid(Genus ~ ., switch = "both") +
  theme_bw() +
  # uncomment to show species names in density plot
  # theme(legend.position = "none",
  #       strip.background = element_rect(fill = NA),
  #       strip.text.y.left = element_text(angle=360)) +
  # theme(strip.text = element_text(face = "italic")) +
  theme(legend.position = "none",
    strip.background = element_blank(),
    strip.text.y = element_blank()
  )

density_plot_tips

# ggsave(paste0(path.to.results, "Climate_TIPS_newramp.svg"), width = 8, height = 10) 


# PC mean value by node
anc_clim_prob <- anc_clim_prob %>%
  group_by(node) %>%
  mutate(mean_PC=mean(value))

# standardize meanPC1 values to rage between 1 and 100, extremes = limits of the PCA
color01 <- 1000* round(rescale(c(-6,unique(anc_clim_prob$mean_PC),4)), digits = 2)

# select color proportional to value using the standardized values
color_anc <- color[color01]

density_plot_anc <- ggplot() +
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
  geom_vline(xintercept = limits, linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 1)) +
  xlim(c(-6, 4)) + # same range as in TIPS FIGURE
  ylab("") +
  xlab("First axis of the PCA") +
  facet_grid(node ~ .) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA),
        strip.text.y = element_text(angle=360))

density_plot_anc

# combine the plots 
plot_grid(barplot, density_plot_anc, ncol = 1, align = "v", axis = "l",
          rel_heights = c(0.3,1))

plot1 <- plot_grid(density_plot_anc, tree, density_plot_tips, 
                   align = "h", 
                   axis = "b", 
                   nrow = 1,
                   rel_widths = c(1,1.5,1))

plot1

# ggsave(paste0(path.to.results, "climate_TIP_ANC.svg"), width = 12, height = 10) 


# generate a plot with the axis of the two subplots aligned
barplot + density_plot_tips + 
  plot_layout(ncol = 1, heights = c(1,40))

# ggsave(paste0(path.to.results, "Climate_TIPS_withbar.svg"), width = 8, height = 12) 

# graphs one by one
n <- length(unique(aral$Genus))
lapply(1:n, function(i){
  col <- color_tip[[i]]
  genus <- unique(aral$Genus)[[i]]
  aral_sel <- aral[aral$Genus == genus,]
  summary_aral_sel <- summary_aral[summary_aral$Genus == genus,]
  
  ggplot() +
    # horizontal boxplots & density plots
    geom_density(data = aral_sel, aes(x= PC1),fill = col) +
    # geom_linerangeh(data = anc_clim,   aes(y = -0.3, xmin = minPC1,
    #                xmax = maxPC1), colour = "lightgrey", size = 2) +
    geom_linerangeh(data = summary_aral_sel,   aes(y = 0, xmin = medianPC1-0.1,
                                               xmax = medianPC1+0.1), colour = "red", size = 6) +
    geom_linerangeh(data = summary_aral_sel,   aes(y = 0, xmin = meanPC1-0.1,
                                               xmax = meanPC1+0.1), colour = "orange", size = 6) +
    scale_y_continuous(breaks = c(0,0.5, 1)) +
    ylim(c(0,1.6)) +
    xlim(c(-6,4)) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      panel.border = element_blank(), # remove panel border
      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "white"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  ggsave(paste0("Results/pdf_tips_", genus, ".png"), width = 10, height = 5)    
})

