##################################################################################
# This script was used to calculate the niche centroids and produce the 
# main results of the peper "Tree diversity increases functional stability of European forests"
# Corresponding author: Marco Patacca, pataccamarco@gmail.com

#------------- Libraries -------------------------------------------------------
library(tidyverse)
library(here)
library(stringr)
library(terra)
library(foreign)
library(BAT)
library(colorspace)
library(cowplot)
library(feather)
library(svglite)
library(sf)
library(ggmap)
#for PCA
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.19")
library(aroma.light)
library(Redmonder)
library(ggrepel)
library(wCorr)
library(pals)
library(RColorBrewer)
library(labeling)

#------------ Functions --------------------------------------------------------
# Load taylormade functions 
source(here("Scripts", "Functions.R"))

#------------- Trait ordination of single trees --------------------------------
# Principal component (Fig.2) of imputed single trees. 
# To have more insights in the imputation method refer to the code in the supplementary material of 
# Maynard et al., 2022 - https://doi.org/10.1038/s41467-022-30888-2 

# Load imputed trait data
read_feather(here("Output","Trait imputed dataset", "Optimized_estimated_trait_table.4.feather"))->imp_tr
colnames(imp_tr)<- c("accepted_bin", "LAT", "LON", "trait1","trait2", "trait3", "trait4","trait5", "trait6","trait7","trait8", "trait14","trait15","trait37",
                     "trait44","trait46","trait48","trait66","trait185","trait281","trait324", "trait719",
                     "trait773","trait1312","trait3110","trait3117","trait24","trait45")

# Read in trait names
use_traits <- read.csv(here("Output", "Intermediate outputs", "Traits", "Trait_lookup.csv")) %>% arrange(TraitID) %>% rename(trait = TraitName) %>%  mutate(trait_label = paste0("trait", TraitID))
use_traits %>% dplyr::select(trait_label, trait_short) %>% filter(trait_label != "trait26")->use_traits

imp_tr %>%
  gather(key=trait, value=value, -accepted_bin, -LAT, -LON) %>% 
  left_join(use_traits, by=c("trait"="trait_label"))->imp_tr

#load angio gymno
ag <- read.csv(here("Data", "ANGIO_GYMNO_lookup.csv")) %>%
  mutate(angio = ifelse(group == "Angiosperms", 1, ifelse(group == "Gymnosperms",0, NA))) %>%
  filter(!is.na(genus)) %>% dplyr::select(accepted_bin, group, angio)

# Tidy and scale for pca
imp_tr %>% 
  dplyr::select(-trait, -LAT, -LON) %>% 
  group_by(trait_short) %>% 
  mutate(value = my_scale(value)) %>%
  ungroup()->imp_tr

# Tidy imputed trait and weight for species number
imp_tr %>% rename(trait=trait_short) %>% 
  spread(trait, value) %>% 
  left_join((.) %>% group_by(accepted_bin) %>% tally() %>% ungroup %>% mutate(spp_n = n) %>% dplyr::select(-n), by = "accepted_bin") %>%
  group_by(accepted_bin) %>% mutate(spp_wt = (1/spp_n)/sum(1/spp_n)) %>% ungroup %>%
  mutate(wt = spp_wt) %>%
  mutate(wt = wt/sum(wt)) %>%
  left_join(ag0 %>% select(accepted_bin, order, family)) %>%
  dplyr::select(accepted_bin, LAT, LON, all_of(use_traits$trait_short), wt) -> imp_tr_weight

# Apply fun PCA1-2
ggpca_weighted_1.2(tr_use = imp_tr_weight, 
                   trait_vars = use_traits$trait_short, 
                   angio_list = ag, line_wd_big = 1,show_labs = TRUE, 
                   num_show = 6, use_vars = use_traits$trait_short[c(4,5,16,9,12,15)], alpha = 0, label_size = 4, flip_coord = FALSE,
                   flip_x = -1, flip_z = -1, flip_zz = -1, flip_y = 1,
                   show_legend = TRUE, title = NULL)->pca.dan12

# Apply fun PCA2-3   
ggpca_weighted_2.3(tr_use = imp_tr_weight, 
                   trait_vars = use_traits$trait_short,
                   angio_list = ag, line_wd_big = 1,show_labs = TRUE, 
                   num_show = 6, use_vars = use_traits$trait_short[c(4,5,16,8,14,2)], alpha = 0, label_size = 4, flip_coord = FALSE,
                   flip_x = -1, flip_z = -1, flip_zz = -1, flip_y = 1,
                   show_legend = TRUE, title = NULL)->pca.dan23

# Plot the two pca together (Fig.2)
plot_grid(pca.dan12$plot,pca.dan23$plot, nrow = 1) -> pca.plots

#save
pdf(here("Output", "Figures", "final", "Figure_2.pdf"), width=15, height=10)
pca.plots
dev.off()

#------------ Compute functional hyperspace ------------------------------------
# The hyperspace is based on trait species averages
#load df combining NFI data with community structure and functional traits data
readRDS(here("Data", "NFIs", "Remeasured_minmax_NFI.RDS"))->rem_NFI

# Column names of rem_NFI (accepted_bin is the species name, uniqueID is an ID unique for each single observation):
# [1] "country"      "LAT"          "LON"          "plot.id"      "census.n"     "plot.obs.id"  "ba"           "accepted_bin" "trait1"       "trait2"       "trait3"       "trait4"      
# [13] "trait5"       "trait6"       "trait7"       "trait8"       "trait14"      "trait15"      "trait37"      "trait44"      "trait46"      "trait48"      "trait66"      "trait185"    
# [25] "trait281"     "trait324"     "trait719"     "trait773"     "trait1312"    "trait3110"    "trait3117"    "trait24"      "trait45"      "uniqueID"     "plot.ba"      "rel.ba" 

# Set the two different methods to build the hyperspace
distance_gower <-c("gower")
distance_euc <-c("euclidean")

# Calculate hyperspaces for all spp
hyperspace(rem_NFI, distance_gower)-> hyperspace_gow
hyperspace(rem_NFI, distance_euc)-> hyperspace_euc

# Calculate the quality of the hyperspace 
## - euclidean distance
distance<- dist(pop.trait)
hyper.quality(distance, hyperspace_euc)->hyperspace_quality_euc     # 0.96
## - gower distance
distance_gow<- gower(pop.trait)
hyper.quality(distance_gow, hyperspace_gow)->hyperspace_quality_gow # 0.97

#------------- Compute community hypervolumes (occupied niche) -----------------
# Apply fun to calculate Hypervolumes for all the communities of remeasured (min-max census) plots
Niche_calc(rem_NFI, hyperspace_gow)->FDis_3d

# Load plot data
readRDS(here("Data", "NFIs", "Remeasured_minmax_NFI.RDS")) %>% 
  dplyr::select(plot.id, country, LAT, LON, census.n, uniqueID) %>% unique()->plot_info

# merge with FDis df
FDis_3d  %>% 
  left_join(plot_info)->FDis_3d

#add standardized census tags 
FDis_3d %>% 
  dplyr::select(-uniqueID) %>% 
  group_by(country) %>% 
  mutate(census.type = case_when(census.n == min(census.n) ~ "first_NFI",
                                 census.n == max(census.n) ~ "last_NFI")) %>% 
  ungroup() -> FDis_3d

# Calculate relative ba of spp per plot
rem_NFI %>%
  dplyr::select(country, LAT, LON, plot.id, census.n, ba, accepted_bin) %>% 
  left_join(ag) %>% 
  group_by(country, census.n, plot.id) %>% 
  mutate(plot.ba = sum(ba)) %>% 
  ungroup() %>% 
  group_by(country, census.n, plot.id, angio) %>% 
  mutate(AG.rel.plot.ba = sum(ba)*100/plot.ba) %>% 
  ungroup()->ag.plot

# Calculate spp diversity per plot x census
ag.plot %>% 
  group_by(country) %>% 
  mutate(census.type = case_when(census.n == min(census.n) ~ "first_NFI",
                                 census.n == max(census.n) ~ "last_NFI")) %>% 
  ungroup() %>% 
  dplyr::select(country, census.type, plot.id, accepted_bin) %>% unique() %>% 
  group_by(country, census.type, plot.id) %>% 
  count()->spp.x.plot

# Merge back with main df 
FDis_3d %>% left_join(spp.x.plot) %>% rename(n_spp_plot = n) ->FDis_3d

#------------ Make maps (fig.3) ------------------------------------------------
# Map centroids
Map_FD_metric(FDis_3d, "C.axis.1", legend.name = "LES")->c.1
Map_FD_metric(FDis_3d, "C.axis.2", legend.name = "SAS")->c.2
Map_FD_metric(FDis_3d, "C.axis.3", legend.name = "RUS")->c.3

# Change axis name
c.1+ labs(tag = "A", fill="", x="", y="Latitude")->c.1
c.2+ labs(tag = "B", fill="", x="Longitude", y="")->c.2
c.3+ labs(tag = "C", fill="", x="", y="")->c.3

# Combine maps together (Fig.3)
plot_grid(c.1,c.2,c.3, nrow =1, ncol = 3,  rel_widths = c(1,1,1), rel_heights = c(1,1,1))->niches_map

# Save
pdf(here("Output", "Figures", "final", "Figure_3.pdf"),width=15, height=10)
niches_map
dev.off()

#----- Plot niche centroid loadings on axes (2nd census) vs n spp --------------
# Axis 1, LES
FDis_3d %>% 
  filter(census.type == "last_NFI") %>% 
  filter(n_spp_plot<=8) %>% 
  ggplot(aes(x=n_spp_plot, y=C.axis.1, col =C.axis.1))+
  geom_jitter(alpha=0.2, show.legend = FALSE)+
  geom_hline(yintercept = 0, col="black")+
  theme(legend.position="none", panel.grid.major = element_blank())+
  scale_colour_continuous_divergingx(palette = "Spectral", na.value = NA)+
  scale_x_continuous(breaks = seq(1, 8, by = 1))+
  labs(x =" ", y= "LES", tag = "A")+
  geom_smooth(method=lm)->axis1

# Axis 2, SAS
FDis_3d %>% 
  filter(census.type == "last_NFI") %>% 
  filter(n_spp_plot<=8) %>% 
  ggplot(aes(x=n_spp_plot, y=C.axis.2, col =C.axis.2))+
  geom_jitter(alpha=0.2, show.legend = FALSE)+
  geom_hline(yintercept = 0, col="black")+
  theme(legend.position="none", panel.grid.major = element_blank())+
  scale_colour_continuous_divergingx(palette = "Spectral", na.value = NA )+
  scale_x_continuous(breaks = seq(1, 8, by = 1))+
  labs(x ="Number of species in the latest census", y= "SAS", tag = "B")+
  geom_smooth(method=lm)->axis2

# Axis 3, RUS
FDis_3d %>% 
  filter(census.type == "last_NFI") %>% 
  filter(n_spp_plot<=8) %>% 
  ggplot(aes(x=n_spp_plot, y=C.axis.3, col =C.axis.3))+
  geom_jitter(alpha=0.2, show.legend = FALSE)+
  geom_hline(yintercept = 0, col="black")+
  theme(legend.position="none", panel.grid.major = element_blank())+
  scale_colour_continuous_divergingx(palette = "Spectral", na.value = NA )+
  scale_x_continuous(breaks = seq(1, 8, by = 1))+
  labs(x =" ", y= "RUS", tag = "C")+
  geom_smooth(method=lm)->axis3

# Combine plots together
plot_grid(axis1, axis2, axis3, nrow = 1, ncol = 3)->Figure_4

# Save
pdf(here("Output", "Figures", "final", "Figure_4.pdf"),width=15, height=10)
Figure_4
dev.off()

#----------- Plot absolute rate of change vs spp n 1st census (Fig.5) ----------
# Add NFI average interval to compute the rate
NFI_interval<-c(16.74, 9.63, 5.93, 5.36, 5, 10.96, 4.74, 10, 7.50, 5, 5, 11.31, 5, 4)
country<-c("Flanders", "Wallonia", "Czech Republic", "Denmark", "France", "Germany", "Ireland", "Italy", "Netherland", "Norway", "Poland", "Spain", "Sweden", "Findland")
cbind(country,NFI_interval) %>%  as.data.frame() -> interval.nfi
interval.nfi$NFI_interval %>% as.numeric()->interval.nfi$NFI_interval

# Calculate niche centroid change over time
FDis_3d %>%  
  dplyr::select(contains("Axis"), census.type, plot.id, country, LAT, LON) %>% 
  pivot_wider(names_from = census.type, values_from = contains("Axis")) %>% left_join(interval.nfi) %>% na.omit() %>% 
  mutate(axis.1.change = (C.axis.1_last_NFI - C.axis.1_first_NFI)/NFI_interval,
         axis.2.change = (C.axis.2_last_NFI - C.axis.2_first_NFI)/NFI_interval,
         axis.3.change = (C.axis.3_last_NFI - C.axis.3_first_NFI)/NFI_interval) -> Centroids.df

# Calculate plot  spp changes over time between census
spp.x.plot %>% 
  pivot_wider(names_from = census.type, values_from = n) %>% 
  rename(spp_NFI_1 = first_NFI, spp_NFI_2 = last_NFI) %>% 
  mutate(spp_change = spp_NFI_2 - spp_NFI_1)->spp_change_df

# Merge centroid change df with spp change df
Centroids.df %>%
  left_join(spp_change_df) -> Centroids.df

# Plot absolute rate of change on LES
Centroids.df %>% filter(spp_NFI_1<=8) %>% filter(spp_change !=0) %>% 
  ggplot(aes(x=spp_NFI_1, y = abs(axis.1.change)))+
  geom_jitter(aes(col=as.factor(spp_change)), alpha=0.7)+
  geom_smooth(method = 'glm', formula =y ~ log(x))+
  scale_colour_discrete_divergingx(palette = "Spectral", na.value = NA)+
  labs(x ="", tag="A")+
  theme(legend.position="none", panel.grid.major = element_blank())+
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.025), name = "LES change")+
  scale_x_continuous(breaks = seq(1, 12, by = 1)) -> change.1

# Plot absolute rate of change on SAS
Centroids.df %>%  filter(spp_NFI_1<=8) %>% filter(spp_change !=0) %>% 
  ggplot(aes(y = abs(axis.2.change),x=spp_NFI_1))+
  geom_jitter(aes( col=as.factor(spp_change)), alpha=0.7)+
  geom_smooth(method = 'glm', formula =y ~ log(x))+
  scale_colour_discrete_divergingx(palette = "Spectral", na.value = NA)+
  labs(x ="Number of species 1st census", tag="B")+
  theme(legend.position="none", panel.grid.major = element_blank())+
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.025), name = "SAS change")+
  scale_x_continuous(breaks = seq(1, 12, by = 1)) -> change.2

# Plot absolute rate of change on RUS
Centroids.df %>%  filter(spp_NFI_1<=8) %>% filter(spp_change !=0) %>% 
  ggplot(aes(y = abs(axis.3.change),x=spp_NFI_1))+
  geom_jitter(aes( col=as.factor(spp_change)), alpha=0.7)+
  geom_smooth(method = 'glm', formula =y ~ log(x))+ 
  scale_colour_discrete_divergingx(palette = "Spectral", na.value = NA)+
  labs(x ="",  tag="C", col = "Species changes")+
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.025), name = "RUS change")+
  scale_x_continuous(breaks = seq(1, 12, by = 1))+
  theme(panel.grid.major = element_blank())-> change.3

# Combine plots together (Fig.5)
plot_grid(change.1, change.2, change.3, nrow = 1, rel_widths = c(0.8, 0.8, 1))->Figure_5

#save
pdf(here("Output", "Figures", "final", "Figure_5.pdf"),width=15, height=10)
Figure_5
dev.off()
