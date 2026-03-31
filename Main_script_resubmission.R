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
library(hypervolume)
library(colorspace)
library(cowplot)
library(feather)
library(svglite)
library(sf)
library(lme4)
library(lmerTest)
library(ncdf4)
library(data.table)
library(broom.mixed)
library(lmerTest)
library(purrr)

# # To install the tools to run the bayseian models
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#cmdstanr::check_cmdstan_toolchain(fix = TRUE)
#cmdstanr::install_cmdstan()
library(brms)
options(brms.backend = "cmdstanr")

# #for weighted PCA
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("aroma.light")
library(aroma.light)

library(ggrepel)
library(RColorBrewer)
library(funspace)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

#------------ Functions --------------------------------------------------------
# Load taylor-made functions 
source(here("Scripts", "4-Final code", "Functions.R"))

#------------- Trait ordination of single trees (Europe scale)------------------
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
  left_join(ag %>% select(accepted_bin, order, family)) %>%
  dplyr::select(accepted_bin, LAT, LON, all_of(use_traits$trait_short), wt) -> imp_tr_weight

# Apply fun PCA1-2
ggpca_weighted_1.2(tr_use = imp_tr_weight, 
                   trait_vars = use_traits$trait_short, 
                   angio_list = ag, line_wd_big = 1,show_labs = TRUE, 
                   num_show = 6, use_vars = use_traits$trait_short[c(4,5,16,9,12,15)], alpha = 0, label_size = 4, flip_coord = FALSE,
                   flip_x = -1, flip_z = -1, flip_zz = -1, flip_y = 1,
                   show_legend = TRUE, title = NULL)->pca.12


#save
pdf(here("Output", "Figures", "final", "Figure_2.pdf"), width=15, height=10)
pca.plots
dev.off()

#------------- Functional space biogeographical regions ------------------------
############################################    BIOREGIONS
#load shapefile with bioregions from https://www.eea.europa.eu/en/analysis/maps-and-charts/biogeographical-regions-in-europe-2
read_sf(here("Data", "Abiotics", "EEA_biogeoregions", "BiogeoRegions2016.shp")) -> biogeo  
st_transform(biogeo, crs = 4326)->biogeo # reproject crs

# Remove biogeo reg not included
biogeo %>% 
  filter(code!= "Anatolian") %>% 
  filter(code!= "Arctic") %>% 
  filter(code!= "BlackSea") %>% 
  filter(code!= "Macronesia") %>% 
  filter(code!= "Steppic") %>% 
  filter(code!= "Outside") %>% 
  filter(code!= "Pannonian") -> biogeo

# make plot coords a spatial obj.
st_as_sf(imp_tr_weight, coords = c("LON", "LAT"), crs=crs(biogeo))->plot_biogeo
# extract plots 
plot_biogeo %>% mutate(uniqueID = 1:length(accepted_bin)) %>%  st_join(biogeo)->plot_biogeo 
# select vars and tidy up
plot_biogeo %>% dplyr::select(uniqueID, accepted_bin, code) %>% st_drop_geometry() %>% unique()->spp_biogeo
  
# plot bio regions
# Example: Define a bounding box manually (xmin, ymin, xmax, ymax)
bbox <- st_bbox(c(xmin = -11, ymin = 35, xmax = 35, ymax = 71), crs = st_crs(4326))
st_crop(biogeo, bbox)->crop_biogeo
#plot
crop_biogeo %>% 
  ggplot(aes(fill = code))+
  geom_sf(color = "black") +
  scale_fill_manual(values = rev(hcl.colors(5, "Roma")))+
  theme_void()+
  labs(fill= "Biogeographical regions")

# Define colors for each region (reversed for better contrast)
region_colors <- rev(hcl.colors(length(unique(crop_biogeo$code)), "Roma"))

# Create a named vector for color mapping
color_mapping <- setNames(region_colors, unique(crop_biogeo$code))

# Plot
for (region in unique(crop_biogeo$code)) {
  print(plot_region(crop_biogeo, region)) # Plot biogeoregions
}

######################### convert wpca output for funspace
# PC 1, 2
# get the mean value to scale the eigenvectors to mean unit variance
sdmean <- sqrt(mean(pca.12[["pca"]][["d"]]^2))
# extract the objects, scale the eigenvalues and x matrix
sdev <- pca.12[["pca"]][["d"]]/sdmean  # d = eigenvalues
rotmat <- t(pca.12[["pca"]][["vt"]])   # vt = eigenvectors
xmat <- pca.12[["pca"]][["pc"]]/sdmean # pc = scores
# Rotate the axes to match pca figure
xmat[,1] <- -xmat[,1]
rotmat[,1] <- -rotmat[,1]
#make the first two axis a df for funspace plot
data.frame(PC1 = xmat[,1], PC2 = xmat[,2]) ->xmat 
# Scale values for a direct comparison with vectors
xmat %>% 
  mutate(PC1 = scale11(PC1),
         PC2 = scale11(PC2)) %>% 
  cbind(spp_biogeo) %>% 
  filter(complete.cases(.), abs(PC1) < 1, abs(PC2) < 1)-> funspace_pc12  # Remove missing & extreme values
######################## This piece of code is embedded in the ggpca_weighted functions 

# Calc functional space bioregions o wpca axis 1 and 2
funspace(funspace_pc12, PCs = c(1,2), group.vec = funspace_pc12$code)->comm_funtest

# Summary of functional metrics
summary(comm_funtest)

############################ compute functional space bioregions on axis 1 and 2 
# Extract and prepare axes to plot:
# Select only relevant top 3 trait variables for PC1, 
 pca.12$loadings %>% as.data.frame() %>% 
  dplyr::select(trait, PC1, PC2) %>%
  mutate(abs_PC1 = abs(PC1)) %>%  # Compute absolute values
  arrange(desc(abs_PC1)) %>%  # Sort by absolute PC1
  slice_head(n = 3) %>% 
   dplyr::select(-abs_PC1)-> main.trait.pc1
# and PC2
pca.12$loadings %>%
  dplyr::select(trait, PC1, PC2, PC3) %>%
  mutate(abs_PC2 = abs(PC2)) %>%
  arrange(desc(abs_PC2)) %>%  # Sort by absolute PC2
  slice_head(n = 3)%>% 
  dplyr::select(-abs_PC2)->main.trait.pc2  # Take top 3 for PC2
# merge together
main.trait.pc1 %>%  
 bind_rows(main.trait.pc2) -> top_loadings
# Extract full loadings (not absolute values) for PC1 & PC2 (since plotting is 2D)
all_loadings <- pca.12$loadings[, 2:3]  
arrow_scaling <- 0.65  # Adjust this value to increase/decrease arrow size
all_loadings_plot <- all_loadings * arrow_scaling  # Scale only for plotting (loadings remain unchanged)
top_loadings_plot <- top_loadings[, 2:3] * arrow_scaling  # Scale only for plotting (loadings remain unchanged)
# Save pca % var ex for plot
xlab<-(paste0("PC1 (", round(pca.12[["axis_perc"]][1]*100,1), "%)"))
ylab<-(paste0("PC2 (", round(pca.12[["axis_perc"]][2]*100,1), "%)"))

# plot for bioregions
plot(x = comm_funtest, #funspace object
     type = "groups", 
     quant.plot = TRUE, #add quantile lines
     arrows = F, #add arrows for PCA loadings
     quant.col = c("gray39" ),
     quant.lwd = 2.5,
     quant.labels = T,
     colors = c('gray85', "#1F28A2"  ),  #Colors are not looped
     globalContour = TRUE,
     globalContour.lwd = 2,
     globalContour.lty = 1,
     globalContour.col = "black",
     plot.panels = TRUE,
     axis.title = F,
     xlim = c(-1,1),
     ylim = c(-1,1))
title(xlab = xlab, cex =4)
title(ylab = ylab, cex =4)
abline(h = 0, v = 0, col = "black", lwd = 1.5, lty = 1) 

#------------ Compute functional hyperspace ------------------------------------
# The hyperspace is based on trait species averages
#load df combining NFI data with community structure and functional traits data
readRDS(here("Data", "NFIs", "Remeasured_minmax_NFI.RDS"))%>% group_by(country) %>% filter(census.n == max(census.n)) %>% ungroup() %>% 
  dplyr::select(uniqueID, accepted_bin, contains("trait"))->species.traits

# Set the two different methods to build the hyperspace
distance_gower <-c("gower")
distance_euc <-c("euclidean")

# Calculate hyperspaces for all spp
hyperspace(species.traits, distance_gower)-> hyperspace_gow
hyperspace(species.traits, distance_euc)-> hyperspace_euc

# # plot species averages with names #Supp fig.
# species.traits$accepted_bin %>% unique() %>% sort() %>%
# cbind(hyperspace_gow) %>% as_data_frame() %>% rename("accepted_bin" = ".") %>%
#   left_join(ag) -> df.plot.try
# df.plot.try %>%
#   ggplot(aes(x=as.numeric(Axis1), y=as.numeric(Axis2), col=(group)))+
#   geom_point()+ labs(col="Group", x = "PCAxis 1 (36.2%)" , y= "PCAxis 2 (14.3%)")+
#   geom_text(aes(label = accepted_bin), check_overlap = TRUE) +
#   geom_hline(yintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
#   geom_vline(xintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
#   theme_bw()+
#   theme(text = element_text(size = 20))-> species_avg_ord
# 
# # Save
# pdf(here("Output", "Figures", "Reworked", "mean_species_ordination.pdf"),width=15, height=10)
# species_avg_ord
# dev.off()

# Calculate the quality of the hyperspace 

#Summarise  per spp (average)
species.traits %>% 
  group_by(accepted_bin) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  as.data.frame() %>% 
  ungroup()->pop.trait

#make it a matrix with spp names on rows
matrix.numeric.traits(pop.trait)->pop.trait
## - euclidean distance
distance<- dist(pop.trait)
## - compare with hyperspace
hyper.quality(distance, hyperspace_euc)->hyperspace_quality_euc     # 0.96
## - gower distance
distance_gow<- gower(pop.trait)
## - compare with hyperspace
hyper.quality(distance_gow, hyperspace_gow)->hyperspace_quality_gow # 0.97

#------------- Compute community hypervolumes (forest functional properties) -----------------
# load remeasured NFI plots
readRDS(here("Data", "NFIs", "Remeasured_minmax_NFI.RDS"))->rem_NFI

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

#Overwrite FR when spp 1 
FDis_3d %>% 
  mutate(FR=case_when(n_spp_plot>1 ~ FR,
                      n_spp_plot<=1 ~ 0))->FDis_3d

#------------ Map forest functional properties ------------------------------------------------
# Divide censuses
FDis_3d %>% filter(census.type == "last_NFI")->FDis_3d_2
FDis_3d %>% filter(census.type == "first_NFI")->FDis_3d_1

# %%%%%%%%%%%%% NFI 1 
# Map centroids First NFI
Map_FD_metric(FDis_3d_1, "C.axis.1", biogeo, legend.name = "LES",  cellsize = 0.5)->c.1_1                      
Map_FD_metric(FDis_3d_1, "C.axis.2", biogeo, legend.name = "SAS",  cellsize = 0.5)->c.2_1

#%%%%%%%%%%%%% NFI 2
# Map centroids Last NFI
Map_FD_metric(FDis_3d_2, "C.axis.1", biogeo, legend.name = "LES", cellsize = 0.5)->c.1_2                      
Map_FD_metric(FDis_3d_2, "C.axis.2", biogeo, legend.name = "SAS",cellsize = 0.5)->c.2_2

# Change axis name
c.1_2$plot+ labs(tag = "A", fill="", x="", y="")+ theme(text = element_text(size = 20))->c.1_2
c.2_2$plot+ labs(tag = "B", fill="", x="", y="")+ theme(text = element_text(size = 20))->c.2_2

# Combine maps together (Fig.3)
cowplot::plot_grid(c.1_2,c.2_2, ncol = 2)->niches_map
niches_map
# Save
pdf(here("Output", "Figures", "Reworked", "Figure_3.pdf"),width=15, height=10)
niches_map
dev.off()

#------------ Calculate community rates of centorid's change -------------------
# Add NFI average interval to compute the rate
NFI_interval<-c(16.74, 9.63, 5.93, 5.36, 5, 10.96, 4.74, 10, 7.50, 5, 5, 11.31, 5, 4)
country<-c("Flanders", "Wallonia", "Czech Republic", "Denmark", "France", "Germany", "Ireland", "Italy", "Netherlands", "Norway", "Poland", "Spain", "Sweden", "Finland")
last_sampling_yr <- c(2016, 2011, 2015, 2021, 2020, 2013, 2022, 2015, 2021, 2021, 2019, 2018, 2017, 1995)
cbind(country,NFI_interval,last_sampling_yr) %>%  as.data.frame() -> interval.nfi
interval.nfi$NFI_interval %>% as.numeric()->interval.nfi$NFI_interval
interval.nfi$last_sampling_yr %>% as.numeric()->interval.nfi$last_sampling_yr

# Calculate niche centroid change over time
FDis_3d %>%  
  dplyr::select(contains("Axis"), census.type, plot.id, country, LAT, LON) %>% 
  pivot_wider(names_from = census.type, values_from = contains("Axis")) %>% left_join(interval.nfi) %>% na.omit() %>% 
  mutate(axis.1.change = (C.axis.1_last_NFI - C.axis.1_first_NFI)/NFI_interval,
         axis.2.change = (C.axis.2_last_NFI - C.axis.2_first_NFI)/NFI_interval) -> Centroids.df_env


# Calculate plot  spp changes over time between census
spp.x.plot %>% 
  pivot_wider(names_from = census.type, values_from = n) %>% 
  rename(spp_NFI_1 = first_NFI, spp_NFI_2 = last_NFI) %>% 
  mutate(spp_change = spp_NFI_2 - spp_NFI_1)->spp_change_df

# Merge centroid change df with spp change df
Centroids.df_env %>%
  left_join(spp_change_df) -> Centroids.df_env

# Plot absolute rate of change on LES
Centroids.df_env %>% filter(spp_NFI_1<=8) %>% 
  ggplot(aes(x=spp_NFI_1, y = (axis.1.change)))+
  geom_jitter(aes(col=as.factor(spp_change)), alpha=0.5)+
  facet_wrap(~ bioreg)+
  geom_hline(yintercept = 0, col="black")+
  scale_colour_discrete_divergingx(palette = "Spectral", na.value = NA)+
  labs(x ="Number of species 1st census", tag="A")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major = element_blank(), text = element_text(size=20))+
  scale_y_continuous(breaks = seq(-0.1, 0.1, by = 0.025), name = "LES change")+
  scale_x_continuous(breaks = seq(1, 12, by = 1)) -> change.1

# Plot absolute rate of change on SAS
Centroids.df_env %>%  filter(spp_NFI_1<=8) %>% 
  ggplot(aes(y = (axis.2.change),x=spp_NFI_1))+
  geom_jitter(aes( col=as.factor(spp_change)), alpha=0.5)+
  facet_wrap(~ bioreg)+  geom_hline(yintercept = 0, col="black")+
  scale_colour_discrete_divergingx(palette = "Spectral", na.value = NA)+
  labs(x ="Number of species 1st census", tag="B", col = "Species change")+
  theme_bw()+
  theme(legend.position=c(.85, .22),  text = element_text(size=20))+
  scale_y_continuous(breaks = seq(-0.1, 0.1, by = 0.025), name = "SAS change")+
  scale_x_continuous(breaks = seq(1, 12, by = 1))+
  theme(panel.grid.major = element_blank())-> change.2

# Combine plots together (Fig.5)
cowplot::plot_grid(change.1, change.2,  nrow = 1)->regional_changes

#save
pdf(here("Output", "Figures", "Reworked", "regional_axis_changes.pdf"),width=15, height=10)
regional_changes
dev.off()

# try to plot latitudinal changes <----------------- provvisorio
# LES
Centroids.df_env %>%  
  ggplot(aes(y = LAT, x = axis.1.change))+
  geom_jitter(aes( col=as.factor(spp_change)), alpha=0.5)+
  #facet_wrap(~ bioreg)+  geom_hline(yintercept = 0, col="black")+
  scale_colour_discrete_divergingx(palette = "PiYG", na.value = NA)+
  labs(x ="LES change", tag="A", col = "Species change")+
  theme_bw()+
  theme(legend.position="none", 
    text = element_text(size=20))+
  scale_y_continuous(name = "Latitude")+
  scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.025))+
  theme(panel.grid.major = element_blank())->LAT_delta_LES

# SAS
Centroids.df_env %>%  
  ggplot(aes(y = LAT, x = axis.2.change))+
  geom_jitter(aes( col=as.factor(spp_change)), alpha=0.5)+
  #facet_wrap(~ bioreg)+  geom_hline(yintercept = 0, col="black")+
  scale_colour_discrete_divergingx(palette = "PiYG", na.value = NA)+
  labs(x ="SAS change", tag="B", col = "Species change")+
  theme_bw()+
  theme(#legend.position=c(.85, .22), 
        text = element_text(size=20))+
  scale_y_continuous(name = "Latitude")+
  scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.025))+
  theme(panel.grid.major = element_blank())->LAT_delta_SAS

# Combine plots together (Fig.5)
cowplot::plot_grid(LAT_delta_LES, LAT_delta_SAS,  nrow = 1, rel_widths = c(0.75,1))->latitudinal_changes

#save
pdf(here("Output", "Figures", "Reworked", "latitudinal_axis_changes.pdf"),width=15, height=10)
latitudinal_changes
dev.off()

#---------------- Load Climatic water balance ----------------------------------
rast(here("Data", "Abiotics", "CWB_2_selfmade", "ERA5", "CWB_monthly.tif"))->CWB

# # Z-transform
mean_r <- app(CWB, mean, na.rm = TRUE)
sd_r   <- app(CWB, sd,   na.rm = TRUE)
z_transformed_CWB <- (CWB - mean_r) / sd_r


# Extract dates from time() or names()
dates_new <- time(z_transformed_CWB)
if (is.null(dates_new)) {
  # tif lost time metadata
  dstr <- sub("^CWB_", "", names(z_transformed_CWB))
  dates_new <- as.Date(paste0(dstr, "_01"), format="%Y_%m_%d")
}
years_new <- as.numeric(format(dates_new, "%Y"))

# Trend calculation
trend_rasters_drought <- calc_trend(z_transformed_CWB, years_new)

# Rename trend layers
names(trend_rasters_drought) <- c("slope", "sd", "std_slope")

# Plot results
plot(trend_rasters_drought[["slope"]], main="Slope of anomaly trend")
plot(trend_rasters_drought[["sd"]], main="Standard deviation of anomalies")
plot(trend_rasters_drought[["std_slope"]], main="Slope normalized by SD")

# Calculate the raster timesries slope of each pixel up to the last plot remeasurment
# Apply to all plots

Centroids.df_env$CWB_std_slope_dynamic_z <- purrr::pmap_dbl(
  list(Centroids.df_env$LON, Centroids.df_env$LAT, Centroids.df_env$last_sampling_yr),
  ~ get_dynamic_pixel_slope(z_transformed_CWB, ..1, ..2, ..3)
)

#---------------- Load disturbances from Viana-Soto & Senf 2024 ----------------
#%%%%%%%%%%%%%%%%%%%%% WIND & BB
read.csv(here("Data", "disturbances_aggregated", "disturbance_0.5deg_barkbeetle-wind.csv"))->dist_wind_bb
dist_wind_bb %>% rename("years" = year)->dist_wind_bb
# Create a grid template
# Take one year subset
d <- subset(dist_wind_bb, years == 2000)

# Define the extent
e <- ext(range(d$long), range(d$lat))

# Define resolution (assuming your grid is regular, e.g. 0.25 degree)
res_x <- min(diff(sort(unique(d$long))))
res_y <- min(diff(sort(unique(d$lat))))

# Example: 0.25 degree grid
r_template <- rast(
  extent = ext(range(dist_wind_bb$long), range(dist_wind_bb$lat)),
  resolution = c(0.5, 0.5),
  crs = "EPSG:4326"
)

# Convert to SpatVector for feeding to Rasterize
v <- vect(d, geom = c("long", "lat"), crs = "EPSG:4326")

# Function to make a raster per year and stack them together (to apply to wind/bb and fire dist)
r_list <- lapply(split(dist_wind_bb, dist_wind_bb$years), function(d) {
  v <- vect(d, geom = c("long", "lat"), crs = "EPSG:4326")
  rasterize(v, r_template, field = "disturbance_rate")})

# Apply function to make rast stack
r_stack_wind <- rast(r_list)
names(r_stack_wind) <- sort(unique(dist_wind_bb$years)) 

# Apply dynamic trends calcs to all plots
Centroids.df_env$WBB_std_slope_dynamic <- purrr::pmap_dbl(
  list(Centroids.df_env$LON, Centroids.df_env$LAT, Centroids.df_env$last_sampling_yr),
  ~get_dynamic_pixel_slope(r_stack_wind, ..1, ..2, ..3)
)

# Extract years
as.numeric(names(r_stack_wind)) -> years

# Apply function
trend_raster_wind <- calc_trend(r_stack_wind, years)

# Inspect results
plot(trend_raster_wind[["slope"]], main="Slope of anomaly trend")
plot(trend_raster_wind[["sd"]], main="Standard deviation of anomalies")
plot(trend_raster_wind[["std_slope"]], main="Slope normalized by SD")

#%%%%%%%%%%%%%%%%%%%%% FIRE
# Load fire data
read.csv(here("Data", "disturbances_aggregated", "disturbance_0.5deg_fire.csv"))->dist_fire
dist_fire %>% rename("years" = year)->dist_fire

# Function to make a raster per year and stack them together
r_list <- lapply(split(dist_fire, dist_fire$years), function(d) {
  v <- vect(d, geom = c("long", "lat"), crs = "EPSG:4326")
  rasterize(v, r_template, field = "disturbance_rate")})

# Apply function
r_stack_fire<- rast(r_list)
names(r_stack_fire) <- sort(unique(dist_fire$years)) 

# Extract years
as.numeric(names(r_stack_fire)) -> years


# Calculate the raster timesries slope of each pizel up to the last plot remeasurment
# Apply dynamic trends calcs to all plots
Centroids.df_env$Fire_std_slope_dynamic <- purrr::pmap_dbl(
  list(Centroids.df_env$LON, Centroids.df_env$LAT, Centroids.df_env$last_sampling_yr),
  ~get_dynamic_pixel_slope(r_stack_fire, ..1, ..2, ..3)
)

# Apply function
trend_raster_fire <- calc_trend(r_stack_fire, years)

# Inspect results
plot(trend_raster_fire[["slope"]], main="Slope of anomaly trend")
plot(trend_raster_fire[["sd"]], main="Standard deviation of anomalies")
plot(trend_raster_fire[["std_slope"]], main="Slope normalized by SD")

# Place 0 where is NA (where there is no trend)
Centroids.df_env[is.na(Centroids.df_env)] <- 0

#--------------------- Plot disturbance trends together ------------------------
# Ensure 'biogeo' is sf
if (!inherits(biogeo, "sf")) stop("biogeo must be an sf object")

# Convert biogeo to the same CRS as your raster
biogeo <- st_transform(biogeo, crs(trend_raster_fire))

# Convert to SpatVector for terra
biogeo_vect <- vect(biogeo)
# --- Crop & mask each raster to Europe
r_fire_eu    <- mask(crop(trend_raster_fire[["std_slope"]], biogeo_vect), biogeo_vect)
r_wbb_eu     <- mask(crop(trend_raster_wind[["std_slope"]], vect(biogeo)), vect(biogeo))
r_drought_eu <- mask(crop(trend_rasters_drought[["std_slope"]], vect(biogeo)), vect(biogeo))

# Make sure extents are clean and zoom tighter to Europe
ext_eu <- ext(-10, 31, 30, 70.2)  # adjust as needed (W, E, S, N)
r_fire_eu <- crop(r_fire_eu, ext_eu)
r_wbb_eu <- crop(r_wbb_eu, ext_eu)
r_drought_eu <- crop(r_drought_eu, ext_eu)

# Convert rasters to data frames for ggplot
df_fire    <- r_to_df(r_fire_eu)
df_wbb     <- r_to_df(r_wbb_eu)
df_drought <- r_to_df(r_drought_eu) 

names(df_fire)[3]    <- "value"
names(df_wbb)[3]     <- "value"
names(df_drought)[3] <- "value"

 # multiply CWB by -1 to have it as measure of drought
df_drought %>% mutate(value = value*(-1))->df_drought

# Make sure your Europe boundary is an sf object
if (!inherits(europe, "sf")) europe <- st_as_sf(europe)

# Create individual maps
p_fire    <- plot_raster_gg(df_fire, "value", biogeo, "Fire", scale_low = "darkgreen", scale_high = "red",
                            res = 0.5, agg_fun = mean, na.rm = TRUE)+ labs(tag = "A", x="", y="Latitude (°N)")
p_wbb     <- plot_raster_gg(df_wbb, "value", biogeo, "Wind & Bark Beetles", scale_low = "darkgreen", scale_high = "red",
                            res = 0.5, agg_fun = mean, na.rm = TRUE)+ labs(tag = "B",  x="Longitude (°E)", y="")
p_drought <- plot_raster_gg(df_drought, "value", biogeo,  "Drought", scale_low = "darkgreen", scale_high = "red",
                            res = 0.5, agg_fun = mean, na.rm = TRUE)+ labs(tag = "C",  x="", y="")

# Combine them with patchwork
combined_plot <- (p_fire + p_wbb + p_drought) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
combined_plot

#save
pdf(here("Output", "Figures", "Reworked", "Disturbance_trends.pdf"), width=15, height=10)
combined_plot
dev.off()

#---------------- Stat Models --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%% Apply GLM Models on axis changes %%%%%%%%%%%%%%%%%%%%%%%%
# formula no interaction c------------------------------------
formula<-c("~ WBB_std_slope_dynamic + CWB_std_slope_dynamic_z + Fire_std_slope_dynamic + (1|country)")

# Apply model on axis 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
drivers_assessment(Centroids.df_env, response_variable = axis.1.change, formula, grouping = TRUE,
                   auto_scale = TRUE)-> Axis_1_DELTA_model

# Model performance & significance# Model performance & significance# Model performance & significance
# performance variables for comparing across formulas (i.e. interaction orno interaction)
Axis_1_DELTA_model$performance 

# Model fit to the data (i.e. R2)
map_df(
  Axis_1_DELTA_model$models,
  ~ as.data.frame(performance::r2(.x)),
  .id = "bioregion"
)

# check model significance % over all plots (its basically by bioreg)
Axis_1_DELTA_model$significance

# Check scaled variables
Axis_1_DELTA_model$scaled_vars 

# Apply model on axis 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
drivers_assessment(Centroids.df_env, response_variable = axis.2.change, formula, grouping = TRUE,
                   auto_scale = TRUE)-> Axis_2_DELTA_model

# Model performance & significance
Axis_2_DELTA_model$performance
Axis_2_DELTA_model$significance

# Check scaled variables
Axis_2_DELTA_model$scaled_vars


# add together all the df %%%%%%%%%%%%%%%%%%%%%%%%%
Axis_1_DELTA_model$model_results %>% 
  rbind(Axis_2_DELTA_model$model_results) %>% 
  mutate(term = case_when(term == "WBB_std_slope_dynamic" ~ "Wind & B.Beetles",
                          term == "Fire_std_slope_dynamic" ~ "Fire",
                          term == "CWB_std_slope_dynamic_z" ~ "Drought"),
         axis = case_when(axis == "axis.1.change" ~ "PC1 LES change",
                          axis == "axis.2.change" ~ "PC2 SAS change"))->DELTA_models

# Produce plots
forest_plot(DELTA_models)->forest_plot_DELTA
forest_plot_DELTA
# Save
pdf(here("Output", "Figures", "Reworked", "Forest_plots_DELTA.pdf"),width=15, height=10)
forest_plot_DELTA
dev.off()

# formula interaction -------------------------------------------------
formula_interaction<-c("~ (WBB_std_slope_dynamic * CWB_std_slope_dynamic_z * Fire_std_slope_dynamic) + (1|country)")

# Apply model on axis 1
drivers_assessment(Centroids.df_env, response_variable = axis.1.change, formula_interaction, grouping = TRUE,
                   auto_scale = TRUE)-> Axis_1_DELTA_model_int
Axis_1_DELTA_model_int$performance

# Apply model on axis 2
drivers_assessment(Centroids.df_env, response_variable = axis.2.change, formula_interaction, grouping = TRUE,
                   auto_scale = TRUE)-> Axis_2_DELTA_model_int
Axis_2_DELTA_model_int$performance


# add together all the df
Axis_1_DELTA_model_int$results %>% 
  rbind(Axis_2_DELTA_model_int$results) %>% 
  mutate(term = case_when(term == "WBB_std_slope_dynamic" ~ "Wind & B.Beetles",
                          term == "Fire_std_slope_dynamic" ~ "Fire",
                          term == "CWB_std_slope_dynamic_z" ~ "Drought",
                          term == "WBB_std_slope_dynamic:CWB_std_slope_dynamic_z" ~ "Wind & B.Beetles : Drought",
                          term == "WBB_std_slope_dynamic:Fire_std_slope_dynamic" ~ "Wind & B.Beetles : Fire",
                          term == "CWB_std_slope_dynamic_z:Fire_std_slope_dynamic" ~ " Drought : Fire",
                          term == "WBB_std_slope_dynamic:CWB_std_slope_dynamic_z:Fire_std_slope_dynamic" ~ "Wind & B.Beetles : Drought : Fire"),
         axis = case_when(axis == "axis.1.change" ~ "PC1 LES change",
                          axis == "axis.2.change" ~ "PC2 SAS change"))->DELTA_models_int


# Produce plots
forest_plot(DELTA_models_int)->forest_plot_DELTA_int
forest_plot_DELTA_int



#%%%%%%%%%%%%%%%%%%%%% Apply Bayesian Models on axis changes %%%%%%%%%%%%%%%%%%%
# formula 
formula<-c(" ~ WBB_std_slope_dynamic +
                   CWB_std_slope_dynamic_z +
                   Fire_std_slope_dynamic +
                   (1 + WBB_std_slope_dynamic + CWB_std_slope_dynamic_z +
                      Fire_std_slope_dynamic | bioreg) +
                   (1 | country)")

drivers_assessment_bayes(Centroids.df_env, response_variable = axis.1.change, formula)-> Axis_1_DELTA_Bayes

drivers_assessment_bayes(Centroids.df_env, response_variable = axis.2.change, formula)-> Axis_2_DELTA_Bayes
