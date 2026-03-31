##################################################################################

# This script contains the functions applied to calculate the niche centroids and produce the 
# main results of the paper "Tree diversity increases functional stability of European forests"
# Corresponding author: Marco Patacca, pataccamarco@gmail.com

#-------------------- scale functions -------------------------------------------
# Function to scale (mean=0, sd=1) returning a vector, rather than built in scale()
my_scale <- function(x){
    return((x-mean(x, na.rm=T))/sd(x,na.rm=T))
  }

# scale a variable from -1 to 1
scale11<-function(x){
  if(sd(x) == 0){
    return(x)
  }else{
    mym <- max(c(max(abs(x[x>0])), max(abs(x[x<0]))))
    x <- x/mym
    return(x)
  }
}

#---------- Principal Components (Fig.2) ---------------------------------------
# These functions are a modified version of the function developed by Maynard et al., 2022 - https://doi.org/10.1038/s41467-022-30888-2 
#----------------------------- PCA axes 1-2 ----------------
ggpca_weighted_1.2 <- function(tr_use, trait_vars, num_show, alpha = 0.1, flip_coord = FALSE, arrow_head = 6, 
                               flip_x = 1, flip_y = 1, use_vars = NULL, label_size = 5, show_labs = TRUE, alpha_const = 1,
                               alpha_all = 0, line_wd_big = 1.5, all_alpha = 0.2, angio_list = NULL, flip_z = 1, flip_zz = 1,
                               show_legend = TRUE, title = ""){
  
  # set axes
  axes = c(1,2)
  # browser()
  # set seed to ensure coordinates are consistent
  set.seed(10)
  
  # flip the coordinates if needed
  if(flip_coord){
    tmp <- flip_x
    flip_x <- flip_y
    flip_y <- tmp
  }
  
  # keep only the unique observations, and add in some noise to allow for merge later on
  tr_unique <- tr_use %>% dplyr::select(wt, all_of(trait_vars)) %>%
    mutate(rowid = 1:nrow(.)) %>% 
    gather(trait, value, -wt, -rowid) %>% 
    rowwise() %>% 
    mutate(value = value+runif(1, -1e-12,1e-12)) %>%
    ungroup %>% 
    spread(trait, value) %>% arrange(rowid) %>% dplyr::select(-rowid)
  
  # keep only the traits
  tr_use <- tr_use %>% dplyr::select(-all_of(trait_vars), -wt) %>% bind_cols(tr_unique)
  
  # get the trait matrix and weights
  x_pca <- tr_unique %>% dplyr::select(-wt) %>% as.matrix()
  w_pca <- tr_unique$wt
  
  # run the weighted pca
  pca_use <- wpca(x = x_pca,  w = w_pca, center=TRUE, scale=TRUE)
  
  # get the mean value to scale the eigenvectors to mean unit variance
  sdmean <- sqrt(mean(pca_use$d^2))
  
  # extract the objects, scale the eigenvalues and x matrix
  sdev <- pca_use$d/sdmean
  rotmat <- t(pca_use$vt)
  xmat <- pca_use$pc/sdmean
  
  # flip the rotations to align the plots
  if(flip_x == -1){
    xmat[,1] <- -xmat[,1]
    rotmat[,1] <- -rotmat[,1]
  }
  if(flip_y == -1){
    xmat[,2] <- -xmat[,2]
    rotmat[,2] <- -rotmat[,2]		
  }
  if(flip_z == -1){
    xmat[,3] <- -xmat[,3]
    rotmat[,3] <- -rotmat[,3]		
  }
  if(flip_zz == -1){
    xmat[,4] <- -xmat[,4]
    rotmat[,4] <- -rotmat[,4]		
  }
  
  # add in names
  rownames(rotmat) <- colnames(x_pca)
  colnames(rotmat) <- paste0("PC",1:ncol(rotmat))
  colnames(xmat) <- paste0("PC",1:ncol(xmat))
  
  if(num_show<0){
    num_show <- nrow(rotmat)
  }
  
  # make sure we aren't displaying more traits than we have
  num_show <- min(num_show, nrow(rotmat))
  
  # get the percentage explained, for plotting on the axes
  percvar <- round((sdev^2)/sum(sdev^2), 3)*100
  
  # get the dataframe for plotting points
  pto <- tr_use %>% left_join(bind_cols(tr_unique %>% dplyr::select(-wt), xmat %>% data.frame() %>% as_tibble() %>% dplyr::select(contains("PC", ignore.case = FALSE))), by = trait_vars) %>% mutate(genus = word(accepted_bin, 1)) 
  
  # scale the axes to allow direct overlay of the vectors
  pt_scaled <- pto %>% dplyr::select(-wt) %>% left_join(angio_list) %>% 
    gather(pc, value, -accepted_bin, -angio, -genus, -group ,-LAT, -LON) %>% 
    group_by(pc) %>% 
    mutate(value = ifelse(grepl("PC", pc, ignore.case = FALSE), scale11(value), value)) %>% ungroup %>% 
    spread(pc, value) %>% dplyr::select(-LAT, -LON)
  
  # the final dataset for plotting, with angio vs. gymno added
  pt_use <-  pt_scaled %>% 
    rename(Genus = genus) %>% mutate(angio = ifelse(angio == 1, "Angiosperm", "Gymnosperm"))
  
  # get the loadings
  loads <- loads_orig <- rotmat%*%diag(sdev) %>% data.frame() %>%
    setNames(paste0("PC",1:ncol(.))) %>% mutate(trait = rownames(.)) %>% 
    as_tibble() %>% dplyr::select(trait, names(.)) 
  
  # get the data frame for the arrows
  arrow_df <- loads %>% dplyr::select(trait, paste0("PC",axes)) %>% setNames(c("trait","x","y")) %>% mutate(Genus = NA) 
  
  # if we haven't specified which variables to plot calculate the ones with the highest loadings
  if(is.null(use_vars)){
    # variable closest to the x or y axis
    best <- arrow_df %>% rowwise() %>%
      mutate(y_rat = min(dist(rbind(c(x,y), c(0,1))), dist(rbind(c(x,y), c(0,-1))))) %>% 
      mutate(x_rat = min(dist(rbind(c(x,y), c(1,0))), dist(rbind(c(x,y), c(-1,0))))) %>% 
      arrange(y_rat) %>% 
      add_column(axis1 = c(rep(1, num_show), rep(0, nrow(.) - num_show))) %>% 
      arrange(x_rat) %>% 
      add_column(axis2 = c(rep(2, num_show), rep(0, nrow(.) - num_show))) %>% 
      rowwise() %>% 
      mutate(axis = ifelse((axis1 + axis2) == 0, 0, ifelse((axis1 + axis2)%in%c(1,2), max(axis1, axis2), min(axis1, axis2)))) %>% 
      ungroup %>% 
      dplyr::select(axis, trait)
  }else{
    best <- tibble(trait = use_vars, axis = axes[1])
  }
  
  # get the colors
  my_color <- my_fill <- c(redmonder.pal(8,"qPBI")[c(1,3, 5,8)],redmonder.pal(8,"qMSOMed")[c(4)])
  my_shp <- c(15,16,17,18,19)
  
  # create a copy for plotting
  plot_dt <- pt_use 
  
  # scale the transparency
  alpha_scale <- plot_dt  %>% group_by(Genus) %>% tally() %>% ungroup %>% mutate(n = n%/%1000) %>% mutate(al = 0.1) %>% dplyr::select(Genus, al) %>% deframe

  # create the base plot                       
  g1 <- ggplot(
    data = plot_dt %>%
      filter(complete.cases(.), abs(PC1) < 1, abs(PC2) < 1),  # Remove missing & extreme values
    aes(x = PC1, y = PC2)) +
    geom_point(aes(color = angio, shape = angio, alpha = angio), size = 0.35) +  
    theme_bw()
  
  # create the transparancy value based on number of observations (note that Genus is set to angio vs. gymno in the main figure)
  alv <- plot_dt %>% dplyr::select(angio, PC1, PC2) %>% filter(complete.cases(.), abs(PC1)< 1, abs(PC2)<1) %>% rename(Genus = angio) %>% 
    group_by(Genus) %>% tally() %>% ungroup %>%
    rowwise() %>% mutate(n = max(n, 15001)) %>% ungroup %>% 
    mutate(al = 1/(n%/%1000)/alpha_const) %>% arrange(Genus) %>% dplyr::select(al) %>% unlist() %>% as.numeric()
  
  # plotting colors/scales
  my_color <- darken(c("darkorange3", "darkgreen"), 0.3)
  my_fill <- lighten(c("darkorange3", "darkgreen"), 0)
  my_shp <- c(22, 24) #my_shp[c(1,2)]
  alpha_scale <- alv
  names(my_color) <- names(my_fill) <- names(my_shp) <- c("Angiosperm", "Gymnosperm")
  
  # set the transparency depending on if we're plotting both angio and gymno, or just one
  if(length(alpha_scale)==1){
    names(alpha_scale) <- unique(plot_dt$angio)
  }else{
    names(alpha_scale) <- c("Angiosperm", "Gymnosperm")
  }
  
  if(num_show < nrow(rotmat)){
    #plot
    g1 <- g1 + 
      geom_hline(yintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
      geom_vline(xintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
      xlab(paste0("PC Axis ",axes[1]," (",percvar[axes[1]],"%)"))+
      ylab(paste0("PC Axis ",axes[2]," (",percvar[axes[2]],"%)"))+
      geom_segment(data = arrow_df %>% left_join(best, by = "trait") %>% filter(axis%in%axes), aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(arrow_head, "pt")), size = line_wd_big, color = "black") + 
      geom_text_repel(data = arrow_df %>% left_join(best, by = "trait") %>% filter(axis%in%axes),  aes(x = x, y = y, label = trait), fontface = "bold", box.padding = 0.75, segment.size = 0, size = label_size*1)+
      geom_segment(data = arrow_df %>% left_join(best, by = "trait") %>% filter(!axis%in%axes), aes(x = 0, y = 0, xend = x, yend = y), size = 0.7, linetype = 1, color  = "gray40", alpha = 0.4)+
      scale_fill_manual(values = my_fill) +
      scale_color_manual(values = my_color) +
      scale_shape_manual(values = my_shp) +
      scale_alpha_manual(values = alpha_scale) +
      theme_classic()
  }else{
    g1 <- g1 + 
      geom_hline(yintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
      geom_vline(xintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
      xlab(paste0("PC Axis ",axes[1]," (",percvar[axes[1]],"%)"))+
      ylab(paste0("PC Axis ",axes[2]," (",percvar[axes[2]],"%)"))+
      geom_segment(data = arrow_df %>% left_join(best, by = "trait"), aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(arrow_head, "pt")), size = line_wd_big, color = "black") + 
      geom_text_repel(data = arrow_df %>% left_join(best, by = "trait"),  aes(x = x, y = y, label = trait), fontface = "bold", box.padding = 0.75, segment.size = 0, size = label_size*1)+
      geom_segment(data = arrow_df %>% left_join(best, by = "trait"), aes(x = 0, y = 0, xend = x, yend = y), size = 0.7, linetype = 1, color  = "gray40", alpha = 0.4)+
      scale_fill_manual(values = my_fill) +
      scale_color_manual(values = my_color) +
      scale_shape_manual(values = my_shp) +
      scale_alpha_manual(values = alpha_scale) + 
      theme_classic()
  }
  # add in a title if desired
  if(!is.null(title)){
    g1 <- g1 + ggtitle(title)
  }
  
  # add the legend
  if(!show_legend){
    g1 <- g1 + theme(legend.position = "none", plot.title = element_text(face = "bold"))
  }else{
    g1 <- g1 + theme(legend.position = c(0.8,0.1), 
                     legend.title = element_blank(), 
                     legend.text = element_text(size = 11),
                     legend.background = element_rect(colour = 'black', fill = NA, linetype='solid', size = 0.2),
                     legend.spacing.y = unit(0, "mm"))+
      guides(fill = guide_legend(override.aes = list(size = 2)), alpha = FALSE)
    
  }
  
  # flip the coordinates to align the figures
  if(flip_coord){
    g1 <- g1 + coord_flip()
  }
  
  # plot
  show(g1)
  
  return(list(plot = g1, vars = best %>% filter(axis%in%axes) %>% dplyr::select(trait) %>% unlist(), pca = pca_use, tr = pto, loadings = loads_orig, axis_perc = sdev^2/sum(sdev^2)))
}


#---------- Build hyperspace ---------------------------------------------------
#The hyperspace needs to be built with all the species together once as is a measure of distance across species.
#Because the traits of my species differ within the same species (variation due to the environmental covariates in the imputation)
#here each plot represents 1 population of the same species, which in the hyperspace will be treated as a separate species
hyperspace <- function(NFI_df, distance_method){
  
# Column names of NFI_df (accepted_bin is the species name):
  # [1] "country"      "LAT"          "LON"          "plot.id"      "census.n"     "plot.obs.id"  "ba"           "accepted_bin" "trait1"       "trait2"       "trait3"       "trait4"      
  # [13] "trait5"       "trait6"       "trait7"       "trait8"       "trait14"      "trait15"      "trait37"      "trait44"      "trait46"      "trait48"      "trait66"      "trait185"    
  # [25] "trait281"     "trait324"     "trait719"     "trait773"     "trait1312"    "trait3110"    "trait3117"    "trait24"      "trait45"      "uniqueID"     "plot.ba"      "rel.ba" 
  
  # 1) BUILD HYPERSPACE   
  #dplyr::select species and all 25 traits
  NFI_df %>% 
    dplyr::select(uniqueID, accepted_bin, contains("trait"))->species.traits
  
  #Summarise spp per plot (average)
  species.traits %>% 
    group_by(accepted_bin) %>%
    summarise(across(where(is.numeric), mean)) %>% 
    as.data.frame() %>% 
    ungroup()->pop.trait
  
  #make it a matrix with spp names on rows
  matrix.numeric.traits(pop.trait)->pop.trait
  
  #set the seed 
  set.seed(1312)
  #build hyper space (here no weight is needed as we are using species averages per trait - i.e. no unbalance of observations)
  hyper.build(trait = pop.trait, distance = distance_method, 
              weight = NULL, axes = 3, convert = NULL) -> hyper.space
  
  return(hyper.space)
}

#---------- Build hypervolume and extract its centroids ------------------------
Niche_calc<- function(NFI_df, hyperspace) {
  
  # 1) Build spp abundance matrix
  NFI_df %>% 
    dplyr::select(uniqueID, accepted_bin, rel.ba) %>% 
    unique() %>% 
    pivot_wider(names_from= accepted_bin, values_from = rel.ba, values_fn = {mean}) %>% 
    as.data.frame()-> community_df
  
  #substitute Nas
  community_df[is.na(community_df)]<-0
  
  #set matrix names on rows (the matrix has to be only numeric)
  matrix.numeric.comm(community_df) ->community_df
  
  
  # 2) Build the df for final results
  # Some functional diversity metrics are commented out as are not needed for replicating the manuscript results
  NFI_df %>%
    dplyr::select(uniqueID) %>% unique() %>%
    mutate(FDis = rep(NA),                    
           FR = rep(NA),
           #FEv =rep(NA),
           C.axis.1 = rep(NA),
           C.axis.2 = rep(NA),
           C.axis.3 = rep(NA))->FDis_results
  
  # 3) Build the hypervolume of the community 
  
  #build kernel hypervolumes
  set.seed(1312) %>% 
    kernel.build(
      comm = community_df,
      trait = hyperspace,
      distance = "gower",         # keep the same distance measure as the one used to build the hyperspace
      method.hv = "gaussian",
      abund = FALSE,
      weight = NULL,
      axes = 3,
      convert = NULL,
      cores = 0) -> hype
  
  
  # asses functional dispersion
  set.seed(1312)
  for (i in 1:length(hype@HVList)){
    ifelse(hype@HVList[[i]]@Dimensionality != 1,  kernel.dispersion(hype@HVList[[i]], func = "divergence", frac =0.1 ), 0 )} -> FDis_results[i,2]

  # asses functional richness
  set.seed(1312)
  for (i in 1:length(hype@HVList)){
    kernel.alpha(hype@HVList[[i]])} -> FDis_results[i,3]

  # #assess functional eveness
  # set.seed(1312)
  # for (i in 1:length(hype@HVList)){
  #   ifelse(hype@HVList[[i]]@Dimensionality != 1, kernel.evenness(hype@HVList[[i]]),0)} -> FDis_results[i,4]
  
  #get community centroid value
  set.seed(1312)
  for (i in 1:length(hype@HVList)){
    get_centroid(hype@HVList[[i]]) %>% as.list() -> FDis_results[i,4:6]
  }
  #5) dplyr::select and print
  return(list(hype, FDis_results)) 
  
  
}

#---------- Function to make maps (Fig.3) --------------------------------------
Map_FD_metric <- function(Data, variable_to_map, biogeo, legend.name,
                          cellsize = 0.2, raster_filename = NULL) {

  # 1. Prepare points as sf
  points_sf <- Data %>%
    dplyr::select({{ variable_to_map }}, LON, LAT) %>%
    rename(x = {{ variable_to_map }}) %>%
    st_as_sf(coords = c("LON", "LAT"), crs = 4326)

  # Align CRS
  if (st_crs(biogeo) != st_crs(points_sf)) {
    biogeo <- st_transform(biogeo, st_crs(points_sf))
  }

  # 2. Create grid and compute mean value per cell
  grid_sf <- st_make_grid(points_sf, cellsize = cellsize, square = TRUE) %>%
    st_sf(geometry = .) %>%
    mutate(cell_id = row_number())

  shares_df <- st_join(points_sf, grid_sf, join = st_within) %>%
    st_drop_geometry() %>%
    group_by(cell_id) %>%
    summarise(x = mean(x, na.rm = TRUE), .groups = "drop")

  grid_shares_sf <- left_join(grid_sf, shares_df, by = "cell_id")

  # 3. Convert to raster properly
  # Use terra's rasterization rather than centroids
  rast_template <- rast(ext(grid_shares_sf), resolution = cellsize, crs = "EPSG:4326")
  rast_obj <- rasterize(vect(grid_shares_sf), rast_template, field = "x")

  # Optionally save raster
  if (!is.null(raster_filename)) {
    writeRaster(rast_obj, raster_filename, overwrite = TRUE)
  }

  # 4. Build the ggplot map
  plot_obj <- ggplot() +
    geom_sf(data = grid_shares_sf, aes(fill = x), color = NA) +
    geom_sf(data = biogeo, fill = NA, color = "black", linewidth = 0.4) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = 2,
      na.value = "white",
      name = legend.name,
      limits = c(-1, 1) * max(abs(grid_shares_sf$x), na.rm = TRUE)
    ) +
    coord_sf(xlim = c(-10, 35), ylim = c(36, 71), expand = FALSE) +
    scale_x_continuous(breaks = seq(-10, 35, 10)) +
    scale_y_continuous(breaks = seq(35, 70, 5)) +
    labs(x = "Longitude (E)", y = "Latitude (N)") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 12)
    )

  #-------------------------------------------------
  # 5. Return both plot and raster
  #-------------------------------------------------
  return(list(plot = plot_obj, raster = rast_obj))
}


#---------- Extract disturbance trends -----------------------------------------
#Function to compute slope, SD, and slope/SD for each pixel throughout all the 
# study period (i.e. full raster stack) for wind/bb and fire data
calc_trend <- function(r_stack, years) {
  stopifnot(nlyr(r_stack) == length(years))  # check
  
  # Regression slope per pixel
  slope_fun <- function(y) {
    if (all(is.na(y))) return(c(NA, NA, NA))
    
    x <- years
    fit <- lm(y ~ x)
    slope <- coef(fit)[2]
    sd_val <- sd(y, na.rm = TRUE)
    std_slope <- ifelse(sd_val == 0, NA, slope / sd_val)
    c(slope, sd_val, std_slope)
  }
  
  # Apply to raster stack
  r_out <- app(r_stack, slope_fun)
  names(r_out) <- c("slope", "sd", "std_slope")
  return(r_out)
}

# Function to calculate slope, sd, and standardized slope for each pixel 
# This function is dynamic per plot (unlikely calc_slope_sd)
# i.e. calculate the raster timesries slope of each pizel up to the last plot remeasurment
get_dynamic_pixel_slope <- function(r_stack, lon, lat, sampling_year) {
  
  # Extract numerical years from layer names
  years <- as.numeric(gsub("_.*$", "", names(r_stack)))
  
  # Only valid sampling years
  if (is.na(sampling_year)) return(NA_real_)
  
  # select time window up to sampling year
  valid_layers <- which(years <= sampling_year)
  if (length(valid_layers) < 3) return(NA_real_)
  
  # Build proper terra point with geometry specification
  pt <- vect(data.frame(x = lon, y = lat), geom = c("x", "y"), crs = crs(r_stack))
  
  # Extract pixel time series
  vals <- terra::extract(r_stack[[valid_layers]], pt)[1, -1] |> unlist()
  if (all(is.na(vals))) return(NA_real_)
  
  # Linear regression slope
  fit <- lm(vals ~ years[valid_layers])
  slope <- coef(fit)[2]
  sd_val <- sd(vals, na.rm = TRUE)
  std_slope <- (slope / sd_val)
  
  return(std_slope)
}

#--------------- Make biogeographical region maps ------------------------------
plot_region <- function(sf_data, region_code) {
  ggplot() +
    geom_sf(data = sf_data, fill = "grey80", color = "black") +  # Full Europe in grey
    geom_sf(data = sf_data %>% filter(code == region_code), 
            aes(fill = code), color = "black") +  # Highlight selected region
    scale_fill_manual(values = color_mapping[region_code]) +  # Highlight color (change as needed)
    theme_void()+
    theme(legend.position = "none") -> figure
  #save
  pdf(here("Output", "Figures","Reworked", paste0(region_code, ".pdf")), width=15, height=10) 
  print(figure)
  dev.off()
}


# --------- Convert rasters to data frames for ggplot -------------------------
r_to_df <- function(r) {
  as.data.frame(r, xy = TRUE, na.rm = TRUE)
}

# --- Define a plotting function for the disturbance trend mapping
plot_raster_gg <- function(df, value, biogeo, title,
                           scale_low = "darkgreen", scale_high = "red",
                           res = 0.5,
                           agg_fun = mean, na.rm = TRUE) {

  # Check data structure
  if (!all(c("x", "y", value) %in% names(df))) {
    stop("Dataframe must contain columns 'x', 'y', and the specified value column.")
  }
  
  # Ensure biogeo is sf and CRS matches raster
  if (!inherits(biogeo, "sf")) {
    stop("'biogeo' must be an sf object (multipolygon geometry).")
  }
  
  # Assuming raster df is in lon/lat (WGS84)
  raster_crs <- st_crs(4326)
  
  # Align CRS
  if (st_crs(biogeo) != raster_crs) {
    biogeo <- st_transform(biogeo, raster_crs)
  }
  
  # Clean invalid geometries (important for multipolygons)
  biogeo <- suppressWarnings(st_make_valid(biogeo))
  
  # ---- SNAP TO A FIXED 0.5° GRID ----
  # Snap x/y to nearest grid center at 'res' degrees
  df_fixed <- df %>%
    mutate(
      x = round(x / res) * res,
      y = round(y / res) * res
    ) %>%
    group_by(x, y) %>%
    summarise(
      value_fixed = agg_fun(.data[[value]], na.rm = na.rm),
      .groups = "drop"
    )
  
  # Core map
  ggplot() +
    # Raster layer with explicit cell size
    geom_tile(
      data = df_fixed,
      aes(x = x, y = y, fill = value_fixed),
      width = res, height = res
    ) +
    # Biogeographical region outlines
    geom_sf(data = biogeo, fill = NA, color = "black", linewidth = 0.4) +
    # Continuous diverging fill scale
    scale_fill_gradient2(
      low = scale_low,
      mid = "white",
      high = scale_high,
      name = "",
      midpoint = 0,
      na.value = "grey90"
    ) +
    # Geographic frame and axis labels
    coord_sf(xlim = c(-10, 35), ylim = c(36, 71), expand = FALSE) +
    scale_x_continuous(breaks = seq(-10, 35, 10)) +
    scale_y_continuous(breaks = seq(35, 70, 5)) +
    labs(
      x = "Longitude (E)",
      y = "Latitude (N)",
      title = title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10)
    )
}


#--------------- Statistical analyses ------------------------------------------
# Function that test a glm model on a response variable (in this case the 
#community centroids position on the axis, by biogeographical region), 
# The predictive variables tested are extreme drought, wind +barkbeetles and fire frequencies 

drivers_assessment <- function(data = Centroids.df_env, response_variable, formula, grouping = FALSE,
                               auto_scale = TRUE)     # If true will check variables if rescaling is needed
                               {

  # --------------------------------------------------------------------
  # 0. Drop geometry if present
  # --------------------------------------------------------------------
  if ("geometry" %in% names(data)) {
    data <- sf::st_drop_geometry(data)
  }
  
  # --------------------------------------------------------------------
  # 1. Response variable and formula construction
  # --------------------------------------------------------------------
  response_var <- deparse(substitute(response_variable))
  model_formula <- as.formula(paste(response_var, formula))
  
  # --------------------------------------------------------------------
  # 2. Extract predictor names safely (fixed effects only)
  # --------------------------------------------------------------------
  fixed_formula <- lme4::nobars(model_formula)
  predictor_terms <- attr(terms(fixed_formula), "term.labels")
  
  predictors <- predictor_terms %>%
    gsub("\\(.*?\\)", "", .) %>%    
    gsub("[:*]", " ", .) %>%        
    trimws() %>%
    strsplit(" ") %>% 
    unlist() %>% 
    unique()
  
  # keep predictors that exist in dataset
  predictors <- predictors[predictors %in% names(data)]
  
  # only numeric predictors
  numeric_predictors <- predictors[sapply(data[predictors], is.numeric)]
  
  # --------------------------------------------------------------------
  # 3. Auto-scale numeric predictors except z-scores
  # --------------------------------------------------------------------
  scaled_vars <- data.frame(variable = character(),
                            mean_before = numeric(),
                            sd_before = numeric(),
                            scaled = logical())
  
  if (auto_scale && length(numeric_predictors) > 0) {
    for (v in numeric_predictors) {
      m <- mean(data[[v]], na.rm = TRUE)
      s <- sd(data[[v]], na.rm = TRUE)
      
      # detect z-score
      is_z <- abs(m) < 0.05 & abs(s - 1) < 0.05
      
      if (!is_z) {
        data[[v]] <- scale(data[[v]])[, 1]
        scaled <- TRUE
      } else {
        scaled <- FALSE
      }
      
      scaled_vars <- rbind(
        scaled_vars,
        data.frame(variable = v, mean_before = m, sd_before = s, scaled = scaled)
      )
    }
  }
  
  # --------------------------------------------------------------------
  # 4. Fit models (grouped or not)
  # --------------------------------------------------------------------
  if (grouping) {
    model_list <- data %>%
      group_by(bioreg) %>%
      group_split() %>%
      map(~ lmerTest::lmer(model_formula, data = .x)) %>%
      set_names(unique(data$bioreg))
    
  } else {
    model <- lmerTest::lmer(model_formula, data = data)
    model_list <- list(All_Data = model)
  }
  
  # --------------------------------------------------------------------
  # 5. Model performance (AIC, BIC, logLik)
  # --------------------------------------------------------------------
  performance <- model_list %>%
    imap_dfr(~{
      data.frame(
        group = .y,
        AIC = AIC(.x),
        BIC = BIC(.x),
        logLik = as.numeric(logLik(.x))
      )
    })
  
  # --------------------------------------------------------------------
  # 6. Extract coefficient results
  # --------------------------------------------------------------------
  results <- model_list %>%
    map_dfr(~ broom.mixed::tidy(.x, conf.int = TRUE, p.value = TRUE), .id = "group") %>%
    filter(term %in% predictors) %>%
    mutate(axis = response_var)
  
  # --------------------------------------------------------------------
  # 7. Overall significance (plots affected per driver)
  # --------------------------------------------------------------------
  region_sig <- results %>%
    mutate(
      significant = p.value < 0.05,
      direction = ifelse(estimate > 0, "positive", "negative")
    ) %>%
    select(group, term, significant, direction) %>%
    rename(bioreg = group)
  
  plots_per_region <- data %>% count(bioreg, name = "n_plots")
  
  region_sig_with_plots <- region_sig %>%
    left_join(plots_per_region, by = "bioreg")
  
  significance_overall <- region_sig_with_plots %>%
    group_by(term) %>%
    summarise(
      total_plots = sum(n_plots),
      sig_plots = sum(ifelse(significant, n_plots, 0)),
      percent_sig = 100 * sig_plots / total_plots,
      
      pos_plots = sum(ifelse(significant & direction == "positive", n_plots, 0)),
      neg_plots = sum(ifelse(significant & direction == "negative", n_plots, 0)),
      percent_pos = 100 * pos_plots / total_plots,
      percent_neg = 100 * neg_plots / total_plots,
      .groups = "drop"
    )
  
  # --------------------------------------------------------------------
  # 8. Return all results INCLUDING MODELS
  # --------------------------------------------------------------------
  return(list(
    model_results = results,
    performance = performance,
    scaled_vars = scaled_vars,
    significance = significance_overall,
    models = model_list      
  ))
}


# The same concept but with bayesian models. This is because contries are too little within bioregion to explain random variatin.
# Bayesian models can deal with it because they have posterior distributon to sample from.

drivers_assessment_bayes <- function(
    data,
    response_variable,
    formula,
    auto_scale = TRUE
) {
  # --------------------------------------------------------------------
  # 0. Drop geometry if present
  # --------------------------------------------------------------------
  if ("geometry" %in% names(data)) {
    data <- sf::st_drop_geometry(data)
  }
  
  # --------------------------------------------------------------------
  # 1. Response variable and model formula
  # --------------------------------------------------------------------
  response_var  <- deparse(substitute(response_variable))
  model_formula <- as.formula(paste(response_var, formula))
  
  # --------------------------------------------------------------------
  # 2. Extract predictor names (only fixed effects)
  # --------------------------------------------------------------------
  fixed_formula   <- lme4::nobars(model_formula)  # remove (1 | country) if present
  predictor_terms <- attr(terms(fixed_formula), "term.labels")
  
  predictors <- predictor_terms %>%
    gsub("\\(.*?\\)", "", .) %>%   # remove functions
    gsub("[:*]", " ", .) %>%       # remove interactions
    trimws() %>%
    strsplit(" ") %>%
    unlist() %>%
    unique()
  
  # Keep only predictors that exist in data
  predictors <- predictors[predictors %in% names(data)]
  
  # Only numeric predictors
  numeric_predictors <- predictors[sapply(data[predictors], is.numeric)]
  
  # --------------------------------------------------------------------
  # 3. Auto-scale numeric predictors except already z-scored vars
  # --------------------------------------------------------------------
  scaled_vars <- data.frame(
    variable = character(),
    mean_before = numeric(),
    sd_before = numeric(),
    scaled = logical(),
    stringsAsFactors = FALSE
  )
  
  if (auto_scale && length(numeric_predictors) > 0) {
    for (v in numeric_predictors) {
      m <- mean(data[[v]], na.rm = TRUE)
      s <- sd(data[[v]], na.rm = TRUE)
      is_z <- abs(m) < 0.05 & abs(s - 1) < 0.05  # detect z-score
      
      if (!is_z) {
        data[[v]] <- scale(data[[v]])[, 1]
        scaled <- TRUE
      } else {
        scaled <- FALSE
      }
      
      scaled_vars <- rbind(
        scaled_vars,
        data.frame(variable = v, mean_before = m, sd_before = s, scaled = scaled)
      )
    }
  }
  
  # --------------------------------------------------------------------
  # 4. Fit one Bayesian model per bioreg
  # --------------------------------------------------------------------
  model_list <- data %>%
    group_by(bioreg) %>%
    group_split() %>%
    map(~ brms::brm(
      model_formula,
      data = .x,
      family = gaussian(),
      chains = 4,
      cores = 4,
      refresh = 0       # suppress output
    )) %>%
    set_names(unique(data$bioreg))
  
  # --------------------------------------------------------------------
  # 5. Model performance (LOOIC + Bayes R²)
  # --------------------------------------------------------------------
  performance <- imap_dfr(model_list, ~ {
    loo_est <- brms::loo(.x)
    r2_est  <- brms::bayes_R2(.x)
    
    data.frame(
      bioreg      = .y,
      LOOIC       = loo_est$estimates["looic", "Estimate"],
      SE_LOOIC    = loo_est$estimates["looic", "SE"],
      Bayes_R2    = mean(r2_est),
      stringsAsFactors = FALSE
    )
  })
  
  # --------------------------------------------------------------------
  # 6. Extract model results (credible intervals, no p-values)
  # --------------------------------------------------------------------
  results <- model_list %>%
    map_dfr(~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE), .id = "bioreg") %>%
    filter(term %in% predictors) %>%
    mutate(axis = response_var) %>%
    mutate(
      significant = case_when(
        conf.low > 0 ~ TRUE,
        conf.high < 0 ~ TRUE,
        TRUE ~ FALSE
      ),
      direction = ifelse(estimate > 0, "positive", "negative")
    )
  
  # --------------------------------------------------------------------
  # 7. Significance (% of total plots affected per driver)
  # --------------------------------------------------------------------
  plots_per_bioreg <- data %>% count(bioreg, name = "n_plots")
  
  significance <- results %>%
    select(bioreg, term, significant, direction) %>%
    left_join(plots_per_bioreg, by = "bioreg") %>%
    group_by(term) %>%
    summarise(
      total_plots = sum(n_plots),
      sig_plots = sum(ifelse(significant, n_plots, 0)),
      percent_sig = 100 * sig_plots / total_plots,
      pos_plots = sum(ifelse(significant & direction == "positive", n_plots, 0)),
      neg_plots = sum(ifelse(significant & direction == "negative", n_plots, 0)),
      percent_pos = 100 * pos_plots / total_plots,
      percent_neg = 100 * neg_plots / total_plots,
      .groups = "drop"
    )
  
  # --------------------------------------------------------------------
  # 8. Return everything
  # --------------------------------------------------------------------
  return(list(
    model_results = results,
    performance   = performance,
    scaled_vars   = scaled_vars,
    significance  = significance,
    models        = model_list
  ))
}


# Function that create the forest plot from the model results produced by the 
#previous function
forest_plot<- function(model_results){
  model_results %>% 
    ggplot(aes(x = estimate, y = term, fill = group)) +
    geom_point(position = position_dodge(width = 0.5), # Plot the point estimates (coefficients)
               size = 5, shape = 21, color = "black") +  # filled circle with black border
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), # Add horizontal error bars (CI)
                   position = position_dodge(width = 0.5), 
                   height = 0.2, color = "black") +        # black error bars
    geom_text(aes(label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "")),,# Add p-values only if significant
      nudge_y = 0.2, size = 5) + 
    geom_vline(xintercept = 0, col = "black") + # vertical line on 0
    facet_grid(axis ~ group, scales = "free") + # faceting with 2 terms
    scale_fill_manual(values = rev(hcl.colors(5, "Roma"))) +  # fill dots with biome colours
    labs(
      x = "Effect size", 
      y = "") +
    theme_minimal() +
    theme(legend.position = "none",
          text = element_text(size = 20),
          strip.background = element_blank(),
          strip.text= element_text(),
          panel.spacing.x = unit(1, "lines"),
          axis.text.y = element_text(),  # Keep for now
          axis.ticks.y = element_line(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) ->plot        # Remove y-axis ticks from all panels
  
  return(plot)
}

