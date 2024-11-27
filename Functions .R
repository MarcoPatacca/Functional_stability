##################################################################################

# This script contains the functions applied to calculate the niche centroids and produce the 
# main results of the peper "Tree diversity increases functional stability of European forests"
# Corresponding author: Marco Patacca, pataccamarco@gmail.com

#-------------------- scale function -------------------------------------------
# Function to scale (mean=0, sd=1) returning a vector, rather than built in scale()
my_scale <- function(x){
    return((x-mean(x, na.rm=T))/sd(x,na.rm=T))
  }

#---------- Principal Components (Fig.2) ---------------------------------------
# These functions are a modified version of the function developed by Maynard et al., 2022 - https://doi.org/10.1038/s41467-022-30888-2 
#----------------------------- PCA axes 1-2 ----------------
ggpca_weighted_1.2 <- function(tr_use, trait_vars, num_show, alpha = 0.1, flip_coord = FALSE, arrow_head = 6, 
                               flip_x = 1, flip_y = 1, use_vars = NULL, label_size = 5, show_labs = TRUE, alpha_const = 1,
                               alpha_all = 0, line_wd_big = 1.5, all_alpha = 0.2, angio_list = NULL, flip_z = 1, flip_zz = 1,
                               show_legend = TRUE, title = ""){
  
  
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
  
  # make sure wearen't displaying more traits than we have
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
  
  # if we haven't specified which variables to plot (Fig 2a), calculate the ones with the highest loadings
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
  
  # create the base plot                          NB : here change manually the axis of points coord. maybe discuss w dan to streamline?
  g1 <- ggplot(data = plot_dt %>% 
                 dplyr::select(angio, PC1, PC2) %>% filter(complete.cases(.), abs(PC1)< 1, abs(PC2)<1) %>% rename(Genus = angio), aes(x = PC1, y = PC2))+
    geom_point(aes(color = Genus, fill = Genus, shape = Genus, alpha = Genus), size = 0.35)+
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
      scale_alpha_manual(values = alpha_scale)
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
      scale_alpha_manual(values = alpha_scale)
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

#----------------------------- PCA axes 2-3 ----------------
ggpca_weighted_2.3 <- function(tr_use, trait_vars, num_show, alpha = 0.1, flip_coord = FALSE, arrow_head = 6, 
                               flip_x = 1, flip_y = 1, use_vars = NULL, label_size = 5, show_labs = TRUE, alpha_const = 1,
                               alpha_all = 0, line_wd_big = 1.5, all_alpha = 0.2, angio_list = NULL, flip_z = 1, flip_zz = 1,
                               show_legend = TRUE, title = ""){
  
  
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
  pca_use <- wpca(x = x_pca, w = w_pca, center=TRUE, scale=TRUE)
  
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
  
  # make sure wearen't displaying more traits than we have
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
  
  # if we haven't specified which variables to plot (Fig 2a), calculate the ones with the highest loadings
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
  
  # create the base plot                          NB : here change manually the axis of points coord. maybe discuss w dan to streamline?
  g1 <- ggplot(data = plot_dt %>% 
                 dplyr::select(angio, PC2, PC3) %>% filter(complete.cases(.), abs(PC2)< 1, abs(PC3)<1) %>% rename(Genus = angio), aes(x = PC2, y = PC3))+
    geom_point(aes(color = Genus, fill = Genus, shape = Genus, alpha = Genus), size = 0.35)+
    theme_bw()
  
  # create the transparancy value based on number of observations (note that Genus is set to angio vs. gymno in the main figure)
  alv <- plot_dt %>% dplyr::select(angio, PC2, PC3) %>% filter(complete.cases(.), abs(PC2)< 1, abs(PC3)<1) %>% rename(Genus = angio) %>% 
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
      scale_alpha_manual(values = alpha_scale)
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
      scale_alpha_manual(values = alpha_scale)
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
    mutate(#FDis = rep(NA),                    
           #FR = rep(NA),
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
  
  
  # # asses functional dispersion  
  # set.seed(1312)
  # for (i in 1:length(hype@HVList)){
  #   ifelse(hype@HVList[[i]]@Dimensionality != 1,  kernel.dispersion(hype@HVList[[i]], func = "divergence", frac =0.1 ), 0 )} -> FDis_results[i,2]
  # 
  # # asses functional richness  
  # set.seed(1312)
  # for (i in 1:length(hype@HVList)){
  #   kernel.alpha(hype@HVList[[i]])} -> FDis_results[i,3]
  # 
  # #assess functional eveness
  # set.seed(1312)
  # for (i in 1:length(hype@HVList)){
  #   ifelse(hype@HVList[[i]]@Dimensionality != 1, kernel.evenness(hype@HVList[[i]]),0)} -> FDis_results[i,4]
  
  #get community centroid value
  set.seed(1312)
  for (i in 1:length(hype@HVList)){
    get_centroid(hype@HVList[[i]]) %>% as.list() -> FDis_results[i,5:7]
  }
  #5) dplyr::select and print
  return(FDis_results) 
  
  
}

#---------- Function to make maps (Fig.3) --------------------------------------
Map_FD_metric <- function(Data, variable_to_map, legend.name){

# create sf object from point
Data %>% group_by(country) %>% filter(census.n == max(census.n)) %>% ungroup() %>% 
   dplyr::select(variable_to_map, LON, LAT) %>% 
   st_as_sf( coords = c("LON", "LAT")) %>% rename(x = variable_to_map)-> points_sf 

 # create a grid covering all points, cell size 100000 units,
# add unique cell_id
grid_sf <- st_make_grid(points_sf, cellsize = 0.2) |>  
  st_sf(geometry = _) %>% 
  mutate(cell_id = row_number(), .before = 1)
print(grid_sf, n = 3)


# spatial join to metch cell_id-s to points,
shares_df <- st_join(points_sf, grid_sf) %>% 
  st_drop_geometry() %>% 
  # number of points per cell (mutate, row count remains the same)
  group_by(cell_id) %>% 
  mutate(x = mean(x)) %>%   
  ungroup()

shares_df

# join grid to shares for ploting
gird_shares_sf <- right_join(grid_sf, shares_df)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-
#get EU map
# Get country contours 
my_world_map <- map_data('world')
# To check the country names (in case necessary)
as.character(unique(unlist(my_world_map$region))) 
# Specify the countries you want in a character vector 
my_world_map1 <- my_world_map[my_world_map$region %in% c(  "Portugal", "Spain", "France", "Switzerland", "Germany",
                                                           "Austria", "Belgium", "UK", "Netherlands",
                                                           "Denmark", "Poland", "Italy", 
                                                           "Croatia", "Slovenia", "Hungary", "Slovakia",
                                                           "Czech Republic", "Finland", "Sweden", "Norway", "Greece", 
                                                           "Serbia", "Bosnia and Herzegovina", "Ireland",
                                                           "Estonia", "Latvia", "Lithuania", "Albania", "Romania" ,  
                                                           "Bulgaria", "Luxembourg", "North Macedonia", "Kosovo", "Montenegro"),]

#plot fdis
gird_shares_sf %>% 
  ggplot() +
  geom_sf(aes(fill = x), colour = NA)+ 
  geom_sf(data = grid_sf, fill = NA, color = NA) +
  scale_fill_distiller(palette = "Spectral", direction = 2, na.value = NA, name = legend.name)+
  geom_polygon(data = my_world_map1, aes(x = long, y = lat, group=group) ,alpha = 0.75, linewidth=0.75,  colour = 'black', fill=  NA) +
  ylim(36, 71)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+theme(panel.grid.major = element_blank())
}






