#############
## Custom functions - Thornback ray
#############

##### Percent change 
percent_log_se_change <- function(data, model_var = "model", year_var = "year", se_var = "se") {
  
  # Mean change in annual log-SE for each model
  base_dat <- data |>
    filter(model == "base") |> 
    select(year, se) |>
    mutate(
      log_se = log(se),
      log_se_lag = lag(log_se),
      pct_change_log_se = 100 * (log_se - log_se_lag) / abs(log_se_lag)
    )
  base_res <- mean(base_dat$pct_change_log_se, na.rm = TRUE)
  
  cov_dat <- data |>
    filter(model != "base") |> 
    select(year, se) |>
    mutate(
      log_se = log(se),
      log_se_lag = lag(log_se),
      pct_change_log_se = 100 * (log_se - log_se_lag) / abs(log_se_lag)
    )
  cov_res <- mean(cov_dat$pct_change_log_se, na.rm = TRUE)
  
  # Mean log-SE for each model
  base_res2 <- mean(base_dat$log_se, na.rm = TRUE)
  cov_res2 <- mean(cov_dat$log_se, na.rm = TRUE)
  res2 <- c(base_res2, cov_res2)
  
  # Make an output object
  output <- paste0(
    "===========================\n",
    "       Report Summary      \n",
    "===========================\n",
    "Mean percent change in log SE of index: ", "\n",
    "Base model: ", base_res, "\n",
    "Covariate model: ", cov_res, "\n",
    "Status: ", 
    ifelse(cov_res < base_res, 
           "Value of covariate model is larger than base model",
           "Value of covariate model is smaller than base model"), "\n",
    "\n",
    "Mean log SE of index: ", "\n",
    "Base model: ", base_res2, "\n",
    "Covariate model: ", cov_res2, "\n",
    "===========================\n"
  )
  cat(output)
  return(output)
  
}

##### Make deviance report
deviance_report <- function(full_cov_model = "full_cov_model", 
                            null_model = "null_model", 
                            base_model = "base_model", 
                            reduced_models) {
  
  dev_null <- deviance(null_model)
  
  # Get deviance explained from base model
  dev_base <- deviance(base_model)
  pde_base <- 1 - dev_base / dev_null
  
  # Get deviance explained from covariate model
  dev_full <- deviance(full_cov_model)
  pde_full <- 1 - dev_full / dev_null
  
  cat("Base model deviance:", dev_base, "\n")
  cat("Full covariate model deviance:", dev_full, "\n")
  cat("Null model deviance:", dev_null, "\n")
  cat("Proportion deviance explained (base):", round(pde_base, 4), "\n\n")
  cat("Proportion deviance explained (full covariate):", round(pde_full, 4), "\n\n")
  
  # Loop over reduced models
  for (i in seq_along(reduced_models)) {
    model <- reduced_models[[i]]
    dev_reduced <- deviance(model)
    pde_reduced <- 1 - dev_reduced / dev_null
    drop_in_pde <- pde_full - pde_reduced
    
    name <- deparse(substitute(model))
    cat("Model", i, "\n")
    cat(" Model formula:", paste(as.character(model$formula[[1]][3])), "\n")
    cat("  Reduced deviance:", dev_reduced, "\n")
    cat("  Proportion deviance explained:", round(pde_reduced, 4), "\n")
    cat("  Drop in PDE from full model:", round(drop_in_pde, 4), "\n\n")
  }
  
}

##### Crossvalidation

# Updated functions for repeated crossvalidation
# Create partitions
holdout <- function(model, training_cutoff, seeds) {
  require(dplyr)
  data <- model$data
  n <- nrow(model$data)
  #test_cutoff <- 1 - training_cutoff # determine proportional cut-off
  training_percent <- (n * training_cutoff) %>% floor
  
  output <- list()
  
  for(s in seeds) {
    set.seed(s)
    train_sample <- sample(1:n, training_percent) # randomly pick rows for training
    test_sample <- setdiff(1:n, train_sample) # get the remaining n% of the rows
    
    output[[as.character(s)]] <- list(
      training = model$data[train_sample, ], 
      test = model$data[test_sample, ]
    )
  }
  
  return(output)
}

# Run repeated 70:30 crossvalidation
rep_cv <- function(test_training = test_training, # supply either base or cov
                   formula = c("base", "full", "full2", "depth_substrate")) {
  
  # Run models using training data
  mesh <- make_mesh(test_training$training, c("X", "Y"), # original mesh settings
                    fmesher_func = fmesher::fm_mesh_2d_inla,
                    cutoff = 10, max.edge = c(75, 100), offset = 10) 
  
  # Define formula
  if (formula == "base") {
    f <- count ~ fgear + fyear
  } else if (formula == "depth_substrate") {
    f <- count ~ fgear + fyear + s(middepth_scaled, m = 1) + substrate
  } else if (formula == "full2") {
    f <- count ~ fgear + fyear + s(bottomT_scaled, m = 1) + 
      s(middepth_scaled, m = 1, by = fgear) + substrate
  } else {
    f <- count ~ fgear + fyear + s(bottomT_scaled, m = 1) + 
      s(middepth_scaled, m = 1) + substrate
  }
  
  cv <- sdmTMB(
    data = test_training$training,
    list(f,
         f),
    mesh = mesh,
    family = nbinom2(),
    offset = test_training$training$log_sweptareakmsqadj,
    spatial = "on",
    time = "year", 
    spatiotemporal = "IID",
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  # Predict onto test data
  cv_preds <- predict(cv, newdata = test_training$test, type = "response")
  
  # RMSE
  rmse <- sqrt(mean((cv_preds$count - cv_preds$est)^2, na.rm = TRUE))
  
  # MAE
  mae <- mean(abs(cv_preds$count - cv_preds$est), na.rm = TRUE)
  
  tibble::tibble(
    RMSE = rmse,
    MAE = mae
  )
}

##### Spatial autocorrelation functions
library(dbscan)
library(spdep)
library(kableExtra)
library(deldir)
library(sf)

# Function to iterate over years to find the lowest k at full connectivity,
# and make a spatial weights list for global Moran's
make_kn_dist_obj <- function(coords) {
  
  kn_dist_obj <- list()
  years <- unique(coords$year)
  max_k = 10
  
  for (y in years) {
    y_coords <- subset(coords, year == y)
    min_k <- list()
    
    for (k in 1:max_k) {
      crds <- cbind(y_coords$lon, y_coords$lat)
      knn <- knearneigh(crds, k = k) 
      nb <- knn2nb(knn)
      
      comps <- n.comp.nb(nb)
      
      if (comps$nc == 1) {
        
        # make and save graph-based k-nearest neighbour connections 
        
        plot(nb, crds, main = paste("Connected KNN graph (k =", k, ") in year", y))
        
        # specify spatial weights for neighbours
        
        min_k[[length(min_k) + 1]] <- list(k = k, nb = nb, nblistw = nb2listw(nb), 
                                           coordsdat = y_coords)
        #min_k <- list(k = k, nb = nb, nblistw = nb2listw(nb))
        #return(nb_plot)
        
        break
        
      }
    }
    
    kn_dist_obj[[as.character(paste0("y", y))]] <- min_k
    
  }
  
  return(kn_dist_obj)
  
}

# Plotting out global Moran's results
globm_plot_fun <- function(moran_data, model_comp = c(1, 2)) {
  
  if (1 %in% model_comp) {
    xname <- "resids"
    p <- ggplot(moran_data, aes(x=x, y=wx)) + 
      geom_point(shape=1, size = 0.8) + 
      geom_smooth(formula=y ~ x, method="lm", linewidth = 0.2) + 
      geom_hline(yintercept=mean(moran_data$wx), lty=2) + 
      geom_vline(xintercept=mean(moran_data$x), lty=2) + 
      #theme_minimal() + 
      cowplot::theme_cowplot() +
      theme(plot.margin = margin(-1, -1, -1, -1),
            aspect.ratio = 1,
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8)) +
      geom_point(data=moran_data[moran_data$is_inf,], aes(x=x, y=wx), 
                 shape=9, size = 0.8) +
      geom_text(data=moran_data[moran_data$is_inf,], 
                aes(x=x, y=wx, label=labels), vjust=1.5 , size = 3) +
      xlab(xname) + ylab(paste0("Spatially lagged ", xname))
    
    return(p)
  }
  
  if (2 %in% model_comp) {
    xname <- "resids2"
    p <- ggplot(moran_data, aes(x=x, y=wx)) + 
      geom_point(shape=1, size = 0.8) + 
      geom_smooth(formula=y ~ x, method="lm", linewidth = 0.2) + 
      geom_hline(yintercept=mean(moran_data$wx), lty=2) + 
      geom_vline(xintercept=mean(moran_data$x), lty=2) + 
      #theme_minimal() + 
      cowplot::theme_cowplot() +
      theme(plot.margin = margin(-1, -1, -1, -1),
            aspect.ratio = 1,
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8)) +
      geom_point(data=moran_data[moran_data$is_inf,], aes(x=x, y=wx), 
                 shape=9, size = 0.8) +
      geom_text(data=moran_data[moran_data$is_inf,], 
                aes(x=x, y=wx, label=labels), vjust=1.5 , size = 3) +
      xlab(xname) + ylab(paste0("Spatially lagged ", xname))
    
    return(p)
  }
  
}

# Calculate local Moran's I for each data point
# function to calculate local Moran's based on conditional permutations
make_local_morans_obj <- function(.x) {
  # Calculate conditional Local Moran's and adjust p-values
  locm_p <- localmoran_perm(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, 
                            nsim = 999, iseed = 123)
  .x[[1]]$coordsdat <- cbind(.x[[1]]$coordsdat, locm_p[, c("Ii", "E.Ii", 
                                                           "Var.Ii", "Z.Ii",
                                                           "Pr(z != E(Ii))",
                                                           "Pr(z != E(Ii)) Sim")
  ])
  .x[[1]]$coordsdat$locm_p_pv <- p.adjust(locm_p[, "Pr(z != E(Ii)) Sim"],
                                          "fdr")
  .x[[1]]$coordsdat$locm_significance <- as.factor(ifelse(.x[[1]]$coordsdat$locm_p_pv <= 0.01, 
                                                          print("sig"), print("insig")))
  
  # Get quadrant for hotspots
  .x[[1]]$coordsdat$quad_mean <- attr(locm_p, "quadr")$mean 
  locm_dat <- .x[[1]]$coordsdat
  return(locm_dat)
}

# function to adjust p-values and check against alpha threshold 

#f <- function(x) sum(x < 0.005) # choose 0.005 as "interesting" alpha threshold
check_psig_locm <- function(x) {
  pv <- x |>
    subset(select = "Pr(z != E(Ii)) Sim", drop = TRUE) 
  
  pva <- function(pv) cbind("none" = pv, 
                          "FDR" = p.adjust(pv, "fdr"), "BY" = p.adjust(pv, "BY"),
                          "Bonferroni" = p.adjust(pv, "bonferroni"))
  pvsp <- pva(pv)
  
  f <- function(x) sum(x < 0.01)
  
  check <- apply(pvsp, 2, f)
  return(check)
}

### Calculate local Getis-Ord 
make_getis_ord_obj <- function(.x) {
  # Calculate conditional Getis-Ord and adjust p-values
  locg_p <- localG_perm(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, 
                        nsim = 999, iseed = 123)
  locg_p_df <- attr(locg_p, "internals")
  .x[[1]]$coordsdat <- cbind(.x[[1]]$coordsdat, locg_p_df[, c("Gi", "E.Gi", 
                                                              "Var.Gi", "StdDev.Gi",
                                                              "Pr(z != E(Gi))",
                                                              "Pr(z != E(Gi)) Sim")
  ])
  
  .x[[1]]$coordsdat$locg_p_pv <- p.adjust(locg_p_df[, "Pr(z != E(Gi)) Sim"],
                                          "fdr")
  .x[[1]]$coordsdat$locg_significance <- ifelse(.x[[1]]$coordsdat$locg_p_pv <= 0.01, 
                                                print("sig"), print("insig"))
  
  # Get hotspot clusters
  locg_p_clust <- attr(locg_p, "cluster")
  .x[[1]]$coordsdat <- cbind(.x[[1]]$coordsdat, locg_p_clust)
  
  locg_dat <- .x[[1]]$coordsdat
  return(locg_dat)
}

check_psig_locg <- function(x) {
  pv <- x |>
    subset(select = "Pr(z != E(Gi)) Sim", drop = TRUE)
  
  pva <- function(pv) cbind("none" = pv, 
                            "FDR" = p.adjust(pv, "fdr"), "BY" = p.adjust(pv, "BY"),
                            "Bonferroni" = p.adjust(pv, "bonferroni"))
  pvsp <- pva(pv)
  
  f <- function(x) sum(x < 0.01)
  
  check <- apply(pvsp, 2, f)
  return(check)
}

# Plotting out local Moran's and Getis Ord results
locm_plot_fun <- function(full_obj, ire, coords, hull, Ii) {
  ggplot() +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0  # Center of your scale
    ) +
    #scale_fill_viridis_c() + 
    geom_sf(data = full_obj, aes(fill = .data[[Ii]])) +
    geom_sf(data = ire, colour = "grey45", fill = "grey45") +
    geom_sf(data = coords, colour = "black", size = 0.001) +
    #geom_sf(data = coords, aes(shape = .data[[locm_significance]]), 
    #        colour = "black", size = 0.5) +
    #scale_shape_manual(values = c(1, 19)) +
    geom_sf(data = hull, colour = "black", fill = NA) +
    theme(
      legend.box.spacing = unit(1, "mm"),
      legend.key.width = unit(0.5, "mm"),
      legend.key.height = unit(1.5, "mm"),
      legend.spacing = unit(1, "mm"),
      legend.title = element_text(size = 4),
      legend.text = element_text(size = 3),
      legend.margin = margin(-1, 1, -1, -1),
      plot.margin = margin(-1, -1, -1, -1),
      aspect.ratio = 1,
      axis.text = element_blank(), 
      axis.ticks = element_blank()
    ) #+
  #labs(shape = "Significance", fill = "Local Moran's I")
  #labs(fill = "Local Moran's I")
}

locg_plot_fun <- function(full_obj, ire, coords, hull, Gi) {
  ggplot() +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0  # Center of your scale
    ) +
    #scale_fill_viridis_c() + 
    geom_sf(data = full_obj, aes(fill = .data[[Gi]])) +
    geom_sf(data = ire, colour = "grey45", fill = "grey45") +
    geom_sf(data = coords, colour = "black", size = 0.001) +
    #geom_sf(data = coords, aes(shape = .data[[locg_significance]]), 
    #        colour = "black", size = 0.5) +
    #scale_shape_manual(values = c(1, 19)) +
    geom_sf(data = hull, colour = "black", fill = NA) +
    theme(
      legend.box.spacing = unit(1, "mm"),
      legend.key.width = unit(0.5, "mm"),
      legend.key.height = unit(1.5, "mm"),
      legend.spacing = unit(1, "mm"),
      legend.title = element_text(size = 4),
      legend.text = element_text(size = 3),
      legend.margin = margin(-1, 1, -1, -1),
      plot.margin = margin(-1, -1, -1, -1),
      aspect.ratio = 1,
      axis.text = element_blank(), 
      axis.ticks = element_blank()
      #  legend.position = c(.09, .85),
      #  legend.background = element_rect(linetype="solid", colour = "black")
      #legend.justification = c("left", "top")#,
      #legend.box.just = "left",
      #legend.margin = margin(2, 2, 2, 2)
    ) #+
  #labs(fill = "Getis-Ord G")
}

locg_hotspot_fun <- function(full_obj, ire, coords, hull, locg_p_clust) {
  ggplot() +
    scale_fill_manual(values = c("blue", "red")) +
    #scale_fill_viridis_c() + 
    geom_sf(data = full_obj, aes(fill = .data[[locg_p_clust]])) +
    geom_sf(data = ire, colour = "grey45", fill = "grey45") +
    geom_sf(data = coords, colour = "black", size = 0.001) +
    #geom_sf(data = coords, aes(shape = .data[[locg_significance]]), 
    #        colour = "black", size = 0.5) +
    #scale_shape_manual(values = c(1, 19)) +
    geom_sf(data = hull, colour = "black", fill = NA) +
    theme(
      legend.box.spacing = unit(1, "mm"),
      legend.key.width = unit(0.5, "mm"),
      legend.key.height = unit(1.5, "mm"),
      legend.spacing = unit(1, "mm"),
      legend.title = element_text(size = 4),
      legend.text = element_text(size = 3),
      legend.margin = margin(-1, 1, -1, -1),
      plot.margin = margin(-1, -1, -1, -1),
      aspect.ratio = 1,
      axis.text = element_blank(), 
      axis.ticks = element_blank())
}

# function to make convex hull and voronoi tesselation
make_hull_and_vortess <- function(dat, shape) {
  years <- unique(dat$year)
  hull_vortess <- list() 
  
  for (y in years) {
    subdat <- subset(dat, year == y)
    
    # Build convex hull 
    y_coords <- subdat[, c("lon", "lat")] 
    y_hull <- chull(y_coords)
    
    # Extract coords defining hull
    y_hull_coords <- y_coords[y_hull, ]
    # add first point coordinates to close polygon
    y_hull_coords <- rbind(y_hull_coords, y_hull_coords[1, ])  
    
    # Make Voronoi tesselation of year subset coordinates
    y_vortess <- with(y_coords, deldir(lon, lat, 
                                       rw = c(min(lon), max(lon), min(lat), max(lat))))
    y_tiles <- tile.list(y_vortess)
    # make hull an sf object
    y_tiles_sf <- st_as_sf.deldir(y_vortess, extract = c("tiles", "triangles"))
    
    # Make a polygon out of 
    y_hull_poly <- st_sfc(st_polygon(list(as.matrix(y_hull_coords[, c("lon", "lat")]))))
    y_hull_poly_sf <- st_sf(geometry = y_hull_poly, crs = 4326)
    y_hull_poly_utm <- st_transform(y_hull_poly_sf, 32629)
    
    # Crop voronoi tesselation to hull boundaries
    y_intersect <- st_intersection(y_tiles_sf, y_hull_poly_utm)
    
    # Transform Irish landmass for plotting
    ire_utm <- st_transform(ire, st_crs(y_hull_poly_utm))
    
    # Make coords to sf object for plotting
    subdat_sf <- st_as_sf(subdat, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
    subdat_utm <- st_transform(subdat_sf, 32629)
    
    # Join test values for each point to voronoi tesselation
    vt_hull_locm <- st_join(y_intersect, subdat_utm)
    
    hull_vortess[[as.character(paste0("y", y))]] <- list(
      hull = y_hull_poly_utm, #vortess = y_intersect, 
      coords = subdat_utm, full_obj = vt_hull_locm,  
      ire = ire_utm)
  } 
  
  return(hull_vortess)
  
}

st_as_sf.deldir <- function(dd, extract = c("tiles", "triangles")) {
  # from https://gist.github.com/joeroe/43bb6efa233cb5c8994b6a6d41a4911a
  extract <- match.arg(extract)
  
  if (extract == "tiles") {
    ddlist <- deldir::tile.list(dd)
  }
  else if (extract == "triangles") {
    ddlist <- deldir::triang.list(dd)
  }
  
  ddlist %>%
    purrr::map(~{cbind(x = .$x, y = .$y)} %>%
                 rbind(.[1,]) %>%
                 list() %>%
                 sf::st_polygon()) %>%
    sf::st_sfc() %>%
    sf::st_sf() %>%
    st_set_crs(4326) %>%
    st_transform(32629) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    return()
}

# Make the Voronoi tesselation and hull
#st_as_sf.deldir2 <- function(dd, extract = c("tiles", "triangles")) {
# from https://gist.github.com/joeroe/43bb6efa233cb5c8994b6a6d41a4911a
#  extract <- match.arg(extract)

#  if (extract == "tiles") {
#    ddlist <- deldir::tile.list(dd)
#  }
#  else if (extract == "triangles") {
#    ddlist <- deldir::triang.list(dd)
#  }

#  ddlist %>%
#    purrr::map(~{cbind(x = .$x, y = .$y)} %>%
#                 rbind(.[1,]) %>%
#                 list() %>%
#                 sf::st_polygon()) %>%
#    sf::st_sfc() %>%
#    #sf::st_sf() %>%
#    st_set_crs(4326) %>%
#    st_transform(32629) -> dd_sfc
#  
#  dd_pts <- as.vector(tidytable::map_int(tiles, ~ .$ptNum[1]))
#  dd_sfc$pt <- dd_pts
#  return(dd_sfc)
#}

##### Plotting out effects predictions and their uncertainty 
plot_map <- function(dat, column) {
  ggplot(dat, aes(round(X, 0), round(Y, 0), fill = {{ column }})) +
    geom_tile(width = 5, height = 5) +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map2 <- function(dat, column) {
  ggplot(dat, aes(round(X, 0), round(Y, 0), fill = {{ column }})) +
    geom_tile(width = 5, height = 5) +
    #facet_wrap(~year) +
    coord_fixed()
}


