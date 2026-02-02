#############
## Custom functions - Horse mackeral
#############

##### Percent change 
percent_log_se_change <- function(data, model_var = "model", 
                                  year_var = "year", se_var = "se") {
  
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
  res <- c(base_res, cov_res)
  
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

##### External consistency plots 
# flip direction of diagonal panels in plots
flip_panel <- function(x, y, ...) {
  # Perform your custom plotting logic here
  # For example, reversing the order of points or axes
  plot(rev(x), rev(y), ...) # Simple example of flipping
}



##### Format marginal effects plots
format_mareff_plots <- function(mareff1 = "mareff1", mareff2 = "mareff2") {
  
  mareff1$theme <- theme_bw()
  mareff1$labels$y <- "Encounter probability"
  mareff1$theme$axis.title <- element_text(size = 16, face = "bold")
  mareff1$theme$axis.text <- element_text(size = 15)
  mareff2$theme <- theme_bw()
  mareff2$labels$y <- "Biomass"
  mareff2$theme$axis.title <- element_text(size = 16, face = "bold")
  mareff2$theme$axis.text <- element_text(size = 15)
  
  res <- mareff1 / mareff2 
  
  return(res)
  
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

##### Function to create test/training data sets cross validation tests

holdout <- function(model, training_cutoff) {
  require(dplyr)
  data <- model$data
  n <- nrow(model$data)
  test_cutoff <- 1 - training_cutoff

  training_percent <- (n * training_cutoff) %>% floor
  train_sample <- sample(1:n, training_percent) # randomly pick rows for training
  test_sample <- setdiff(1:n, train_sample) # get the remaining n% of the rows
  
  train <- model$data[train_sample, ] 
  test <- model$data[test_sample, ] 
  
  output <- list(training = train, test = test)
  return(output)
}


##### Function to report outcomes of cross validation tests
#crossval_report <- function(base_cv = "base_cv", cov_cv = "cov_cv") {
#  cat("BASE MODEL", "\n")
#  cat(length(base_cv$fold_loglik), "-fold crossvalidation)", "\n")
#  cat(" Model validation: predictive accuracy with fitted model", "\n")
#  cat("  RMSE between original abundance and CV predictions across entire dataset:", 
#      sqrt(mean((base_cv$data$numkm - base_cv$data$cv_predicted)^2)), "\n")
#  cat("  MAE between original abundance and CV predictions across entire dataset:",
#      mean(abs(base_cv$data$numkm - base_cv$data$cv_predicted)), "\n")
#  cat("  All models converged?", base_cv$converged, "\n")
#  cat("\n")
#  cat("BEST FIT COVARIATE MODEL", "\n")
#  cat(length(cov_cv$fold_loglik), "-fold crossvalidation)", "\n")
#  cat(" Model validation: predictive accuracy with fitted model", "\n")
#  cat("  RMSE between original abundance and CV predictions across entire dataset:",
#      sqrt(mean((cov_cv$data$numkm - cov_cv$data$cv_predicted)^2)), "\n")
#  cat("  MAE between original abundance and CV predictions across entire dataset:", 
#      mean(abs(cov_cv$data$numkm - cov_cv$data$cv_predicted)), "\n")
#  cat("  All models converged?", cov_cv$converged, "\n")
#  cat("\n")
#  cat("BASE MODEL VS COVARIATE MODEL", "\n")
#  cat("Base model:", "\n")
#  cat(" Log-liks of each fold:", "\n")
#  cat(base_cv$fold_loglik, "\n")
#  cat(" Sum of all model loglikelihoods (on same mesh):", 
#      base_cv$sum_loglik, "\n")
#  cat("Covariate model:", "\n")
#  cat(" Log-liks of each fold:", "\n")
#  cat(cov_cv$fold_loglik, "\n")
#  cat(" Sum of all model loglikelihoods (on same mesh):", 
#      cov_cv$sum_loglik, "\n")
#  cat("Sum of covariate model logliks > Sum of base model logliks?", 
#      cov_cv$sum_loglik > base_cv$sum_loglik, "\n")
  
#}


crossval_report <- function(base_preds, cov_preds) {
  cat("BASE MODEL", "\n")
  cat("70:30 Holdout crossvalidation", "\n")
  cat(" Model validation: predictive accuracy with fitted model", "\n")
  cat("  RMSE between original abundance and CV predictions across entire dataset:", 
      sqrt(mean((base_preds$numkm - base_preds$est)^2, na.rm = TRUE)), "\n")
  cat("  MAE between original abundance and CV predictions across entire dataset:",
      mean(abs(base_preds$numkm - base_preds$est), na.rm = TRUE), "\n")
  cat("\n")
  cat("BEST FIT COVARIATE MODEL", "\n")
  cat("70:30 Holdout crossvalidation", "\n")
  cat(" Model validation: predictive accuracy with fitted model", "\n")
  cat("  RMSE between original abundance and CV predictions across entire dataset:",
      sqrt(mean((cov_preds$numkm - cov_preds$est)^2, na.rm = TRUE)), "\n")
  cat("  MAE between original abundance and CV predictions across entire dataset:", 
      mean(abs(cov_preds$numkm - cov_preds$est), na.rm = TRUE), "\n")
  cat("\n")
}





model_auc <- function(cv) {
  cv_base <- cv$basemodel_cv
  cv_cov <- cv$covmodel_cv
  
  roc_base_m1 <- pROC::roc(cv_base$data$present, plogis(cv_base$data$cv_predicted))
  roc_base_m2 <- pROC::roc(cv_base$data$numkm, log(cv_base$data$cv_predicted))
  roc_cov_m1 <- pROC::roc(cv_cov$data$present, plogis(cv_cov$data$cv_predicted))
  roc_cov_m2 <- pROC::roc(cv_cov$data$numkm, log(cv_cov$data$cv_predicted))
  
  cat("BASE MODEL", "\n")
  cat(length(cv_base$fold_loglik), "-fold crossvalidation", "\n")
  cat(" Binomial predictor", "\n")
  cat(" Model performance: predictive accuracy with fitted model", "\n")
  cat("  AUC between original abundance and CV predictions across entire dataset:", 
      pROC::auc(roc_base_m1), "\n")
  cat(" Lognormal predictor", "\n")
  cat(" Model performance: predictive accuracy with fitted model", "\n")
  cat("  AUC between original abundance and CV predictions across entire dataset:", 
      pROC::auc(roc_base_m2), "\n")
  cat("  All models converged?", cv_base$converged, "\n")
  cat("\n")
  cat("BEST FIT COVARIATE MODEL", "\n")
  cat(length(cv_cov$fold_loglik), "-fold crossvalidation", "\n")
  cat(" Binomial predictor", "\n")
  cat(" Model performance: predictive accuracy with fitted model", "\n")
  cat("  AUC between original abundance and CV predictions across entire dataset:", 
      pROC::auc(roc_cov_m1), "\n")
  cat(" Lognormal predictor", "\n")
  cat(" Model performance: predictive accuracy with fitted model", "\n")
  cat("  AUC between original abundance and CV predictions across entire dataset:", 
      pROC::auc(roc_cov_m2), "\n")
  cat("  All models converged?", cv_cov$converged, "\n")
  cat("\n")
  
}

##### DHARMa functions
make_dharma <- function(model = "model", model_comp = c(1, 2)) {
  out <- list()
  
  # More efficient - Predict once!
  p <- predict(model, type = "response")
  
  if (1 %in% model_comp) {
    s1 <- simulate(model, nsim = 1000, model = 1, type = "mle-mvn")
    res1 <- DHARMa::createDHARMa(
      simulatedResponse = s1,
      observedResponse = model$data$present,
      fittedPredictedResponse = p$est1
    )
    out[["m1"]] <- res1
  }
  
  if (2 %in% model_comp) {
    pos <- model$data$numkm > 0
    s2 <- simulate(model, nsim = 1000, model = 2, type = "mle-mvn")
    res2 <- DHARMa::createDHARMa(
      simulatedResponse = s2[pos, ],
      observedResponse = model$data$numkm[pos],
      fittedPredictedResponse = p$est2[pos]
    )
    out[["m2"]] <- res2
  }
  
  # Return single item directly if only one component was selected
  if (length(out) == 1) return(out[[1]])
  
  return(out)
  
}

##### Spatial autocorrelation functions
library(dbscan)
library(spdep)
library(kableExtra)

# Function to iterate over years to find the lowest k at full connectivity,
# and make a spatial weights list for global Moran's
make_kn_dist_obj_binom <- function(coords) {
  kn_dist_obj <- list()
  years <- unique(coords$year)
  max_k = 40
  # Note: spdep also allows distance-based clustering - I just went with k-based method 
  
  # helper function to find the smallest k with connected graph
  get_min_k_connected <- function(y_coords, label, y) {
    min_k <- list()
    crds <- cbind(y_coords$lon, y_coords$lat)
    
    for (k in 1:max_k) {
      crds <- cbind(y_coords$lon, y_coords$lat)
      knn <- knearneigh(crds, k = k) 
      nb <- knn2nb(knn)
      
      comps <- n.comp.nb(nb)
      
      if (comps$nc == 1) {
        # make and save graph-based k-nearest neighbour connections 
        plot(nb, crds, 
             main = paste("Connected KNN graph (", label, ", k =", k, ") in year", y))
        # specify spatial weights for neighbours
        min_k[[length(min_k) + 1]] <- list(k = k, nb = nb, nblistw = nb2listw(nb), 
                                           coordsdat = y_coords)
        break
        
      }
    }
    return(min_k)
  } 
  for (y in years) {
    y_coords <- subset(coords, year == y)
    binom <- y_coords[, c("year", "lon", "lat", "resids")]
    #pos <- y_coords[, c("year", "lon", "lat", "resids2")]
    
    min_k_binom <- get_min_k_connected(binom, "binom", y)
    #min_k_pos <- get_min_k_connected(pos, "pos", y)
    
    kn_dist_obj[[as.character(paste0("y", y))]] <- min_k_binom
    #kn_dist_obj[[as.character(paste0("y", y))]] <- list(min_k_binom, 
    #                                                    min_k_pos)
  }
  return(kn_dist_obj)
}

make_kn_dist_obj_pos <- function(coords) {
  kn_dist_obj <- list()
  years <- unique(coords$year)
  max_k = 40
  
  # helper function to find the smallest k with connected graph
  get_min_k_connected <- function(y_coords, label, y) {
    min_k <- list()
    crds <- cbind(y_coords$lon, y_coords$lat)
    
    for (k in 1:max_k) {
      crds <- cbind(y_coords$lon, y_coords$lat)
      knn <- knearneigh(crds, k = k) 
      nb <- knn2nb(knn)
      
      comps <- n.comp.nb(nb)
      
      if (comps$nc == 1) {
        # make and save graph-based k-nearest neighbour connections 
        plot(nb, crds, 
             main = paste("Connected KNN graph (", label, ", k =", k, ") in year", y))
        
        # specify spatial weights for neighbours
        min_k[[length(min_k) + 1]] <- list(k = k, nb = nb, 
                                           # zero.policy accounts for no catch entries
                                           nblistw = nb2listw(nb),  
                                           coordsdat = y_coords)
        break
      }
    }
    return(min_k)
  } 
  
  for (y in years) {
    y_coords <- subset(coords, year == y)
    #binom <- y_coords[, c("year", "lon", "lat", "resids")]
    pos <- y_coords[, c("year", "lon", "lat", "resids2")]
    pos <- pos[! is.na(pos$resids2), ] # remove NAs
    
    #min_k_binom <- get_min_k_connected(binom, "binom", y)
    min_k_pos <- get_min_k_connected(pos, "pos", y)
    
    kn_dist_obj[[as.character(paste0("y", y))]] <- min_k_pos
    #kn_dist_obj[[as.character(paste0("y", y))]] <- list(min_k_binom, 
    #                                                    min_k_pos)
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
    y_vortess <- with(y_coords, deldir::deldir(lon, lat, 
                                               rw = c(min(lon), max(lon), min(lat), max(lat))))
    y_tiles <- deldir::tile.list(y_vortess)
    # make hull an sf object
    y_tiles_sf <- st_as_sf.deldir(y_vortess, extract = c("tiles", "triangles"))
    
    # Make a polygon out of 
    y_hull_poly <- st_sfc(st_polygon(list(as.matrix(y_hull_coords[, c("lon", "lat")]))))
    y_hull_poly_sf <- st_sf(geometry = y_hull_poly, crs = 4326)
    y_hull_poly_utm <- st_transform(y_hull_poly_sf, 32629)
    
    # Crop voronoi tesselation to hull boundaries
    y_intersect <- st_intersection(y_tiles_sf, y_hull_poly_utm)
    
    # Transform Irish landmass for plotting
    shape_utm <- st_transform(shape, st_crs(y_hull_poly_utm))
    
    # Make coords to sf object for plotting
    subdat_sf <- st_as_sf(subdat, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
    subdat_utm <- st_transform(subdat_sf, 32629)
    
    # Join test values for each point to voronoi tesselation
    vt_hull_locm <- st_join(y_intersect, subdat_utm)
    
    hull_vortess[[as.character(paste0("y", y))]] <- list(
      hull = y_hull_poly_utm, #vortess = y_intersect, 
      coords = subdat_utm, full_obj = vt_hull_locm,  
      ire = shape_utm)
  } 
  
  return(hull_vortess)
  
}

### Calculate local spatial autocorrelation
# function to calculate local Moran's based on conditional permutations
make_local_morans_obj_binom <- function(.x) {
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
  .x[[1]]$coordsdat$locm_significance <- ifelse(.x[[1]]$coordsdat$locm_p_pv <= 0.01, 
                                                print("sig"), print("insig"))
  
  # Get quadrant for hotspots
  .x[[1]]$coordsdat$quad_mean <- attr(locm_p, "quadr")$mean 
  locm_dat <- .x[[1]]$coordsdat
  return(locm_dat)
}

make_local_morans_obj_pos <- function(.x) {
  # Calculate conditional Local Moran's and adjust p-values
  locm_p <- localmoran_perm(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, 
                            nsim = 999, iseed = 123)
  .x[[1]]$coordsdat <- cbind(.x[[1]]$coordsdat, locm_p[, c("Ii", "E.Ii", 
                                                           "Var.Ii", "Z.Ii",
                                                           "Pr(z != E(Ii))",
                                                           "Pr(z != E(Ii)) Sim")
  ])
  .x[[1]]$coordsdat$locm_p_pv <- p.adjust(locm_p[, "Pr(z != E(Ii)) Sim"],
                                          "fdr")
  .x[[1]]$coordsdat$locm_significance <- ifelse(.x[[1]]$coordsdat$locm_p_pv <= 0.01, 
                                                print("sig"), print("insig"))
  
  # Get quadrant for hotspots
  .x[[1]]$coordsdat$quad_mean <- attr(locm_p, "quadr")$mean 
  locm_dat <- .x[[1]]$coordsdat
  return(locm_dat)
}

# function to adjust p-values and check against alpha threshold 
check_psig_locm <- function(x) {
  pv <- x |>
    subset(select = "Pr(z != E(Ii)) Sim", drop = TRUE) 
  
  pva <- function(pv) cbind("none" = pv, 
                            "FDR" = p.adjust(pv, "fdr"), "BY" = p.adjust(pv, "BY"),
                            "Bonferroni" = p.adjust(pv, "bonferroni"))
  pvsp <- pva(pv)
  
  f <- function(x) sum(x < 0.005) # choose 0.005 as "interesting" alpha threshold
  #f <- function(x) sum(x < 0.01)
  check <- apply(pvsp, 2, f)
  return(check)
}

### Calculate local Getis-Ord 
make_getis_ord_obj_binom <- function(.x) {
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

make_getis_ord_obj_pos <- function(.x) {
  # Calculate conditional Getis-Ord and adjust p-values
  dat <- .x[[1]]$coordsdat
  locg_p <- localG_perm(dat$resids2, .x[[1]]$nblistw, 
                        nsim = 999, iseed = 123)
  locg_p_df <- attr(locg_p, "internals")
  dat <- cbind(dat, locg_p_df[, c("Gi", "E.Gi", "Var.Gi", "StdDev.Gi",
                                  "Pr(z != E(Gi))", "Pr(z != E(Gi)) Sim")
  ])
  
  dat$locg_p_pv <- p.adjust(locg_p_df[, "Pr(z != E(Gi)) Sim"],
                            "fdr")
  dat$locg_significance <- ifelse(dat$locg_p_pv <= 0.01, 
                                  print("sig"), print("insig"))
  
  # Get hotspot clusters
  locg_p_clust <- attr(locg_p, "cluster")
  dat <- cbind(dat, locg_p_clust)
  
  locg_dat <- dat
  return(locg_dat)
}

check_psig_locg <- function(x) {
  pv <- x |>
    subset(select = "Pr(z != E(Gi)) Sim", drop = TRUE) 
  
  pva <- function(pv) cbind("none" = pv, 
                            "FDR" = p.adjust(pv, "fdr"), "BY" = p.adjust(pv, "BY"),
                            "Bonferroni" = p.adjust(pv, "bonferroni"))
  pvsp <- pva(pv)
  
  f <- function(x) sum(x < 0.005) # choose 0.005 as "interesting" alpha threshold
  #f <- function(x) sum(x < 0.01)
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
    #geom_sf(data = coords, colour = "black", size = 0.001) +
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
    #geom_sf(data = coords, colour = "black", size = 0.001) +
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
    #geom_sf(data = coords, colour = "black", size = 0.001) +
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

