#############
## whiting fits
#############

library(sdmTMB)
library(sdmTMBextra)
library(mgcv)
library(dplyr)

load("whg_biomass_sdmTMB.RData")
names(dat)

summary(dat)
dat |> group_by(year,age) |> 
  summarise(n=sum(biomass)) |> 
  tidyr::pivot_wider(names_from=age,values_from=n) |>
  print(n=50) 

# Remove the entries with habitat "outliers"
table(dat$fsubstrate)
dat <- dat[!dat$fsubstrate %in% "f4",]
dat <- droplevels(dat)
dat$substrate2 <- as.factor(dat$substrate_chr)

# Transform coordinates to UTM (ensures constant distances)
dat <- add_utm_columns(dat, c("lon", "lat"), units = "km") # CRS = 32629

# Subset and order age
dat <- subset(dat, age < 3)
dat <- dat[order(dat$age),]


# Check distances between points
dist_list <- list()

for(a in 0:2) {
  sub_dat <- subset(dat, age == a)
  dist_mat <- dist(cbind(sub_dat$X, sub_dat$Y))
  dist_list[[as.character(paste0("age", a))]] <- dist_mat
}
png(filename = "site_distance.png")
par(mfrow = c(3,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(dist_list$age0, 
       freq = TRUE,
       main = "Age 0", 
       xlab = "Distance between sites (km)",
       ylab = "Frequency")
plot(x = sort(dist_list$age0), 
       y = (1:length(dist_list$age0))/length(dist_list$age0), 
       type = "l",
       xlab = "Distance between sites (km)",
       ylab = "Cumulative proportion")
hist(dist_list$age1, 
       freq = TRUE,
       main = "Age 1", 
       xlab = "Distance between sites (km)",
       ylab = "Frequency")
plot(x = sort(dist_list$age1), 
       y = (1:length(dist_list$age1))/length(dist_list$age1), 
       type = "l",
       xlab = "Distance between sites (km)",
       ylab = "Cumulative proportion")
hist(dist_list$age2, 
       freq = TRUE,
       main = "Age 2", 
       xlab = "Distance between sites (km)",
       ylab = "Frequency")
plot(x = sort(dist_list$age2), 
       y = (1:length(dist_list$age2))/length(dist_list$age2), 
       type = "l",
       xlab = "Distance between sites (km)",
       ylab = "Cumulative proportion")
dev.off()

################
## Univariate sdmTMB
################

ages <- unique(dat$age)

for (ages in 0:2) {
  
  #ages = 0
  sub_dat <- subset(dat, age == ages)
  
  # Define paths
  sdmTMB_dir <- paste0(getwd(), "/fits")
  if (!file.exists(sdmTMB_dir)) {
    dir.create(sdmTMB_dir)
    print(paste("Created base directory:", sdmTMB_dir))
  }
  
  output_dir <- paste0(sdmTMB_dir, "/", ages)
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
    print(paste("Created output directory for group", ages))
  }
  
  
  ### Build mesh
  #mesh <- make_mesh(sub_dat, c("X", "Y"), cutoff = 20)
  mesh2 <- make_mesh(sub_dat, c("X", "Y"),
                     fmesher_func = fmesher::fm_mesh_2d_inla,
                     cutoff = 6, # cutoff for similar sampling locations
                     max.edge = c(75, 100), # inner and outer max triangle lengths
                     #max.edge = 75, # inner and outer max triangle lengths
                     #max.edge = c(50, 100), # inner and outer max triangle lengths
                     #offset = 150 # roughly equivalent to spatial range
                     offset = 10 # roughly equivalent to spatial range
                     )
  
  png(filename = paste0(output_dir, "_SPDE_mesh.png"))
  plot(mesh2)
  points(sub_dat[, c("X", "Y")], col = 2, pch = 16, cex = 0.5)
  dev.off()
  
  ### Build hurdle model
  # GLM-like 
  #  m_dlogn <- sdmTMB(
  #  data = sub_dat,
  #  list(biomass ~ fyear,
  #       biomass ~ fyear),
  #  mesh = mesh2,
  #  family = delta_lognormal(),
  #  offset = log(sub_dat$areakmsqadj),
  #  spatial = list("off","off"),
  #  time = "year", 
  #  spatiotemporal = list("off", "off"),
  #  #anisotropy = TRUE,
  #  share_range = TRUE,
  #  silent = FALSE
  #  )
    #sanity(m_dlogn)
  
  # Spatial model
  #  m_dlogn_spat <- sdmTMB(
  #    data = sub_dat,
  #    #list(biomass ~ fyear,
  #    #     biomass ~ fyear),
  #    list(biomass ~ fyear,
  #         biomass ~ fyear),
  #    mesh = mesh2,
  #    family = delta_lognormal(),
  #    offset = log(sub_dat$areakmsqadj),
  #    spatial = list("on","on"),
  #    time = "year", 
  #    spatiotemporal = list("off", "off"),
  #    anisotropy = TRUE,
  #    share_range = TRUE,
  #    silent = FALSE
  #  )
    #sanity(m_dlogn_spat)
    
  # Spatiotemporal models
    m_dlogn_spattemp <- sdmTMB(
      data = sub_dat,
      list(biomass ~ fyear,
           biomass ~ fyear),
      mesh = mesh2,
      family = delta_lognormal(),
      offset = log(sub_dat$areakmsqadj),
      spatial = list("on","on"),
      time = "year", 
      spatiotemporal = list("IID", "IID"),
      anisotropy = TRUE,
      share_range = TRUE,
      silent = FALSE
    )
    #sanity(m_dlogn_spattemp2)
    
    #m_dlogn_spattemp_iso <- sdmTMB(
    #  data = sub_dat,
    #  list(biomass ~ fyear,
    #     biomass ~ fyear),
    #  mesh = mesh2,
    #  family = delta_lognormal(),
    #  offset = log(sub_dat$areakmsqadj),
    #  spatial = list("on","on"),
    #  time = "year", 
    #  spatiotemporal = list("IID", "IID"),
    #  anisotropy = FALSE,
    #  share_range = TRUE,
    #  silent = FALSE
    #)
  
  # Environmental covariate models
    m_dlogn3 <- sdmTMB(
      data = sub_dat,
      list(biomass ~ fyear + s(middepth_scaled, m = 1),
           biomass ~ fyear + s(middepth_scaled, m = 1)),
      mesh = mesh2,
      family = delta_lognormal(),
      offset = log(sub_dat$areakmsqadj),
      spatial = list("on","on"),
      time = "year", 
      spatiotemporal = list("IID", "IID"),
      anisotropy = TRUE,
      share_range = TRUE,
      silent = FALSE
    )
    #sanity(m_dlogn3)
    
    m_dlogn4 <- sdmTMB( # does not converge
      data = sub_dat,
      list(biomass ~ fyear + s(sbt_scaled, m = 1),
           biomass ~ fyear + s(sbt_scaled, m = 1)),
      mesh = mesh2,
      family = delta_lognormal(),
      offset = log(sub_dat$areakmsqadj),
      spatial = list("on","on"),
      time = "year", 
      spatiotemporal = list("IID", "IID"),
      anisotropy = TRUE,
      share_range = TRUE,
      silent = FALSE
    )
    #sanity(m_dlogn4)
    
    m_dlogn5 <- sdmTMB(
      data = sub_dat,
      list(biomass ~ fyear + substrate2,
           biomass ~ fyear + substrate2),
      mesh = mesh2,
      family = delta_lognormal(),
      offset = log(sub_dat$areakmsqadj),
      spatial = list("on","on"),
      time = "year", 
      spatiotemporal = list("IID", "IID"),
      anisotropy = TRUE,
      share_range = TRUE,
      silent = FALSE
    )
    #sanity(m_dlogn5)
    
    m_dlogn6 <- sdmTMB( 
      data = sub_dat,
      list(biomass ~ fyear + s(middepth_scaled, m = 1) + substrate2,
           biomass ~ fyear + s(middepth_scaled, m = 1) + substrate2),
      mesh = mesh2,
      family = delta_lognormal(),
      offset = log(sub_dat$areakmsqadj),
      spatial = list("on","on"),
      time = "year", 
      spatiotemporal = list("IID", "IID"),
      anisotropy = TRUE,
      share_range = TRUE,
      silent = FALSE
    )
    #sanity(m_dlogn6)
    
  #  m_dlogn7 <- sdmTMB( 
  #    data = sub_dat,
  #    list(biomass ~ fyear + s(sbt_scaled, m = 1, by = fyear),
  #         biomass ~ fyear + s(sbt_scaled, m = 1, by = fyear)),
  #    mesh = mesh2,
  #    family = delta_lognormal(),
  #    offset = log(sub_dat$areakmsqadj),
  #    spatial = list("on","on"),
  #    time = "year", 
  #    spatiotemporal = list("IID", "IID"),
  #    anisotropy = TRUE,
  #    share_range = TRUE,
  #    silent = FALSE
  #  )
    #sanity(m_dlogn7)
    
    m_dlogn8 <- sdmTMB( 
      data = sub_dat,
      list(biomass ~ fyear + s(middepth_scaled, m = 1) + s(sbt_scaled, m = 1),
           biomass ~ fyear + s(middepth_scaled, m = 1) + s(sbt_scaled, m = 1)),
      mesh = mesh2,
      family = delta_lognormal(),
      offset = log(sub_dat$areakmsqadj),
      spatial = list("on","on"),
      time = "year", 
      spatiotemporal = list("IID", "IID"),
      anisotropy = TRUE,
      share_range = TRUE,
      silent = FALSE
    )
    #sanity(m_dlogn8)
    
  #  m_dlogn9 <- sdmTMB( 
  #    data = sub_dat,
  #    list(biomass ~ fyear + s(middepth_scaled, sbt_scaled, m = 1),
  #         biomass ~ fyear + s(middepth_scaled, sbt_scaled, m = 1)),
  #    mesh = mesh2,
  #    family = delta_lognormal(),
  #    offset = log(sub_dat$areakmsqadj),
  #    spatial = list("on","on"),
  #    time = "year", 
  #    spatiotemporal = list("IID", "IID"),
  #    anisotropy = TRUE,
  #    share_range = TRUE,
  #    silent = FALSE
  #  )
    
    m_dlogn10 <- sdmTMB( 
      data = sub_dat,
      list(biomass ~ fyear + s(middepth_scaled, m = 1) + s(sbt_scaled, m = 1) + substrate2,
           biomass ~ fyear + s(middepth_scaled, m = 1) + s(sbt_scaled, m = 1) + substrate2),
      mesh = mesh2,
      family = delta_lognormal(),
      offset = log(sub_dat$areakmsqadj),
      spatial = list("on","on"),
      time = "year", 
      spatiotemporal = list("IID", "IID"),
      anisotropy = TRUE,
      share_range = TRUE,
      silent = FALSE
    )
    #sanity(m_dlogn8)
    
 #   m_dlogn11 <- sdmTMB( 
 #     data = sub_dat,
 #     list(biomass ~ fyear + s(middepth_scaled, sbt_scaled, m = 1) + substrate2,
 #          biomass ~ fyear + s(middepth_scaled, sbt_scaled, m = 1) + substrate2),
 #     mesh = mesh2,
 #     family = delta_lognormal(),
 #     offset = log(sub_dat$areakmsqadj),
 #     spatial = list("on","on"),
 #     time = "year", 
 #     spatiotemporal = list("IID", "IID"),
 #     anisotropy = TRUE,
 #     share_range = TRUE,
 #     silent = FALSE
 #   )
    
  ## Quick model comparison
  # remove m_dgamma_spattemp1 for speed
    #aic <- AIC(m_dlogn, m_dlogn_spat, m_dlogn_spattemp, m_dlogn_spattemp_iso,
    #           m_dlogn3, m_dlogn4, m_dlogn5, m_dlogn6, m_dlogn7, m_dlogn8, 
    #           m_dlogn9, m_dlogn10, m_dlogn11
    #           )
    #save(aic, file = paste0(output_dir, "/", ages, "_aic.RData"))
    #models <- list(m_dlogn, m_dlogn_spat, m_dlogn_spattemp, m_dlogn_spattemp_iso,
    #               m_dlogn3, m_dlogn4, m_dlogn5, m_dlogn6, m_dlogn7, m_dlogn8, 
    #               m_dlogn9, m_dlogn10, m_dlogn11)
    
    aic <- AIC(m_dlogn_spattemp, #m_dlogn_spattemp_iso, 
               m_dlogn3, m_dlogn4, 
               m_dlogn5, m_dlogn6, m_dlogn8, m_dlogn10)
    save(aic, file = paste0(output_dir, "/", ages, "_aic.RData"))
    models <- list(m_dlogn_spattemp, #m_dlogn_spattemp_iso, 
                   m_dlogn3, m_dlogn4, 
                   m_dlogn5, m_dlogn6, m_dlogn8, m_dlogn10)
    
    
    ######## Model summaries ######## 
    model_summaries <- purrr::map(models, ~ list(
      mod <- .x,
      aic = AIC(.x),
      #caic = cAIC(.x),
      loglik = logLik(.x)
    ))
    #names(model_summaries) <- c("m_dlogn", "m_dlogn_spat", "m_dlogn_spattemp",
    #                            "m_dlogn_spattemp_iso", "m_dlogn3", "m_dlogn4", 
    #                            "m_dlogn5", "m_dlogn6", "m_dlogn7", "m_dlogn8", 
    #                            "m_dlogn9", "m_dlogn10", "m_dlogn11"
    #                            )
    #names(model_summaries) <- c("m_dlogn_spattemp",
    #                            "m_dlogn3", "m_dlogn6")
    
    names(model_summaries) <- c("m_dlogn_spattemp", #"m_dlogn_spattemp_iso", 
                                "m_dlogn3", "m_dlogn4", "m_dlogn5", "m_dlogn6", 
                                "m_dlogn8", "m_dlogn10")
    saveRDS(model_summaries, file = paste0(output_dir, "/", ages, "_model_object.rds"))
  
}

# time-varying tests

ages <- unique(dat$age)

for (ages in 0:2) {
  
  ages = 2
  sub_dat <- subset(dat, age == ages)
  
  # Define paths
  sdmTMB_dir <- paste0(getwd(), "/fits")
  if (!file.exists(sdmTMB_dir)) {
    dir.create(sdmTMB_dir)
    print(paste("Created base directory:", sdmTMB_dir))
  }
  
  output_dir <- paste0(sdmTMB_dir, "/", ages)
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
    print(paste("Created output directory for group", ages))
  }
  
  
  ### Build mesh
  #mesh <- make_mesh(sub_dat, c("X", "Y"), cutoff = 20)
  mesh2 <- make_mesh(sub_dat, c("X", "Y"),
                     fmesher_func = fmesher::fm_mesh_2d_inla,
                     cutoff = 6, # cutoff for similar sampling locations
                     max.edge = 75, # inner and outer max triangle lengths
                     #max.edge = c(50, 100), # inner and outer max triangle lengths
                     offset = 150 # roughly equivalent to spatial range
  ) 
  
  png(filename = paste0(output_dir, "_SPDE_mesh.png"))
  #plot(mesh)
  plot(mesh2)
  points(sub_dat[, c("X", "Y")], col = 2, pch = 16, cex = 0.5)
  dev.off()
  
  # Spatiotemporal models
  m_dlogn_spattemp <- sdmTMB(
    data = sub_dat,
    list(biomass ~ fyear,
         biomass ~ fyear),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  # Environmental covariate models
  m_dlogn4a <- sdmTMB( # does not converge
    data = sub_dat,
    list(biomass ~ fyear + s(sbt_scaled, by = fyear),
         biomass ~ fyear + s(sbt_scaled, by = fyear)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  sanity(m_dlogn4a)
  
  m_dlogn4b <- sdmTMB( # does not converge
    data = sub_dat,
    list(biomass ~ fyear + s(sbt_scaled, year),
         biomass ~ fyear + s(sbt_scaled, year)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  sanity(m_dlogn4b)
  
  m_dlogn4c <- sdmTMB( # does not converge
    data = sub_dat,
    list(biomass ~ fyear,
         biomass ~ fyear),
    time_varying = ~ 0 + sbt_scaled,
    time_varying_type = "rw",
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  sanity(m_dlogn4c)
  
  models <- list(m_dlogn_spattemp, m_dlogn_spattemp_iso, m_dlogn3, m_dlogn4, 
                 m_dlogn5, m_dlogn6, m_dlogn8, m_dlogn10)
  models <- list(m_dlogn4c)
  
  
  ######## Model summaries ######## 
  model_summaries <- purrr::map(models, ~ list(
    mod <- .x,
    aic = AIC(.x),
    #caic = cAIC(.x),
    loglik = logLik(.x)
  ))
  
  names(model_summaries) <- c("m_dlogn4")
  saveRDS(model_summaries, file = paste0(output_dir, "/", ages, "_timevarying_model_object.rds"))
  
}

# excluding intercepts

ages <- unique(dat$age)

for (ages in 0:2) {
  
  ages = 0
  sub_dat <- subset(dat, age == ages)
  
  # Define paths
  sdmTMB_dir <- paste0(getwd(), "/fits")
  if (!file.exists(sdmTMB_dir)) {
    dir.create(sdmTMB_dir)
    print(paste("Created base directory:", sdmTMB_dir))
  }
  
  output_dir <- paste0(sdmTMB_dir, "/", ages)
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
    print(paste("Created output directory for group", ages))
  }
  
  
  ### Build mesh
  #mesh <- make_mesh(sub_dat, c("X", "Y"), cutoff = 20)
  mesh2 <- make_mesh(sub_dat, c("X", "Y"),
                     fmesher_func = fmesher::fm_mesh_2d_inla,
                     cutoff = 6, # cutoff for similar sampling locations
                     max.edge = 75, # inner and outer max triangle lengths
                     #max.edge = c(50, 100), # inner and outer max triangle lengths
                     offset = 150 # roughly equivalent to spatial range
  ) 
  
  png(filename = paste0(output_dir, "_SPDE_mesh.png"))
  #plot(mesh)
  plot(mesh2)
  points(sub_dat[, c("X", "Y")], col = 2, pch = 16, cex = 0.5)
  dev.off()
  
  # Spatiotemporal models
  m_dlogn_spattemp <- sdmTMB(
    data = sub_dat,
    list(biomass ~ 0 + fyear,
         biomass ~ 0 + fyear),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = TRUE
  )
  
  m_dlogn_spattemp2 <- sdmTMB(
    data = sub_dat,
    list(biomass ~ fyear,
         biomass ~ fyear),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = TRUE
  )
  
  m_dlogn3 <- sdmTMB(
    data = sub_dat,
    list(biomass ~ 0 + fyear + s(middepth_scaled, m = 1),
         biomass ~ 0 + fyear + s(middepth_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = TRUE
  )
  #sanity(m_dlogn3)
  
  m_dlogn5 <- sdmTMB(
    data = sub_dat,
    list(biomass ~ 0 + fyear + substrate2,
         biomass ~ 0 + fyear + substrate2),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = TRUE
  )
  #sanity(m_dlogn5)
  
  m_dlogn6 <- sdmTMB( 
    data = sub_dat,
    list(biomass ~ 0 + fyear + s(middepth_scaled, m = 1) + substrate2,
         biomass ~ 0 + fyear + s(middepth_scaled, m = 1) + substrate2),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(sub_dat$areakmsqadj),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = TRUE
  )

  models <- list(m_dlogn_spattemp, m_dlogn3, m_dlogn5, m_dlogn6)

  ######## Model summaries ######## 
  model_summaries <- purrr::map(models, ~ list(
    mod <- .x,
    aic = AIC(.x)
  ))
  
  names(model_summaries) <- c("m_dlogn_spattemp_noint", "m_dlogn3_noint", 
                              "m_dlogn5_noint", "m_dlogn6_noint")
  saveRDS(model_summaries, 
          file = paste0(output_dir, "/", ages, "_noint_model_object.rds"))
  
}









