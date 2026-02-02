#############
## Horse mackeral sdmTMB fit
#############

library(sdmTMB)
library(sdmTMBextra)
library(mgcv)

load("hom_biomass_sdmTMB.RData")
names(dat3)

# Transform coordinates to UTM (ensures constant distances)
dat3 <- add_utm_columns(dat3, c("lon", "lat"), units = "km")

# Check distances between points
dist_mat <- dist(cbind(dat3$X, dat3$Y))

png(filename = paste0(getwd(), "/fits/", "site_distance.png"))
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(dist_mat, 
     freq = TRUE,
     main = "Age 0", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
plot(x = sort(dist_mat), 
     y = (1:length(dist_mat))/length(dist_mat), 
     type = "l",
     xlab = "Distance between sites (km)",
     ylab = "Cumulative proportion")
dev.off()

rm(dist_mat)

################
## Univariate sdmTMB w/ depth as confounding effect
################

sdmTMB_dir <- paste0(getwd(), "/fits")
dir.create(sdmTMB_dir)
  
  ######## Build hurdle model ########
  
  ## Build mesh
  #mesh <- make_mesh(dat3, c("X", "Y"), cutoff = 20)
  
  mesh2 <- make_mesh(dat3, c("X", "Y"),
                     fmesher_func = fmesher::fm_mesh_2d_inla,
                     cutoff = 20, # minimum triangle edge length
                     #max.edge = 500, # inner and outer max triangle lengths
                     offset = 300) # inner and outer border widths
  # think about building other mesh with higher (more coarse) inner cutoff
  
  png(filename = paste0(sdmTMB_dir, "/SPDE_mesh.png"))
  #par(mfrow=c(1,2))
  #plot(mesh)
  #points(dat3[, c("X", "Y")], col = 2, pch = 16, cex = 0.5)
  plot(mesh2)
  points(dat3[, c("X", "Y")], col = 2, pch = 16, cex = 0.5)
  dev.off()
  
  ### NOT RUN ###
  
  ### Build hurdle model
  # GLM-like 
  #m_dlogn <- sdmTMB(
  #  data = dat3,
  #  list(numkm ~ fyear + fgear,
  #       numkm ~ fyear + fgear),
  #  mesh = mesh,
  #  family = delta_lognormal(),
  #  offset = log(dat3$sweptareakm2),
  #  spatial = list("off","off"),
  #  time = "year", 
  #  spatiotemporal = list("off", "off"),
  #  #anisotropy = TRUE,
  #  share_range = TRUE,
  #  silent = FALSE
  #)
  #sanity(m_dlogn)
  
  #m_dlogn_spat <- sdmTMB(
  #  data = dat3,
  #  list(numkm ~ fyear + fgear,
  #       numkm ~ fyear + fgear),
  #  mesh = mesh,
  #  family = delta_lognormal(),
  #  offset = log(dat3$sweptareakm2),
  #  spatial = list("on","on"),
  #  time = "year", 
  #  spatiotemporal = list("off", "off"),
  #  anisotropy = TRUE,
  #  share_range = TRUE,
  #  silent = FALSE
  #)
  #sanity(m_dlogn_spat)
  
  ################
  
  m_dlogn_spattemp <- sdmTMB(
    data = dat3,
    list(numkm ~ fyear + fgear,
         numkm ~ fyear + fgear),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_d <- sdmTMB(
    data = dat3,
    list(numkm ~ fyear + fgear + s(depth_scaled, m = 1),
         numkm ~ fyear + fgear + s(depth_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_sst <- sdmTMB( 
    data = dat3,
    list(numkm ~ fyear + fgear + s(sst_scaled, m = 1),
         numkm ~ fyear + fgear + s(sst_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_sst_t <- sdmTMB( 
    data = dat3,
    list(numkm ~ fyear + fgear + s(sst_scaled, m = 1, by = fyear),
         numkm ~ fyear + fgear + s(sst_scaled, m = 1, by = fyear)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_sst_q <- sdmTMB( 
    data = dat3,
    list(numkm ~ fyear + fgear + s(sst_scaled, m = 1, by = as.factor(quarter)),
         numkm ~ fyear + fgear + s(sst_scaled, m = 1, by = as.factor(quarter))),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_mld <- sdmTMB(
    data = dat3,
    list(numkm ~ fyear + fgear + s(mld_scaled, m = 1),
         numkm ~ fyear + fgear + s(mld_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_mld_t <- sdmTMB(
    data = dat3,
    list(numkm ~ fyear + fgear + s(mld_scaled, m = 1, by = fyear),
         numkm ~ fyear + fgear + s(mld_scaled, m = 1, by = fyear)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_dxmld <- sdmTMB( 
    data = dat3,
    #list(numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(mld_scaled, m = 1) +
    #       s(depth_scaled, mld_scaled, m = 1),
    #     numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(mld_scaled, m = 1) +
    #     s(depth_scaled, mld_scaled, m = 1)),
    list(numkm ~ fyear + fgear + s(depth_scaled, mld_scaled, m = 1),
         numkm ~ fyear + fgear + s(depth_scaled, mld_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  #sanity(m_dlogn_dxmld)
  
  m_dlogn_dxmld2 <- sdmTMB( 
    data = dat3,
    #list(numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(mld_scaled, m = 1) +
    #       s(depth_scaled, mld_scaled, m = 1),
    #     numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(mld_scaled, m = 1) +
    #     s(depth_scaled, mld_scaled, m = 1)),
    list(numkm ~ fyear + fgear + t2(depth_scaled, mld_scaled),
         numkm ~ fyear + fgear + t2(depth_scaled, mld_scaled)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_dxsst <- sdmTMB( 
    data = dat3,
    #list(numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(sst_scaled, m = 1) +
    #     s(depth_scaled, sst_scaled, m = 1),
    #     numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(sst_scaled, m = 1) +
    #       s(depth_scaled, sst_scaled, m = 1)),
    list(numkm ~ fyear + fgear + s(depth_scaled, sst_scaled, m = 1),
         numkm ~ fyear + fgear + s(depth_scaled, sst_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_dxsst2 <- sdmTMB( 
    data = dat3,
    #list(numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(sst_scaled, m = 1) +
    #     s(depth_scaled, sst_scaled, m = 1),
    #     numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(sst_scaled, m = 1) +
    #       s(depth_scaled, sst_scaled, m = 1)),
    list(numkm ~ fyear + fgear + t2(depth_scaled, sst_scaled),
         numkm ~ fyear + fgear + t2(depth_scaled, sst_scaled)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_mldxsst <- sdmTMB( # does not converge
    data = dat3,
    #list(numkm ~ fyear + fgear + s(mld_scaled, m = 1) + s(sst_scaled, m = 1) +
    #     s(mld_scaled, sst_scaled, m = 1),
    #     numkm ~ fyear + fgear + s(mld_scaled, m = 1) + s(sst_scaled, m = 1) +
    #       s(mld_scaled, sst_scaled, m = 1)),
    list(numkm ~ fyear + fgear + s(mld_scaled, sst_scaled, m = 1),
         numkm ~ fyear + fgear + s(mld_scaled, sst_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_mldxsst2 <- sdmTMB( # does not converge
    data = dat3,
    #list(numkm ~ fyear + fgear + s(mld_scaled, m = 1) + s(sst_scaled, m = 1) +
    #     s(mld_scaled, sst_scaled, m = 1),
    #     numkm ~ fyear + fgear + s(mld_scaled, m = 1) + s(sst_scaled, m = 1) +
    #       s(mld_scaled, sst_scaled, m = 1)),
    list(numkm ~ fyear + fgear + t2(mld_scaled, sst_scaled),
         numkm ~ fyear + fgear + t2(mld_scaled, sst_scaled)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  #m_dlogn_dxmld_sst <- sdmTMB( 
    #Error in solve.default(h, g) : 
    #system is computationally singular: reciprocal condition number = 4.38011e-17
  #  data = dat3,
  #  list(numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(mld_scaled, m = 1) +
  #       s(depth_scaled, mld_scaled, m = 1) + s(sst_scaled, m = 1),
  #       numkm ~ fyear + fgear + s(depth_scaled, m = 1) + s(mld_scaled, m = 1) +
  #       s(depth_scaled, mld_scaled, m = 1) + s(sst_scaled, m = 1)),
  #  mesh = mesh,
  #  family = delta_lognormal(),
  #  offset = log(dat3$sweptareakm2),
  #  spatial = list("on","on"),
  #  time = "year", 
  #  spatiotemporal = list("IID", "IID"),
  #  anisotropy = TRUE,
  #  share_range = TRUE,
  #  silent = FALSE
  #)
  
  m_dlogn_dxmld_sst <- sdmTMB( 
    data = dat3,
    list(numkm ~ fyear + fgear + s(depth_scaled, mld_scaled, m = 1) + 
           s(sst_scaled, m = 1),
         numkm ~ fyear + fgear + s(depth_scaled, mld_scaled, m = 1) + 
           s(sst_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
  m_dlogn_dxmld_sst2 <- sdmTMB( 
    data = dat3,
    list(numkm ~ fyear + fgear + t2(depth_scaled, mld_scaled) + 
           s(sst_scaled, m = 1),
         numkm ~ fyear + fgear + t2(depth_scaled, mld_scaled) + 
           s(sst_scaled, m = 1)),
    mesh = mesh2,
    family = delta_lognormal(),
    offset = log(dat3$sweptareakm2),
    spatial = list("on","on"),
    time = "year", 
    spatiotemporal = list("IID", "IID"),
    anisotropy = TRUE,
    share_range = TRUE,
    silent = FALSE
  )
  
 
  ## Quick model comparison
  aic <- AIC(m_dlogn_spattemp, m_dlogn_d, 
             m_dlogn_sst, m_dlogn_mld, m_dlogn_sst_q,
             m_dlogn_sst_t, m_dlogn_mld_t,
             m_dlogn_dxmld, m_dlogn_dxmld2, m_dlogn_dxsst, m_dlogn_dxsst2, 
             m_dlogn_mldxsst, m_dlogn_mldxsst2,
             m_dlogn_dxmld_sst, m_dlogn_dxmld_sst2)
  save(aic, file = paste0(sdmTMB_dir, "/", "aic.RData"))
  
  ######## Model summaries ########
  models <- list(m_dlogn_spattemp, m_dlogn_d, 
                 m_dlogn_sst, m_dlogn_mld, m_dlogn_sst_q,
                 m_dlogn_sst_t, m_dlogn_mld_t,
                 m_dlogn_dxmld, m_dlogn_dxmld2, m_dlogn_dxsst, m_dlogn_dxsst2, 
                 m_dlogn_mldxsst, m_dlogn_mldxsst2,
                 m_dlogn_dxmld_sst, m_dlogn_dxmld_sst2)
 
  model_summaries <- purrr::map(models, ~ list(
    mod <- .x,
    rf = tidy(.x, "ran_pars", conf.int = TRUE),
    ff = tidy(.x, "fixed", conf.int = TRUE),
    loglik = logLik(.x)
  ))
  
  names(model_summaries) <- c("m_dlogn_spattemp", 
                              "m_dlogn_d", "m_dlogn_sst", "m_dlogn_mld",
                              "m_dlogn_sst_q",
                              "m_dlogn_sst_t", "m_dlogn_mld_t",
                              "m_dlogn_dxmld", "m_dlogn_dxmld2", 
                              "m_dlogn_dxsst", "m_dlogn_dxsst2", 
                              "m_dlogn_mldxsst", "m_dlogn_mldxsst2",
                              "m_dlogn_dxmld_sst", "m_dlogn_dxmld_sst2")
  saveRDS(model_summaries, file = paste0(sdmTMB_dir, "/", "final_model_object.rds"))

  
  
  




