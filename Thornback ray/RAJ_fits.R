#############
## Thornback ray sdmTMB fit
#############

library(sdmTMB)
library(sdmTMBextra)
library(dplyr)

load("raj_count_sdmTMB.RData")
names(dat)

summary(dat)
dat |> group_by(survey, year,fquarter,sexmat) |> 
  summarise(n=sum(count)) |> 
  tidyr::pivot_wider(names_from=sexmat,values_from=n) |>
  print(n=50) # not many zeros except in first NIGFS years

# Remove years with 0 catches in all categories
dat <- dat[!(dat$survey == "NIGFS" & dat$year %in% c(2005:2007)), ]

# Proportion of zeros per sex-maturity category
sink(file = paste0(getwd(), "/fits", "zero_proportions.txt"))
dat |> group_by(sexmat) |> 
  summarise(n=sum(count == 0)) |>
  mutate(prop_zero = n/4663) 
sink(file = NULL)

with(dat, table(sexmat, count > 0))

# Transform coordinates to UTM (ensures constant distances)
dat <- add_utm_columns(dat, c("lon", "lat"), units = "km")

# Change order of substrate factor levels
str(dat$substrate)
levels(dat$substrate)
dat$substrate <- factor(dat$substrate, levels = c("MudSand", "CoarseMixed"))

# Check distances between points
dist_list <- list()
groups <- unique(dat$sexmat)

for(smat in groups) {
  sub_dat <- subset(dat, sexmat == smat)
  dist_mat <- dist(cbind(sub_dat$X, sub_dat$Y))
  dist_list[[as.character(paste0("sexmat", smat))]] <- dist_mat
}
png(filename = paste0("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/RAJ/fits/", "site_distance.png"))
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(dist_list$sexmatFimm, 
     freq = TRUE,
     main = "Fimm", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
plot(x = sort(dist_list$sexmatFimm), 
     y = (1:length(dist_list$sexmatFimm))/length(dist_list$sexmatFimm), 
     type = "l",
     xlab = "Distance between sites (km)",
     ylab = "Cumulative proportion")
dev.off()


################
## Univariate sdmTMB w/ depth as confounding effect
################

groups <- unique(dat$sexmat)

for (smat in groups) {
  
  sub_dat <- subset(dat, sexmat == smat)
  
  # Define paths
  sdmTMB_dir <- paste0(getwd(), "/fits")
  if (!file.exists(sdmTMB_dir)) {
    dir.create(sdmTMB_dir)
    print(paste("Created base directory:", sdmTMB_dir))
  }
  
  output_dir <- paste0(sdmTMB_dir, "/", smat)
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
    print(paste("Created output directory for group", smat, ":", output_dir))
  }
  
  ######## Build hurdle model ########
  
  ## Build mesh
  #mesh <- make_mesh(sub_dat, c("X", "Y"), cutoff = 20)
  mesh2 <- make_mesh(sub_dat, c("X", "Y"),
                     fmesher_func = fmesher::fm_mesh_2d_inla,
                     cutoff = 10, # for closely located stations
                     max.edge = c(75, 100), # inner triangle max lengths
                     #max.edge = 300, # inner triangle max lengths
                     #offset = 150 # prevents edge effects
                     offset = 10 # prevents edge effects
                     
                     ) 
  
  
  png(filename = paste0(output_dir, "_SPDE_mesh.png"))
  #plot(mesh)
  plot(mesh2)
  points(sub_dat[, c("X", "Y")], col = 2, pch = 16, cex = 0.5)
  dev.off()
  
  ## negative binomial2
  #m_d_2nb <- sdmTMB(
  #  data = sub_dat,
  #  formula = count ~ fgear + fyear,
  #  #mesh = mesh,
  #  mesh = mesh2,
  #  family = nbinom2(),
  #  offset = sub_dat$log_sweptareakmsqadj,
  #  spatial = "off",
  #  time = "year", 
  #  spatiotemporal = "off",
  #  silent = FALSE
  #)
  
  #m_d_2nbc <- update(m_d_2nb, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1))
  
  #m_d_2nbd <- update(m_d_2nb, formula. = count ~ fgear + fyear + substrate)
  
  #m_d_2nb2 <- sdmTMB(
  #  data = sub_dat,
  #  formula = count ~ fgear + fyear,
  #  #mesh = mesh,
  #  mesh = mesh2,
  #  family = nbinom2(),
  #  offset = sub_dat$log_sweptareakmsqadj,
  #  spatial = "on",
  #  time = "year", 
  #  spatiotemporal = "off",
  #  anisotropy = TRUE,
  #  silent = FALSE
  #)
  
  #m_d_2nb2a <- update(m_d_2nb2, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1))
  
  #m_d_2nb2b <- update(m_d_2nb2, formula. = count ~ fgear + fyear + substrate)
  
  #m_d_2nb2c <- update(m_d_2nb2, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1, by = fgear))
  
  #m_d_2nb2d <- update(m_d_2nb2, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1) + substrate)
  
  #m_d_2nb2e <- update(m_d_2nb2, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1) + s(bottomT_scaled, m = 1))
  
  #m_d_2nb2f <- update(m_d_2nb2, formula. = count ~ fgear + fyear + s(bottomT_scaled, m = 1))
  
  #m_d_2nb2g <- update(m_d_2nb2, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1) + s(bottomT_scaled, m = 1) + substrate)
  
  m_d_2nb3 <- sdmTMB(
    data = sub_dat,
    formula = count ~ fgear + fyear,
    #mesh = mesh,
    mesh = mesh2,
    family = nbinom2(),
    offset = sub_dat$log_sweptareakmsqadj,
    spatial = "on",
    time = "year", 
    spatiotemporal = "IID",
    anisotropy = TRUE,
    silent = FALSE
  )
  
  m_d_2nb4 <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1))
  
  m_d_2nb4a <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1, by = fgear))
  
  #m_d_2nb4b <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1, by = substrate))
  
  m_d_2nb4c <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1) + s(bottomT_scaled, m = 1))
    
  m_d_2nb5 <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(bottomT_scaled, m = 1))
  
  m_d_2nb6 <- update(m_d_2nb3, formula. = count ~ fgear + fyear + substrate)
  
  m_d_2nb7 <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1) + substrate)
  
  m_d_2nb7a <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(middepth_scaled, m = 1, by = fgear) + substrate)
  
  m_d_2nb8 <- update(m_d_2nb3, formula. = count ~ fgear + fyear + s(bottomT_scaled, m = 1) + substrate)
  
  m_d_2nb9 <- update(m_d_2nb3, 
                     formula. = count ~ fgear + fyear + s(bottomT_scaled, m = 1) + s(middepth_scaled, m = 1) + substrate)
  
  m_d_2nb9a <- update(m_d_2nb3, 
                      formula. = count ~ fgear + fyear + s(bottomT_scaled, m = 1) + 
                        s(middepth_scaled, m = 1, by = fgear) + substrate)
  
  ######## Model summaries ########
  models <- list(m_d_2nb3, m_d_2nb4, m_d_2nb4a, 
                 m_d_2nb4c, m_d_2nb5, m_d_2nb6, m_d_2nb7, m_d_2nb7a, 
                 m_d_2nb8, m_d_2nb9, m_d_2nb9a)
  
  model_summaries <- purrr::map(models, ~ list(
    mod <- .x,
    rf = tidy(.x, "ran_pars", conf.int = TRUE),
    ff = tidy(.x, "fixed", conf.int = TRUE),
    aic = AIC(.x),
    loglik = logLik(.x)
    )
  )
  names(model_summaries) <- c("m_d_2nb3", "m_d_2nb4", "m_d_2nb4b", "m_d_2nb4c",
                              "m_d_2nb5", "m_d_2nb6", "m_d_2nb7", 
                              "m_d_2nb7a", "m_d_2nb8", "m_d_2nb9", "m_d_2nb9a")
  saveRDS(model_summaries, file = paste0(output_dir, "/", smat, "_model_object4.rds"))
  
  
} 
  
  


  
