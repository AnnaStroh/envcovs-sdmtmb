#############
## whiting prediction grid
#############

library(sf)
library(raster)
library(terra)
library(dplyr)
library(ggplot2)

#path <- 
setwd(path)
load("whg_biomass_sdmTMB.RData")


##### Make grid of survey shapefile -----
# Load shapefile of survey
survey_shape <- st_read(paste0(path, "survey shapefiles/Datras_IE-IGFS_-_-_-.shp"))
survey_shape <- st_make_valid(survey_shape)
# had to remove a stratum for the model bc of missing substrate data!
sf_use_s2(FALSE)
survey_shape <- st_transform(survey_shape, 2157)
#survey_shape$geometry <- st_geometry(survey_shape) / 1000
survey_points <- st_as_sf(dat[dat$age == 0,], coords = c("lon", "lat"), crs = 4326)
survey_points <- st_transform(survey_points, 2157)
#survey_points$geometry <- st_geometry(survey_points) / 1000
ggplot() + 
  geom_sf(data = survey_shape, aes(fill = Descriptio), alpha = 0.4) +
  geom_sf(data = survey_points)
survey_sub <- survey_shape |>
  filter(! Descriptio %in% c("VIa_Coast, 1-80m", "VIa_Medium, 81-125m",
                             "VIa_Deep, 126-200m", "VIa_Slope, 201-600m",
                             "VIIb_Deep, 126-200m"))
ggplot() + 
  geom_sf(data = survey_sub, aes(fill = Descriptio), alpha = 0.4) +
  geom_sf(data = survey_points)

# rasterise the survey area
resolution <- 5000
r <- rast(vect(survey_sub), resolution = resolution)
rr <- terra::rasterize(vect(survey_sub), r, background = 0)
plot(rr)

grid <- as.data.frame(rr, xy = TRUE, cell = TRUE)
grid$area <- grid$layer * resolution * resolution
grid <- dplyr::filter(grid, area > 0) |> 
  dplyr::select(-layer)

ggplot(grid, aes(x, y, colour = area)) +
  geom_tile(width = resolution, height = resolution, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.1) +
  coord_fixed()

# make coords 
coords <- vect(grid, geom = c("x", "y"), crs = "EPSG:2157")


##### Get rasters of environmental covariates -----
### TEMPERATURE ----
# Get temperature raster
#nc_path <- 
nc_list <- list.files(nc_path, pattern = "*.nc", full.names = TRUE)

rasts <- rast(nc_list)
r_utm <- terra::project(rasts, "EPSG:2157")

timestamps <- time(r_utm) # original time stamp of SpatRaster
index <- unique(as.POSIXct(as.Date(survey_points$lag_date)))
index <- sort(index) 
r_sub <- r_utm[[time(r_utm) %in% index]]

var_layers <- grep("bottomT", names(r_sub), value = TRUE)
r_var <- r_sub[[names(r_sub) %in% var_layers]]

# since I add several rasters with different time ranges, layer names are duplicated
# this will cause confusion during value extraction 
newnames <- NULL
for (i in 1:nlyr(r_var)) {
  newnames[i] <- paste0("bottomT", i)
}
names(r_var) <- newnames

plot(r_var[[1]])
plot(vect(survey_sub), add = TRUE, border = "white", col = NA)

# Resample temperature raster to survey raster resolution
sbt_resampled <- resample(r_var, rr, method = "bilinear")
plot(sbt_resampled[[1]])
plot(vect(survey_sub), add = TRUE, border = "white", col = NA)

# Extract values for survey grid coordinates 
sbt_extracted <- extract(sbt_resampled, coords, method = "simple",
                         xy = TRUE,
                         cells = TRUE)
sbt_extr_melt <- reshape2::melt(sbt_extracted,
                                id = c("cell", "ID", "x", "y"))
head(sbt_extr_melt)
sbt_lookup <- data.frame(
  variable = names(sbt_resampled),
  lag_date = time(sbt_resampled))
sbt_grid <- merge(sbt_extr_melt, sbt_lookup)
names(sbt_grid)[names(sbt_grid) == "value"] <- "bottomT"
head(sbt_grid)

### DEPTH ----
#nc_path2 <- 
depth <- rast(paste0(nc_path2,
                     "GEBCO/gebco_2024_n65.7642_s39.353_w-17.1826_e10.1074.nc"))
depth_utm <- terra::project(depth, "EPSG:2157")

plot(depth_utm)
plot(vect(survey_sub), add = TRUE, border = "white", col = NA)

# Resample data
depth_resampled <- resample(depth_utm, rr, method = "bilinear")
plot(depth_resampled)
plot(vect(survey_sub), add = TRUE, border = "white", col = NA)

# Extract values for survey grid coordinates 
depth_extracted <- extract(depth_resampled, coords, method = "simple",
                           xy = TRUE,
                           cells = TRUE)
head(depth_extracted)
names(depth_extracted)[names(depth_extracted) == "elevation"] <- "middepth"

range(abs(depth_extracted$middepth))
range(dat$middepth)

### SUBSTRATE ----
#nc_path2 <- 
substr <- st_read(paste0(nc_path2,
                         "Folk/Multiscale - folk 5/seabed_substrate_1m.shp"))
substr_utm <- st_transform(substr, 2157)
ggplot() + geom_sf(data = substr, aes(fill = folk_5cl_t)) + 
  geom_sf(data = survey_sub, colour = "white", fill = NA)

# Crop all files to survey extent
substr_utm <- st_make_valid(substr_utm)
survey_sub <- st_make_valid(survey_sub)
substr_c <- st_intersection(substr_utm, survey_sub)
ggplot() + geom_sf(data = substr_c , aes(fill = folk_5cl_t)) +
  geom_sf(data = survey_points , colour = "black")

# Rasterise vector
substr_r <- rast(vect(substr_c), resolution = resolution)
substr_rr <- rasterize(vect(substr_c), substr_r, 
                       field = "folk_5cl") # take numeric substrate variable

# Resample data
substr_resampled <- resample(substr_rr, rr, method = "bilinear")
plot(substr_resampled)
plot(vect(survey_sub), add = TRUE, border = "white", col = NA)

# Extract values for survey grid coordinates 
substr_extracted <- extract(substr_resampled, coords, method = "simple",
                            xy = TRUE,
                            cells = TRUE)
head(substr_extracted)

# make a lookup table for substrate levels and merge in
substr_lookup <- data.frame(
  folk_5cl = substr_utm$folk_5cl,
  folk_5cl_t = substr_utm$folk_5cl_t, 
  substrate_chr = gsub("[0-9]\\. ", "", substr_utm$folk_5cl_t, perl = TRUE)
)
substr_lookup <- distinct(substr_lookup)

substr_grid <- merge(substr_extracted, substr_lookup)
head(substr_grid)


### MAKE GRID ----
# bring all variable together
all_covs <- merge(sbt_grid, depth_extracted[, c("cell", "middepth")], by = "cell") |>
  merge(substr_grid)
all_covs <- all_covs[complete.cases(all_covs),]
head(all_covs)
all_covs$fyear <- lubridate::year(lubridate::ymd(all_covs$lag_date))
all_covs$year <- all_covs$fyear
all_covs$substrate <- all_covs$folk_5cl
all_covs$fsubstrate <- as.factor(paste0("f", all_covs$folk_5cl))
all_covs <- all_covs[!all_covs$fsubstrate %in% "f4",]
all_covs <- droplevels(all_covs)
all_covs$substrate2 <- as.factor(all_covs$substrate_chr)
all_covs$middepth <- abs(all_covs$middepth)


# get WGS84 coordinates and bind to temp grid
coords_wgs84 <- terra::project(coords, "EPSG:4326") 
grid_coords <- terra::as.data.frame(coords_wgs84, geom = "XY")
head(grid_coords)
colnames(grid_coords)[3:4] <- c("lon", "lat")

tmp_grid <- merge(grid_coords, all_covs)
head(tmp_grid)

# Make an additional grid with average annual temperature 
agg <- tmp_grid |>
  group_by(lat, lon, year, cell) |>
  summarise(meanbottomT = mean(bottomT))

av_grid <- tmp_grid |>
  merge(agg[, c("year", "cell", "meanbottomT")], by = c("cell", "year")) |>
  select(-c(cell, area, ID, x, y, variable, folk_5cl_t, bottomT, lag_date)) |>
  distinct()
head(av_grid)

# Make a grid with all observations for individual lag dates
tmp_grid <- tmp_grid[order(tmp_grid$year),!names(tmp_grid) %in% c("cell", 
                                                                  "area", "ID", 
                                                                  "x", "y", 
                                                                  "variable", 
                                                                  "folk_5cl_t")]
head(tmp_grid)

# Scale covariates to match model inputs
# altered from https://stackoverflow.com/questions/48815209/create-r-function-
#that-standardizes-multiple-variables-and-creates-new-column
add_scaled <- function(data, vars = colnames(data), ...) {
  data.frame(data,
             setNames(data.frame(scale(data[, vars, drop = FALSE],
                                       center = TRUE, scale = TRUE)),
                      paste0(vars, "_scaled")))
}
pred_grid <- add_scaled(tmp_grid, vars = c("middepth","bottomT"))
av_pred_grid <- add_scaled(av_grid, vars = c("middepth","meanbottomT"))

save(pred_grid, file = paste0(getwd(), "/", "prediction_grid.RData"))
save(av_pred_grid, file = paste0(getwd(), "/", "average_prediction_grid.RData"))

# Check grid
load("prediction_grid.RData")
head(pred_grid)
load("whg_biomass_sdmTMB.RData")
head(dat)

ggplot() +
  geom_point(data = pred_grid, aes(lon, lat)) +
  #scale_colour_viridis_c(direction = -1) #+
  geom_point(data = dat[dat$age ==0,], aes(lon, lat), colour = "green") +
  labs(title = "Prediction grid (black) overlayed with whiting age 0 sampling points (green)")
  coord_fixed()
ggsave(filename = paste0(getwd(), "/", "prediction_grid_whg.jpg"), plot = last_plot())


##### Compare covariate distributions of prediction grid and observations
ages <- sort(unique(dat$age))[1:3]
plot_list = list()
for (a in 0:2) {
  obs <- subset(dat, age == a)
  pos_obs <- subset(obs, biomass > 0)
  
  p_temp <- ggplot() +
    geom_histogram(data = pred_grid, aes(x = bottomT_scaled, y = after_stat(density)),
                   colour = 1, fill = "orange", alpha = 0.25) +
    geom_density(data = pred_grid, aes(x = bottomT_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
    geom_histogram(data = obs, aes(x = sbt_scaled, y = after_stat(density)),
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
    geom_density(data = obs, aes(x = sbt_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
    labs(title = paste("Age", a), y = "SBT density")

  p_temp_pos <- ggplot() +
    geom_histogram(data = pred_grid, aes(x = bottomT_scaled, y = after_stat(density)),
                   colour = 1, fill = "orange", alpha = 0.25) +
    geom_density(data = pred_grid, aes(x = bottomT_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
    geom_histogram(data = pos_obs, aes(x = sbt_scaled, y = after_stat(density)),
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
    geom_density(data = pos_obs, aes(x = sbt_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
    labs(title = paste("Positive catches: Age", a), y = "SBT density")
  
  p_depth <- ggplot() +
    geom_histogram(data = pred_grid, aes(x = middepth_scaled, y = after_stat(density)),
                   colour = 1, fill = "orange", alpha = 0.25) +
    geom_density(data = pred_grid, aes(x = middepth_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
    geom_histogram(data = obs, aes(x = middepth_scaled, y = after_stat(density)),
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
    geom_density(data = obs, aes(x = middepth_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
    labs(title = paste("Age", a), y = "Depth density")
  
  p_depth_pos <- ggplot() +
    geom_histogram(data = pred_grid, aes(x = middepth_scaled, y = after_stat(density)),
                   colour = 1, fill = "orange", alpha = 0.25) +
    geom_density(data = pred_grid, aes(x = middepth_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
    geom_histogram(data = pos_obs, aes(x = middepth_scaled, y = after_stat(density)),
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
    geom_density(data = pos_obs, aes(x = middepth_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
    labs(title = paste("Positive catches: Age", a), y = "Depth density")
  
  p_substr <- ggplot() +
    geom_histogram(data = pred_grid, aes(x = fsubstrate), stat="count",
                   colour = 1, fill = "orange", alpha = 0.25) +
    geom_histogram(data = obs, aes(x = fsubstrate), stat="count",
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
    labs(title = paste("Age", a), y = "Substrate density")
 
   p_substr_pos <- ggplot() +
    geom_histogram(data = pred_grid, aes(x = fsubstrate), stat="count",
                   colour = 1, fill = "orange", alpha = 0.25) +
    geom_histogram(data = pos_obs, aes(x = fsubstrate), stat="count",
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
    labs(title = paste("Positive catches: Age", a), y = "Substrate density")
  
  plot_list[[as.character(a)]] <- list("SBT" = list("binom" = p_temp, 
                                                 "pos" = p_temp_pos),
                                       "Depth" = list("binom" = p_depth, 
                                                   "pos" = p_depth_pos),
                                       "Substrate" = list("binom" = p_substr, 
                                                       "pos" = p_substr_pos))

}

obs_pred_comp <- cowplot::plot_grid(
  plot_list[[1]]$Depth$binom, plot_list[[1]]$Depth$pos,
  plot_list[[2]]$Depth$binom, plot_list[[2]]$Depth$pos,
  plot_list[[3]]$Depth$binom, plot_list[[3]]$Depth$pos,
  ncol = 2, nrow = 3)
save_plot(paste0("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG/fits/",
                 "obs_predgrid_depth_comp.jpg"),
          obs_pred_comp, base_asp = 1.6 )

obs_pred_comp2 <- cowplot::plot_grid(
  plot_list[[1]]$Substrate$binom, plot_list[[1]]$Substrate$pos,
  plot_list[[2]]$Substrate$binom, plot_list[[2]]$Substrate$pos,
  plot_list[[3]]$Substrate$binom, plot_list[[3]]$Substrate$pos,
  ncol = 2, nrow = 3)
save_plot(paste0("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG/fits/",
                 "obs_predgrid_substr_comp.jpg"),
          obs_pred_comp2, base_asp = 1.6 )


