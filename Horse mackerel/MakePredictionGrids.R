#############
## whiting prediction grid
#############

library(sf)
library(raster)
library(terra)
library(dplyr)
library(ggplot2)
library(cowplot)

#path <- 

##### Make grid of survey shapefile -----
load("hom_biomass_sdmTMB.RData")
survey_points <- st_as_sf(dat3, coords = c("lon", "lat"), crs = 4326)
survey_points <- st_transform(survey_points, 2157)

# Load combined survey areas boundary and survey points hull
shps_hull <- st_read(paste0(path, "/plotting aids/", "shps_hull_updated.shp"))
ggplot() + geom_sf(data = shps_hull)

# identify holes
shps <- st_read(paste0(path, "/plotting aids/", "shps_union.shp"))
shps_coords <- as.data.frame(st_coordinates(shps))
holes <- shps_coords[shps_coords$L1 > 1 & 
                       shps_coords$Y > 680000 & shps_coords$X < 850010,]
holes_sf <- st_as_sf(holes, coords = c("X", "Y"), crs = 2157)
ggplot() + geom_sf(data = holes_sf, aes(colour = factor(L1)))
table(holes$L1)
ggplot() + geom_sf(data = holes_sf[holes_sf$L1 == 29,], aes(colour = factor(L1)))
ggplot() + geom_sf(data = holes_sf[holes_sf$L1 == 30,], aes(colour = factor(L1)))
ggplot() + geom_sf(data = holes_sf[holes_sf$L1 == 67,], aes(colour = factor(L1)))
ggplot() + geom_sf(data = holes_sf[holes_sf$L1 == 147,], aes(colour = factor(L1)))
ggplot() + geom_sf(data = holes_sf[holes_sf$L1 == 149,], aes(colour = factor(L1)))
ggplot() + geom_sf(data = holes_sf[holes_sf$L1 == 150,], aes(colour = factor(L1)))
intent_holes <- shps_coords[shps_coords$L1 %in% c(29, 30, 67, 147, 149, 150),] 

outer <- shps_coords[shps_coords$L1 == 1,]
ggplot(data = outer, aes(X, Y)) + geom_point()

wanted_coords <- rbind(intent_holes, outer)
ggplot(wanted_coords, aes(X, Y, colour = factor(L1))) + geom_point()

# make a new shapefile that covers the holes (areas that are not sampled)
wanted_coords <- wanted_coords |>
  mutate(type = case_when(L1 == 29 ~ "hole1", L1 == 30 ~ "hole2",
                          L1 == 67 ~ "hole3", L1 == 147 ~ "hole4",
                          L1 == 149 ~ "hole5", L1 == 150 ~ "hole6",
                          # i see points within the area
                          L1 == 1 ~ "outer"))
line_list <- split(wanted_coords[, c("X", "Y")], wanted_coords$type)
coords_list <- lapply(line_list, as.matrix)
line_st <- purrr::map(coords_list, ~ st_linestring(.x))
holes_st <- line_st[names(line_st) != "outer"]
line_sfc <- st_sfc(holes_st, crs = 2157)
attr <- unique(wanted_coords[, "type"])[1:6]
#attr <- unique(grep("hole", wanted_coords$type, value = TRUE))
i <- match(names(holes_st), attr)
lines_sf <- st_as_sf(data.frame(id = attr[i], length = st_length(line_sfc), 
                                line_sfc))
poly_sf <- lines_sf |> st_cast("POLYGON")

ggplot() + geom_sf(data = poly_sf, aes(colour = id), fill = NA)

# remove untrawled areas from hull
uns_hull <- st_difference(shps_hull, st_union(st_combine(poly_sf)))
ggplot() + geom_sf(data = uns_hull)

# this area contains island that should be removed
neatl <- st_read(paste0(path, "/plotting aids/", "neatl_canvas.shp"))
neatl_utm <- st_transform(neatl, 2157)

# remove land and islands
ref_hull <- st_difference(uns_hull, st_union(st_combine(neatl_utm)))
ggplot() + geom_sf(data = ref_hull, fill = NA, colour = "red")

# Now make a regularly spaced grid
# rasterise
resolution = 5000
r <- rast(vect(ref_hull), resolution = resolution)
rr <- rasterize(vect(ref_hull), r, background = 0)
plot(rr)

grid <- as.data.frame(rr, xy = TRUE, cell = TRUE)
grid$area <- grid$layer * resolution * resolution
grid <- dplyr::filter(grid, area > 0) |> 
  dplyr::select(-layer)

ggplot(grid, aes(x, y, colour = area)) +
  geom_tile(width = resolution, height = resolution, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.01) +
  coord_fixed()

# make coords 
coords <- vect(grid, geom = c("x", "y"), crs = "EPSG:2157")
plot(coords, cex = 0.1)

rm(list = c("line_list", "holes", "holes_sf", 
            "holes_st", "intent_holes", "line_sfc", "line_st", 
            "lines_sf", "outer", "poly_sf", "shps_coords", 
            "wanted_coords", "neatl", "coords_list", 
            "shps_hull", "uns_hull"))

##### Get rasters of environmental covariates -----
### TEMPERATURE & MLD ----
# Get temperature raster
#nc_path <- 
nc_list <- list.files(nc_path, pattern = "*.nc", full.names = TRUE)
#nc_path2 <- 
nc_list2 <- list.files(nc_path2, pattern = "*.nc", full.names = TRUE)
full_nclist <- c(nc_list, nc_list2)

process_nc_fun3 <- function(nc_file, grid_coords, rr) {
  ## Read all files
  r <- rast(nc_file)
  #r <- rast(full_nclist)
  r_utm <- terra::project(r, "EPSG:2157")
  
  # Prepare vector (fish) data
  #grid_coords = coords
  vec_utm <- terra::project(grid_coords, "EPSG:2157")
  
  extr_list <- list()
  
  # Multidimensional raster cube - Loop through all variable names
  #vars <- unique(gsub("_.*", "", names(r_utm)))
  vars <- c("mlotst", "thetao")
  for (var in vars) {
    #var = "thetao"
    #var = "mlotst"
    var_layers <- grep(var, names(r_utm))
    r_var <- subset(r_utm, var_layers)
    
    #extr_out <- NULL
    
    if(grepl("mlotst", var)) {
      
      ## Reorder time dimension
      r_time <- time(r_var)
      r_time_sorted <- order(r_time)
      r_sorted <- r_var[[r_time_sorted]] # reorders layers
      time(r_sorted) <- r_time[r_time_sorted] # explicitly reorders time dim
      
      # Resample to resolution of grid
      resampled <- resample(r_sorted, rr, method = "bilinear")
      
      # Calculate annual mean value (for each cell)
      r_y <- tapp(resampled, "years", mean)
      
      ## Bilinear data extraction
      extracted_val <- terra::extract(r_y, vec_utm,
                                      method = "simple",
                                      xy = TRUE,
                                      cells = TRUE)
      extr_df <- terra::as.data.frame(extracted_val) 
      extr_melt <- reshape2::melt(extr_df,
                                  id = c("cell", "ID", "x", "y"))
      # Data wrangling of output
      colnames(extr_melt)[5] <- "year"
      extr_melt$year <- gsub("y_", "", extr_melt$year)
      
      extr_out <- extr_melt
      
      extr_out$env_var <- print(var) 
      print(extr_out)
      
    } else if (grepl("thetao", var)) {
      
      # Resample to resolution of grid
      resampled <- resample(r_var, rr, method = "bilinear")
      
      ## Calculate monthly SST (mean across upper 20m)
      # requires summarising across the depth layers for a given time 
      # this goes by repetition in date
      monthly_sst <- tapp(resampled, "yearmonths", fun = mean)
      timestamps <- unique(time(r_var))
      time(monthly_sst) <- timestamps
      
      ## Reorder time dimension
      unsorted <- time(monthly_sst)
      r_time_sorted <- order(unsorted)
      r_sorted <- monthly_sst[[r_time_sorted]] # reorders layers
      time(r_sorted) <- unsorted[r_time_sorted] # explicitly reorder time dimension
      
      # Calculate annual mean SST (for each cell)
      r_y <- tapp(r_sorted, "years", mean)
      
      ## Bilinear data extraction
      extracted_val <- terra::extract(r_y, vec_utm,
                                      method = "simple",
                                      xy = TRUE,
                                      cells = TRUE)
      extr_df <- terra::as.data.frame(extracted_val) 
      extr_melt <- reshape2::melt(extr_df,
                                  id = c("cell", "ID", "x", "y"))
      # Data wrangling of output
      colnames(extr_melt)[5] <- "year"
      extr_melt$year <- gsub("y_", "", extr_melt$year)
      
      extr_out <- extr_melt
      
      extr_out$env_var <- print(var) 
      print(extr_out)
    }
    
    # Save in list
    extr_list[[paste(var)]] <- extr_out
    #print(extr_list)
    
  }
  return(extr_list)
}

grid_extr <- process_nc_fun3(nc_file = full_nclist, grid_coords = coords,
                             rr = rr)

#grid_cop <- do.call(rbind,grid_extr)
#save(hom_covs, file = "copernicus_extracted_hom.RData")

range(grid_extr$mlotst$value, na.rm = TRUE)
range(dat3$mld, na.rm = TRUE)
range(grid_extr$thetao$value, na.rm = TRUE)
range(dat3$sst, na.rm = TRUE)

par(mfrow = c(2, 2))
hist(grid_extr$mlotst$value)
hist(grid_extr$thetao$value)
hist(dat3$mld)
hist(dat3$sst)


### DEPTH ----
#nc_path2 <-
depth <- rast(paste0(nc_path2,
                     "GEBCO/gebco_2024_n65.7642_s39.353_w-17.1826_e10.1074.nc"))
depth_utm <- terra::project(depth, "EPSG:2157")

plot(depth_utm)
plot(vect(ref_hull), add = TRUE, border = "white", col = NA)

# Resample data
depth_resampled <- resample(depth_utm, rr, method = "bilinear")
plot(depth_resampled)
plot(vect(ref_hull), add = TRUE, border = "white", col = NA)

# Extract values for survey grid coordinates 
depth_extracted <- extract(depth_resampled, coords, method = "simple",
                           xy = TRUE,
                           cells = TRUE)
head(depth_extracted)
names(depth_extracted)[names(depth_extracted) == "elevation"] <- "depth"

depth_extracted <- depth_extracted |> filter_all(any_vars(!is.na(.)))

range(abs(depth_extracted$depth))
range(dat3$depth)

### GEAR ----
shps_ori <- st_read(paste0(path, "/plotting aids/", "survey_shps.shp"))
shps_utm <- st_transform(shps_ori, 2157)

ggplot(data = survey_points, aes(colour = gear)) +
  geom_sf() +
  facet_wrap(~ fyear)
ggplot() + 
  geom_sf(data = shps_utm, aes(fill = DatasetVer), alpha = 0.4)

# Rasterise vector
# repeat this lines since IGFS has been sampling NIGFS area with GOV in 2003 and 2004
shps2 <- shps_utm |>
  mutate(gear = case_when(
    DatasetVer == "SP-PORC" ~ "PORB",
    DatasetVer == "SP-NORTH" ~ "BAK",
    .default = "GOV"
  ))
gear_r <- rast(vect(shps2), resolution = resolution)
gear_rr <- rasterize(vect(shps2), gear_r, 
                       field = "gear") 
plot(gear_rr)
#plot(vect(exclude_islands), add = TRUE, border = "white", col = NA)

# Resample data
gear_resampled <- resample(gear_rr, rr, method = "bilinear")
plot(gear_resampled)
class(names(gear_resampled))

# Extract values for survey grid coordinates 
gear_extracted <- extract(gear_rr, coords, method = "simple",
                          xy = TRUE,
                          cells = TRUE)
head(gear_extracted)
gear_extracted$gear <- as.character(gear_extracted$gear)

ggplot(gear_extracted, aes(x, y, colour = gear)) +
  geom_tile(width = resolution, height = resolution, fill = NA) +
  scale_colour_viridis_d(direction = -1) +
  geom_point(size = 0.1) +
  coord_fixed()


### MAKE GRID ----
# bring all variable together
sst_grid <- grid_extr$thetao[-ncol(grid_extr$thetao)]
rownames(sst_grid) <- NULL
colnames(sst_grid)[6] <- "sst"
mld_grid <- grid_extr$mlotst[-ncol(grid_extr$mlotst)]
rownames(mld_grid) <- NULL
colnames(mld_grid)[6] <- "mld"
temp <- merge(sst_grid, mld_grid) 
temp <- temp[complete.cases(temp), ] # remove NAs
temp[! complete.cases(temp), ]

temp2 <- merge(depth_extracted[, c("cell", "depth", "x", "y")], 
               gear_extracted[, c("cell", "gear")], by = "cell")
temp2 <- temp2[complete.cases(temp2), ] # remove NAs
temp2[! complete.cases(temp2), ]
all_covs <- merge(temp, temp2)
all_covs <- all_covs[complete.cases(all_covs),]
all_covs <- all_covs[order(all_covs$year),]
all_covs$fyear <- as.factor(all_covs$year)
all_covs$depth <- abs(all_covs$depth)
head(all_covs)

all_covs[duplicated(all_covs),] # check for duplicated rows

ggplot(all_covs[all_covs$year == "2003",], aes(x, y, colour = gear)) +
  geom_point() +
  scale_colour_viridis_d(direction = -1) +
  geom_point(size = 0.1) +
  coord_fixed() +
  facet_wrap( ~ year)

ggplot(all_covs[all_covs$year == "2003",], aes(x, y, colour = sst)) +
  geom_point() +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.1) +
  coord_fixed() +
  facet_wrap( ~ year)

ggplot(all_covs[all_covs$year == "2003",], aes(x, y, colour = mld)) +
  geom_point() +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.1) +
  coord_fixed() +
  facet_wrap( ~ year)

# get WGS84 coordinates and bind to temp grid
coords_wgs84 <- project(coords, "EPSG:4326") 
grid_coords <- terra::as.data.frame(coords_wgs84, geom = "XY")
head(grid_coords)
colnames(grid_coords)[3:4] <- c("lon", "lat")

tmp_grid <- merge(grid_coords, all_covs)
head(tmp_grid)

tmp_grid <- tmp_grid[order(tmp_grid$year),!names(tmp_grid) %in% c("cell", 
                                                                  "area", "ID", 
                                                                  "x", "y")]

# Scale covariates to match model inputs
# altered from https://stackoverflow.com/questions/48815209/create-r-function-
#that-standardizes-multiple-variables-and-creates-new-column
add_scaled <- function(data, vars = colnames(data), ...) {
  data.frame(data,
             setNames(data.frame(scale(data[, vars, drop = FALSE],
                                       center = TRUE, scale = TRUE)),
                      paste0(vars, "_scaled")))
}
pred_grid <- add_scaled(tmp_grid, vars = c("depth", "sst", "mld"))

save(pred_grid, file = "prediction_grid.RData")


##### Compare covariate distributions of prediction grid and observations
load("prediction_grid.RData")
load("hom_biomass_sdmTMB.RData")

pos_obs <- subset(dat3, numkm > 0)
  
p_temp <- ggplot() +
  geom_histogram(data = pred_grid, aes(x = sst_scaled, y = after_stat(density)),
                   colour = 1, fill = "orange", alpha = 0.25) +
  geom_density(data = pred_grid, aes(x = sst_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
  geom_histogram(data = dat3, aes(x = sst_scaled, y = after_stat(density)),
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
  geom_density(data = dat3, aes(x = sst_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
  labs(title = "Age 0", y = "SBT density")
  
p_temp_pos <- ggplot() +
  geom_histogram(data = pred_grid, aes(x = sst_scaled, y = after_stat(density)),
                 colour = 1, fill = "orange", alpha = 0.25) +
  geom_density(data = pred_grid, aes(x = sst_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
  geom_histogram(data = pos_obs, aes(x = sst_scaled, y = after_stat(density)),
                 colour = 1, fill = "deepskyblue", alpha = 0.25) +
  geom_density(data = pos_obs, aes(x = sst_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
  labs(title = "Positive catches: Age 0", y = "SST density")
  
p_depth <- ggplot() +
  geom_histogram(data = pred_grid, aes(x = depth_scaled, y = after_stat(density)),
                 colour = 1, fill = "orange", alpha = 0.25) +
  geom_density(data = pred_grid, aes(x = depth_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
  geom_histogram(data = dat3, aes(x = depth_scaled, y = after_stat(density)),
                 colour = 1, fill = "deepskyblue", alpha = 0.25) +
  geom_density(data = dat3, aes(x = depth_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
  labs(title = "Age 0", y = "Depth density")
  
p_depth_pos <- ggplot() +
  geom_histogram(data = pred_grid, aes(x = depth_scaled, y = after_stat(density)),
                   colour = 1, fill = "orange", alpha = 0.25) +
  geom_density(data = pred_grid, aes(x = depth_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
  geom_histogram(data = pos_obs, aes(x = depth_scaled, y = after_stat(density)),
                 colour = 1, fill = "deepskyblue", alpha = 0.25) +
  geom_density(data = pos_obs, aes(x = depth_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
  labs(title = "Positive catches: Age 0", y = "Depth density")
  
p_mld <- ggplot() +
  geom_histogram(data = pred_grid, aes(x = mld_scaled, y = after_stat(density)),
                   colour = 1, fill = "orange", alpha = 0.25) +
  geom_density(data = pred_grid, aes(x = mld_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
  geom_histogram(data = dat3, aes(x = mld_scaled, y = after_stat(density)),
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
  geom_density(data = dat3, aes(x = mld_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
  labs(title = "Age 0", y = "Depth density")
  
p_mld_pos <- ggplot() +
  geom_histogram(data = pred_grid, aes(x = mld_scaled, y = after_stat(density)),
                 colour = 1, fill = "orange", alpha = 0.25) +
  geom_density(data = pred_grid, aes(x = mld_scaled),
                 lwd = 1, colour = "orange",
                 fill = "orange", alpha = 0.25) +
  geom_histogram(data = pos_obs, aes(x = mld_scaled, y = after_stat(density)),
                   colour = 1, fill = "deepskyblue", alpha = 0.25) +
  geom_density(data = pos_obs, aes(x = mld_scaled),
                 lwd = 1, colour = "deepskyblue",
                 fill = "deepskyblue", alpha = 0.25) +
  labs(title = "Positive catches: Age 0", y = "Depth density")
  
obs_pred_comp <- cowplot::plot_grid(
  p_temp, p_temp_pos,
  p_depth, p_depth_pos,
  p_mld, p_mld_pos,
  ncol = 2, nrow = 3)
cowplot::save_plot(paste0("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/HOM/fits/",
                 "obs_predgrid_comp.jpg"),
          obs_pred_comp, base_asp = 1.6 )





