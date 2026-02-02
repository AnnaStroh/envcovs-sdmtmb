#######
## Extract environmental variables
## Objective: get the environmental variable for each haul point observation
## Benthic/demersal species only 
#######

library(ggplot2); theme_set(theme_bw())
library(viridis)
library(sf)
library(lubridate)
library(terra)
library(tidyterra)
library(dplyr)

######################################################
## Prepare fish data
######################################################

## THORNBACK RAY
raj_path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/RAJ"
load(paste0(raj_path, "/", "raj_cpue_base.RData"))
raj <- raj_final
rm("raj_final")

# make grid and sf object
cols2keep <- c("haulid", "survey", "year", "month", "sexmat", "shootlat", "shootlong")
raj_grid <- raj[, cols2keep]
colnames(raj_grid)[5] <- "age_stage"

# calculate time lag (1 month prior to being caught) for lagged environmental impacts
raj_grid$lagmonth <- raj_grid$month - 1
raj_grid$lag_date <- with(raj_grid, make_date(year, lagmonth, 1))
#raj_grid$dummy_date <- with(raj_grid, make_date(year, month, 1))

# make sf object
raj_vec <- st_as_sf(raj_grid, coords = c("shootlong", "shootlat"), crs = 4326, remove = FALSE) |>
  vect() |> arrange(lag_date)

## WHITING
whg_path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG"
load(paste0(whg_path, "/", "whg_biomass_sdmTMB_temp.RData"))
whg <- dat
rm("dat")

# get month from HH and make grid
load("C:/Users/astroh/Desktop/Chapter 2/VAST/RAJ/raw_HH.RData")
all_HH$Survey[all_HH$Survey == "IE-IGFS"] <- "IGFS"
all_HH$haulid <- with(all_HH, paste0(Survey, Year, Quarter, HaulNo))
whg_temp <- merge(whg, all_HH)
cols2keep <- c("haulid", "survey", "year", "Month", "age", "lat", "lon")
whg_grid <- whg_temp[, cols2keep]
colnames(whg_grid)[5] <- "age_stage"

# calculate time lag (1 month prior to being caught) for lagged environmental impacts
whg_grid$lagmonth <- whg_grid$Month - 1
whg_grid$lag_date <- with(whg_grid, make_date(year, lagmonth, 1))
#whg_grid$dummy_date <- with(whg_grid, make_date(year, Month, 1))

# make sf object
whg_vec <- st_as_sf(whg_grid, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  vect() |> arrange(lag_date)

######################################################
## Prepare spatial cropping
######################################################

## Get survey strata (IGFS and NIGFS)
igfs <- st_read("C:/Users/astroh/Desktop/Chapter 1/Plotting canvases/IGFS_Strata_final.shp")
nigfs <- st_read("C:/Users/astroh/Desktop/Chapter 2/Environmental data/NIGFS/Datras_NIGFS_-_-_-.shp")

## Subset and bind survey strata 
igfs_sub <- igfs[, c("Primary", "Secondary", "geometry")]
colnames(igfs_sub)[1:2] <- c("Strata", "Survey")

nigfs_sub <- nigfs[, c("AreaName", "DatasetVer", "geometry")]
colnames(nigfs_sub)[1:2] <- c("Strata", "Survey")

survey_strata <- rbind(igfs_sub, nigfs_sub)
st_is_valid(survey_strata)
sf_use_s2(FALSE)
survey_vec <- vect(survey_strata)
saveRDS(survey_vec, file = "survey_vec.rds") # needs to be read in with vect()

# Clear unnecessary objects
rm(list = c("all_HH", "whg_temp", "igfs", "nigfs", "igfs_sub", "nigfs_sub", 
            "survey_strata", "raj", "raj_grid", "whg", "whg_grid"))

######################################################
## Temperature, mixed layer thickness (Copernicus)
######################################################

#### Function to process nc files and extract data
process_nc_fun <- function(nc_file, area_vec, timestamps, vars, fish_vec) {
  ## Read all files
  r <- rast(nc_file)
  r_utm <- terra::project(r, "EPSG:2157")
  
  ## Crop all files to survey extent
  area_vec_utm <- terra::project(area_vec, "EPSG:2157")
  r_c <- crop(r_utm, area_vec_utm, mask = TRUE)
  
  ## Subset time dimension
  timestamps <- time(r_c) # original time stamp of SpatRaster
  newtime <- as_date(timestamps) # transform from date:time to date format 
  time(r_c) <- newtime # assign date-transformed time dimension to raster
  
  fish_dates <- unique(fish_vec$lag_date)
  
  index <- newtime[newtime %in% fish_dates] # create date index
  
  r_sub <- r_c[[time(r_c) %in% index]] # subset original raster to fish date range
  
  ## Extract location-specific and date-specific values for each raster variable
  
  # Prepare vector (fish) data
  vec <- fish_vec |>
    # can make coordinate read more smart in future with if statememt
    tidyterra::select(haulid, lat, lon, lag_date) |> # whg
    #tidyterra::select(haulid, shootlat, shootlong, lag_date) |> # raj
    tidyterra::distinct() |>
    tidyterra::arrange(lag_date) # for easier indexing later
  
  # Multidimensional raster cube - Loop through all variable names
  #vars <- c("bottomT", "mlotst", "so", "thetao") # indexing variable (raster chunk) names
  vars <- "bottomT" # indexing variable (raster chunk) names
  extr_list <- list()
  
  for (var in vars) {
    var_layers <- grep(var, names(r_sub))
    r_var <- subset(r_sub, var_layers)
    
    # Bilinear data extraction
    extracted_val <- terra::extract(r_var, vec,
                                    method = 'bilinear',
                                    xy = TRUE,
                                    bind = TRUE,
                                    ID = TRUE)
    if(!var %in% c("so", "thetao")) {
      # Data wrangling 
      # Every date layer extracted for each location, no matter whether it's correct
      lookup <- data.frame(layer_index = seq_along(time(r_var)),
                           layer_name = names(r_var), 
                           lag_date = time(r_var))
      extr_df <- terra::as.data.frame(extracted_val) 
      extr_melt <- reshape2::melt(extr_df,
                                  id = c("haulid", "lat", "lon", "lag_date")
                                  #id = c("haulid", "shootlat", "shootlong", "lag_date")
                                  )
      
      colnames(extr_melt)[5] <- "layer_name"
      extr_out <- merge(extr_melt, lookup) # merging by both layer name AND date
      extr_out <- extr_out[, ! names(extr_out) %in% c("layer_index", "layer_name")]
      extr_out$depth <- NA
      extr_out$env_var <- print(var) 
      
      # Save in list
      extr_list[[paste(var)]] <- extr_out
      
    } else {
      # Adding data wrangling steps to accommodate vars with depth levels
      lookup2 <- data.frame(
        #layer_index = seq_along(time(r_var)), # dropping for computational performance
        layer_name = names(r_var), 
        lag_date = time(r_var))
      lookup2$depth <- sub(".*?(\\d+\\.\\d+).*", "\\1", lookup2$layer_name)
      lookup2$layer_name <- gsub("=\\d+\\.\\d+", "", lookup2$layer_name)
      
      extr_df <- terra::as.data.frame(extracted_val) 
      extr_melt <- reshape2::melt(extr_df,
                                  #id = c("haulid", "lat", "lon", "lag_date")
                                  id = c("haulid", "lat", "lon", "lag_date"))
      colnames(extr_melt)[5] <- "layer_name"
      extr_melt$layer_name <- gsub("\\.\\d+\\.\\d+", "", extr_melt$layer_name)
      extr_out <- merge(extr_melt, lookup2)
      extr_out <- extr_out[, ! names(extr_out) %in% "layer_name"]
      extr_out$env_var <- print(var) 
      
      # Save in list
      extr_list[[paste(var)]] <- extr_out
    }
  }
  
  ## Export as dataframe
  extr_vars <- do.call(rbind, extr_list)
  extr_vars <- extr_vars[complete.cases(extr_vars[, "value"]),]
  # may need to add output of NAs for personal record
  return(extr_vars)
}

## Lapply to use function on all nc files
nc_path <- "C:/Users/astroh/Desktop/Chapter 2/Environmental data/copernicus"
nc_list <- list.files(nc_path, pattern = "*.nc", full.names = TRUE)

values_raj <- lapply(nc_list, process_nc_fun, # apply function over list of netCDF files
                     area_vec = survey_vec, # any area of interest
                     timestamps = timestamps,
                     vars = vars, # variable names 
                     fish_vec = raj_vec # any point-sample for species of interest 
) 
raj_covs <- do.call(rbind, values_raj)
row.names(raj_covs) <- NULL
range(raj_covs$value)
save(raj_covs, file = "raj_bottomT_data.RData")

values_whg <- lapply(nc_list, process_nc_fun, # apply function over list of netCDF files
                     area_vec = survey_vec, # any area of interest
                     timestamps = timestamps,
                     vars = vars, # variable names 
                     fish_vec = whg_vec # any point-sample for species of interest 
)
whg_covs <- do.call(rbind, values_whg)
row.names(whg_covs) <- NULL
range(whg_covs$value)
save(whg_covs, file = "whg_bottomT_data.RData")


# Depending on the species data, rbind may create huge df
# Optionally save as RDS for further processing
#ts_vars <- do.call(rbind, values)
head(ts_vars) 
# or
saveRDS(values_raj, file = "copernicus_extracted_raj2.rds")
saveRDS(values_whg, file = "copernicus_extracted_whg2.rds")
#saveRDS(values_whg2, file = "copernicus_extracted_whg2.rds")

######################################################
## HORSE MACKEREL
######################################################

#### Data preparation ----
hom_path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/HOM"
load(paste0(hom_path, "/", "WHM_recruits_by_haul.RData"))
hom <- haulResDF
rm("haulResDF")

# make sf object
hom_vec <- st_as_sf(hom, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  vect() |> arrange(year)


#### Adjusted function for horse mackerel ----

#### Function to process nc files and extract data
process_nc_fun2 <- function(nc_file, fish_vec) {
  ## Read all files
  r <- rast(nc_file)
  #r <- rast(full_nclist)
  r_utm <- terra::project(r, "EPSG:2157")
  
  # Prepare vector (fish) data
  #fish_vec = hom_vec
  vec <- fish_vec |>
    tidyterra::select(haul, ices_eco, ices_area, lat, lon, year) |> 
    tidyterra::distinct() |>
    tidyterra::arrange(year) # for easier indexing later
  vec_utm <- terra::project(vec, "EPSG:2157")
  
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
      r_sorted <- r_var[[r_time_sorted]]
      time(r_sorted) <- r_time[r_time_sorted]
      
      # Calculate annual mean value (for each cell)
      r_y <- tapp(r_sorted, "years", mean)
      
      ## Bilinear data extraction
      extracted_val <- terra::extract(r_y, vec_utm,
                                      method = 'bilinear',
                                      xy = TRUE,
                                      bind = TRUE,
                                      ID = TRUE)
      extr_df <- terra::as.data.frame(extracted_val) 
      extr_melt <- reshape2::melt(extr_df,
                                  id = c("haul", "ices_eco", "ices_area", 
                                         "lat", "lon", "year"))
      # Data wrangling of output
      colnames(extr_melt)[7] <- "layer_name"
      extr_melt$layer_name <- gsub("y_", "", extr_melt$layer_name)
      
      extr_out <- extr_melt[extr_melt$layer_name == extr_melt$year, ]
      
      # Are records correctly matched?
      check <- all(extr_out$layer_name == extr_out$year)
      if (check == TRUE) {
        print("All years correctly matched.")
      }
      
      extr_out$env_var <- print(var) 
      print(extr_out)
      
    } else if (grepl("thetao", var)) {
      
      ## Calculate monthly SST (mean across upper 20m)
      # requires summarising across the depth layers for a given time 
      # this goes by repetition in date
      monthly_sst <- tapp(r_var, "yearmonths", fun = mean)
      timestamps <- unique(time(r_var))
      time(monthly_sst) <- timestamps
      
      ## Reorder time dimension
      unsorted <- time(monthly_sst)
      r_time_sorted <- order(unsorted)
      r_sorted <- monthly_sst[[r_time_sorted]] # reorders layers
      time(r_sorted) <- unsorted[r_time_sorted] # explicitly reorder time dimension
      
      # Save intermediary product for check-up
      sst_list <- list()
      sst_list[[paste(var)]] <- wrap(monthly_sst, proxy=FALSE)
      
      # Calculate annual mean SST (for each cell)
      r_y <- tapp(r_sorted, "years", mean)
      
      ## Bilinear data extraction
      extracted_val <- terra::extract(r_y, vec_utm,
                                      method = 'bilinear',
                                      xy = TRUE,
                                      bind = TRUE,
                                      ID = TRUE)
      extr_df <- terra::as.data.frame(extracted_val) 
      extr_melt <- reshape2::melt(extr_df,
                                  id = c("haul", "ices_eco", "ices_area", "lat", "lon", "year"))
      # Data wrangling of output
      colnames(extr_melt)[7] <- "layer_name"
      extr_melt$layer_name <- gsub("y_", "", extr_melt$layer_name)
      
      extr_out <- extr_melt[extr_melt$layer_name == extr_melt$year, ]
      
      # Are records correctly matched?
      check <- all(extr_out$layer_name == extr_out$year)
      if (check == TRUE) {
        print("All years correctly matched.")
      }
      extr_out$env_var <- print(var) 
      print(extr_out)
    }
    
    # Save in list
    extr_list[[paste(var)]] <- extr_out
    #print(extr_list)
    
  }
  return(extr_list)
}

nc_path <- "C:/Users/astroh/Desktop/Chapter 2/Environmental data/copernicus/HOM downloads"
nc_list <- list.files(nc_path, pattern = "*.nc", full.names = TRUE)
nc_path2 <- "C:/Users/astroh/Desktop/Chapter 2/Environmental data/copernicus/HOM downloads/2021"
nc_list2021 <- list.files(nc_path2, pattern = "*.nc", full.names = TRUE)
full_nclist <- c(nc_list, nc_list2021)

values_hom2 <- process_nc_fun2(nc_file = full_nclist, fish_vec = hom_vec)

hom_covs <- do.call(rbind,values_hom2)
save(hom_covs, file = "copernicus_extracted_hom.RData")
