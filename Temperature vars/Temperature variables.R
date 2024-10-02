#####
## Processing temperature variables
## Date: 26/08/2024 | Author: AS
## 
#####

library(rerddap)
library(terra)
library(tidyterra)
library(sf)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw())

# Read in SST and SBT data from rerddap -----------------------------------------------

urlBase <- "https://erddap.marine.ie/erddap/" 
pars <- c("sea_surface_temperature", "sea_bottom_temperature")
#times <- c("2012-01-15T00:00:00Z", "2024-07-15T00:00:00Z")
times <- c("2012-01-15T00:00:00Z", "2018-12-15T00:00:00Z") # working with Coilin's vast data
lats <- c(50, 56.5)
lons <- c(-13, -4)
dataInfo <- rerddap::info("IMI_Model_Stats", url = urlBase)

dats <- griddap(dataInfo,
                longitude = lons,
                latitude = lats,
                time = times,
                fields = pars,
                url = urlBase, 
                read = FALSE,
                store = disk()) # nc file

file.size(dats$summary$filename)/1e6 ## in MB

r <- rast(dats$summary$filename)

r

# key to understanding this raster
varnames(r)
names(r)

plot(r[[1:84]]) # sst
plot(r[[85:nlyr(r)]]) # sbt

# Cut raster to survey size -----------------------------------------------

path2 <- "C:/Users/astroh/Desktop/Chapter 1/Plotting canvases/"
survey <- read_sf(paste0(path2, 'IGFS_Strata_final.shp')) |>
  vect()
#survey
plot(survey) # weird artefact that looks like a point

# Crop raster
r_cropped <- terra::crop(r, survey, 
                         mask = TRUE)
r_cropped
#plot(r_cropped[[1:84]])
#plot(r_cropped[[85:nlyr(r)]])

rm(list = c("r")) # remove since I am working now with cropped

# Load survey data (VAST input) -----------------------------------------------

# NOTE: read data is Coilin's whiting data subsetted to IGFS

fish <- read.csv("CNAA_perHaul_DS_V3.csv") |>
  filter(Survey == "IGFS") |>
  reshape2::melt(id.vars = c("Survey", "Year", "HaulNo", "HaulDur", "Stratum", "ShootLong", "ShootLat", "AreaKmSq")) |>
  mutate(Age = gsub("X", "", variable), 
         NAge = value, .keep = "unused" )
#names(fish) <- tolower(names(fish))
head(fish)

### no dates (y/m/d) in VAST data, only year -> we need to attach dates

# Bind sampling dates to VAST input ---------------------------------------

library(icesDatras)

## get haul data
#years_igfs <- getSurveyYearList("IE-IGFS")
years_igfs <- 2003:2018 # limit for now due to Coilin's data
hh_igfs_t <- list()
for(y in years_igfs){
  print(y)
  tmp <- getHHdata(survey = "IE-IGFS", year = y, quarter = 4)
  hh_igfs_t[[paste(y)]] <- tmp
  rm(tmp)
}
df_igfs <- do.call(rbind, hh_igfs_t)
row.names(df_igfs) <- NULL
df_igfs <- subset(df_igfs, HaulVal == "V")

names(df_igfs)
hh_dates <- df_igfs |>
  select(Year, Month, Day,
         HaulNo, ShootLong, ShootLat) |>
  mutate(SamplingDate = make_date(Year, Month, Day))
head(hh_dates)

## merge dates with VAST input
fish_dates <- merge(fish, hh_dates, by = c("Year", "HaulNo", "ShootLat", "ShootLong"))
head(fish_dates)

min(fish_dates$SamplingDate) # "2003-10-27"

## calculate "TemperatureDate" -> sampling date of each station - 1 month (to account for time-lag in temperature effects)
fish_dates2 <- fish_dates |>
  mutate(TempMonth = Month - 1, 
         TempMonth = as.integer(TempMonth),
         TempDate = as_date(paste(Year, TempMonth, Day, sep = "-")),
         HaulID = paste0(Survey, Year, HaulNo)) |> # for later merging
  arrange(TempDate) |>
  filter(Year > 2011) # limit for now due to available temperature data

## some NAs -> calculated dates do not exist in the calender
invalid_dates <- fish_dates2 |>
  filter(is.na(TempDate)) |>
  select(SamplingDate) |>
  distinct()
invalid_dates # all 31st day of the October - no 31st in September in these year

## go with nearest day date in the month (Sept 30th) for all NAs in temperature date
fish_dates2$TempDate <- if_else(is.na(fish_dates2$TempDate), 
                                as_date(paste(fish_dates2$Year, fish_dates2$TempMonth, 30, sep = "-")),
                                fish_dates2$TempDate)
min(fish_dates2$TempDate)

# plot for reference
ggplot(data = fish_dates2, aes(ShootLong, ShootLat)) +
  geom_point() +
  facet_wrap( ~ Year) +
  labs(title = "IGFS stations for whiting 2012-2018")

# Reduce time to survey months --------------------------------------------------

times <- ymd(time(r_cropped))
years <- as.factor(year(times))
year_date <- data.frame(years, times)

sept_dec <- year_date |>
  group_by(years) |>
  mutate(month = month(times)) |>
  filter(month >= 9)

survey_dates <- sept_dec$times
time(r_cropped) <- times
new_time <- times[times %in% survey_dates]
r3 <- r_cropped[[time(r_cropped) %in% sept_dec$times]]

# check multi-band raster after subsetting
r3
#varnames(r3) <- c("sea_surface_temperature", "sea_bottom_temperature") #not possible
names(r3)

# raster summaries
nlyr(r3)/2 # equal numbers of layers on both bands
# sea surface temperature
sst_df <- global(r3[[1:28]], "mean", na.rm = TRUE)
sst_df$time <- as.POSIXct(time(r3[[1:28]]))

ggplot(sst_df, aes(x = time, y = mean)) +
  geom_line(col = "blue") +
  geom_point(col = "blue") +
  geom_hline(yintercept = mean(sst_df$mean), colour = "blue3") +
  ggtitle("Mean surface temperature") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  xlab("Time") +
  ylab("Mean temperature (\u00B0C)")

# sea bottom temperature
sbt_df <- global(r3[[29:nlyr(r3)]], "mean", na.rm = TRUE)
sbt_df$time <- as.POSIXct(time(r3[[29:nlyr(r3)]]))

ggplot(sbt_df, aes(x = time, y = mean)) +
  geom_line(col = "blue") +
  geom_point(col = "blue") +
  geom_hline(yintercept = mean(sbt_df$mean), colour = "blue3") +
  ggtitle("Mean bottom temperature") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  xlab("Time") +
  ylab("Mean temperature (\u00B0C)")


# Extract temperature for sampling stations -----------------------------------------------

## to match temperature dates and raster dates, I need dates in %Y-%m format
#disregard the difference of temp dates to the middle of the month for now
#nvm monthly mean!!

# temperature time
fish_dates2$TempDateRaster <- as_date(paste(fish_dates2$Year, 
                                            fish_dates2$TempMonth, 
                                            15, sep = "-"))
min(fish_dates2$TempDateRaster)
#head(fish_dates2)

# make SpatVector (over sf route) for extraction
igfs_sv <- fish_dates2 |>
  #select(Year, HaulID, TempDateRaster, TempDate, Long, Lat) |>
  select(Year, HaulID, ShootLong, ShootLat) |>
  vect(geom = c("ShootLong", "ShootLat"), crs = "wgs84", keepgeom = TRUE) |>
  distinct()
min(igfs_sv$TempDateRaster)
# extract raster data for survey stations 
igfs_temp_extr2 <- terra::extract(r3, igfs_sv,
                                 method = "simple",
                                 bind = TRUE)
head(igfs_temp_extr2)

# from wide to long
igfs_temp_extr2_df <- terra::as.data.frame(igfs_temp_extr2)
igfs_temp_melt2 <- reshape2::melt(igfs_temp_extr2_df, 
                                  id = c("Year", "HaulID", 
                                         "ShootLong", "ShootLat"))
head(igfs_temp_melt2)

colnames(igfs_temp_melt2)[5:6] <- c("Layer", "Temperature")


# Merge extracted data with VAST fish data --------------------------------

# prepare merging dfs
layer_time <- data.frame(Layer = names(r3), # layer information
                         TempDateRaster = time(r3)) 

fish$HaulID <- paste0(fish$Survey, fish$Year, fish$HaulNo) # vast data w/ HaulID

# merge data in
whg_temps <- igfs_temp_melt2 |>
  
  # merge in layer information
  merge(layer_time) |>
  
  # create new variable to track sst and sbt
  mutate(TempCat = ifelse( grepl("bottom", Layer), print("sbt"), print("sst") )) |>
  
  # merge in fish information incl sampling and temperature dates
  merge(fish_dates2, by = c("HaulID", "TempDateRaster")) |>
  
  # tidy up
  select(Survey, HaulID, Year.x, HaulNo, HaulDur, ShootLong.x, ShootLat.x, Stratum,
         SamplingDate, TempDate, TempDateRaster, TempCat, Temperature, 
         AreaKmSq, Age, NAge) |>
  dplyr::rename(Year = Year.x,
                ShootLong = ShootLong.x,
                ShootLat = ShootLat.x) |>
  arrange(TempDate, Age) |>
  distinct() #remove duplicates

# Add depth data  -------------------------------------------

head(df_igfs)

hh_igfs_depth <- df_igfs |>
  mutate(Survey = gsub("IE-", "", Survey),
         HaulID = paste0(Survey, Year, HaulNo),
         SamplingDate = make_date(Year, Month, Day)) |>
  select(HaulID, SamplingDate, Depth) 
head(hh_igfs_depth)

whg_temps_depth <- whg_temps |>
  merge(hh_igfs_depth, by = c("HaulID", "SamplingDate")) |>
  select(Survey, Year, HaulNo, HaulDur, ShootLat, ShootLong, Stratum, 
         SamplingDate, TempDate, TempCat, Temperature, Depth, 
         AreaKmSq, Age, NAge) |>
  arrange(TempDate, Age)
head(whg_temps_depth)


# Data cleanup and export -------------------------------------------------

getwd()
write.csv(whg_temps_depth, 
          file = "whiting_sst_sbt_depth.csv", 
          row.names = FALSE)
























