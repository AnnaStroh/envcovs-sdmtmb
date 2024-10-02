#####
## Processing sea bottom substrate
## Date: 26/08/2024
## 
#####


library(sf)
library(dplyr)
library(terra)
library(tidyterra)
library(lubridate)
library(ggplot2); theme_set(theme_bw())
library(viridis)

# Read in habitat data  --------------------------------------------------

# read file previously cropped to survey extent

# see Vasquez et al 2015 for info on habitat variable calculation (in GitHub)

substrate_sf <- read_sf("MBBT/survey_MBBT.shp") |>
  select(- c(Secondary, ID3, ID4, ZLEVEL, PERIMETER, AREA, 
             Stn_Target, Hist_75, Rand_25, Test)) 
head(substrate_sf)

unique(substrate_sf$Biozone) # env conditions of the seafloor incl depth
unique(substrate_sf$Substrate) # substrate sediments
unique(substrate_sf$MSFD_BBHT) # energy "climate" + biozone + substrate

# this plot in github!! it takes a while to render
ggplot() +
  geom_sf(data = substrate_sf, aes(fill = MSFD_BBHT) ) 

# Read in whiting data ---------------------------------------------------

fish <- read.csv("CNAA_perHaul_DS_V3.csv") |>
  filter(Survey == "IGFS") |>
  reshape2::melt(id.vars = c("Survey", "Year", "HaulNo", "HaulDur", "Stratum", "ShootLong", "ShootLat", "AreaKmSq")) |>
  mutate(Age = gsub("X", "", variable), 
         NAge = value, .keep = "unused" )
#names(fish) <- tolower(names(fish))
head(fish)


# Expand haul info for tow tracks -----------------------------------------

library(icesDatras)

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
hh_coords <- df_igfs |>
  select(Survey, HaulNo, Year,
         HaulLong, HaulLat) |>
  mutate(Survey = gsub("IE-", "", Survey))
head(hh_coords)

# append haul coords to fish data
fish_coords <- merge(fish, hh_coords, by = c("Survey", "HaulNo", "Year"))
head(fish_coords)

rm(list = c("df_igfs", "hh_coords"))

# Rasterise habitat data  --------------------------------------

# sf to SpatVector
mbbt_v <- vect(substrate_sf)
mbbt_v <- project(mbbt_v, "EPSG:4326")

# Rasterisation
mbbt_r <- rast(mbbt_v) 
mbbt_r <- project(mbbt_r, "EPSG:4326", 
          res = 1/200) # increase spatial resolution to 0.005 degrees

mbbt_rast <- rasterize(mbbt_v, mbbt_r, 
                       field = "MSFD_BBHT", # only use benthic habitat classification (not substrate vars, etc)
                       touches = TRUE)

plot(mbbt_rast, main = "touch = TRUE") # variability across sampling region

# for GitHub upload
#writeRaster(mbbt_rast, filename = "BBHT_raster.tif")
test <- rast("BBHT_raster.tif")
plot(test)

# all mbbt variables 
#nams <- names(mbbt_v)
#allrast <- lapply(nams, function(x) {
 # rasterize(mbbt_v, mbbt_r,
  #          field = x,
   #         touches = TRUE
  #)
#})

# bind all objects
#mbbt_allrast <- do.call("c", allrast)
#mbbt_allrast


# Create tow track objects ------------------------------------------------

# Create haul ids
fish_coords$HaulID <- paste0(fish_coords$Survey, fish_coords$Year, fish_coords$HaulNo)
lines_df <- fish_coords |>
  select(HaulID, ShootLong, ShootLat, HaulLong, HaulLat) |>
  distinct()

# Make lines
# code lines 105-109 from Hans's vignette on spatial basics
lines_list <- split(lines_df[,c('ShootLong','HaulLong','ShootLat','HaulLat')],lines_df$HaulID)
coords_list <- lapply(lines_list,function(x) matrix(as.numeric(x),ncol=2))
lines_ls <- lapply(coords_list,st_linestring)
lines_sfc <- st_as_sfc(lines_ls,crs=4326)
#lines_sfc_irenet <- st_transform(lines_sfc, crs = crs(substrate_sf))
lines_sf <- st_as_sf(data.frame(HaulID=names(lines_ls),lines_sfc))
#lines_sf <- st_as_sf(data.frame(HaulID=names(lines_ls),lines_sfc_irenet))

lines_sf

# Plot lines 
fish_coords_lines <- merge(lines_sf, fish_coords)
head(fish_coords_lines)

ggplot() +
  geom_sf(data = fish_coords_lines) +
  facet_wrap( ~ Year)


# Extract raster data for tow tracks --------------------------------------

lines_v <- vect(lines_sf)
lines_v

### For run: Replace "mbbt_rast" with "test"
                      
# What can I expect?
plot(mbbt_rast, 
     main = "Whiting towing tracks (2003-2018) over broad benthic habitats") 
plot(lines_v, add = TRUE) # some lines without raster coverage

# coverage better on shapefile? (rasterisation forces data loss?)
ggplot() +
  geom_sf(data = substrate_sf, aes(fill = MSFD_BBHT)) + 
  scale_fill_viridis_d() +
  geom_sf(data = lines_sf, colour = "red") +
  labs(title = "Whiting tow tracks (2003-2018) over benthic broad habitats",
       subtitle = "Polygon data")

# Extract raster data along tow tracks
lines_extr <- extractAlong(mbbt_rast, lines_v,
                           ID = TRUE, 
                           xy = TRUE,
                           online = TRUE)
head(lines_extr)


# Merge fish and extracted data -------------------------------------------

# Dealing with categorical raster values in terra
# got number values -> terra transforms categorical layer values to integer
# see https://cran.r-project.org/web/packages/terra/terra.pdf#page=108.08
categories <- as.data.frame(levels(mbbt_rast))
colnames(categories) <- c("MSFD_BBHT", "BBHT") # using "MSFD_BBHT" only as a temporary ID label for merging

#lines_extr2 <- merge(lines_extr, categories)
#head(lines_extr2)
#colnames(lines_extr2)[1] <- "Level_ID"

length(unique(lines_extr$ID)) # less IDs than original vector -> no values to extract?
length(unique(lines_v$HaulID))

lines_df <- fish_coords_lines |>
  # make raster ID ("ID") that connects vector element and raster values extracted for such vector element
  group_by(HaulID) |>
  mutate(ID = cur_group_id()) |>
  # merge extracted values by raster ID into fish sf with line geometry
  merge(lines_extr) |>
  # specify to avoid confusion between coordinates
  dplyr::rename(Rast_Long = x,
                Rast_Lat = y) |>
  # merge in raster value categories
  merge(categories) |>
  # for easier comprehension
  arrange(HaulID) 

head(lines_df)

## Now I need to determine the maximum frequency of habitat occurrence for each haul id 

str(lines_df)

fish_coords_habitats <- lines_df |> 
  group_by(HaulID) |> 
  mutate(MaxBBHT = names(which.max(table(BBHT)))) |>
  select( - c(ID, Rast_Lat, Rast_Long) ) |> 
  distinct()

head(fish_coords_habitats)

table(fish_coords_habitats$MaxBBHT) # 2061 entries without habitat information


# Inspect results ---------------------------------------------------------

# habitat occurrence
# remove ages to avoid counting of duplicates
fish_coords_habitats |>
  select(Year, MaxBBHT) |>
  distinct() |>
  ggplot(aes(MaxBBHT, fill = MaxBBHT)) +
  geom_histogram(stat = "count", position = "dodge", binwidth = 1/10) +
  scale_fill_viridis_d() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("Count of dominant habitats in IGFS whiting 2003-2018")

# biomass over habitat
fish_coords_habitats |>
  select(Year, MaxBBHT, NAge) |>
  distinct() |>
  ggplot(aes(x = MaxBBHT, y = NAge, fill = MaxBBHT)) +
  geom_col(position = "dodge", binwidth = 1/10) +
  scale_fill_viridis_d() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("Dominant habitats x IGFS whiting biomass 2003-2018")


# Export data -------------------------------------------------------------

names(fish_coords_habitats)

fish_habitats <- fish_coords_habitats |>
  select( - c(MSFD_BBHT, HaulID, HaulLong, HaulLat, BBHT) ) |>
  st_drop_geometry()

head(fish_habitats)

write.csv(fish_habitats, file = "whiting_bbht.csv", row.names = FALSE)













