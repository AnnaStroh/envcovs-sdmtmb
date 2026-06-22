######
### Remove benthic/demersal sampling points outside observed habitat
### 
######

library(sf)
library(ggplot2); theme_set(theme_bw())
library(dplyr)


## Get fish data
whg_path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG"
load(paste0(whg_path, "/", "whg_cnaa.RData")) # update data below
#whg <- read.csv(paste0(whg_path, "/", "CNAA_perHaul_WGCSE2025.csv"))
whg <- dat
rm("dat")
whg_sf <- st_as_sf(whg, coords = c("shootlong", "shootlat"), 
                   crs = 4326, remove = FALSE)
whg_utm <- st_transform(whg_sf, 2157)

raj_path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/RAJ"
load(paste0(raj_path, "/", "RAJ_CPUE_HH_SexMat.RData"))
raj <- raj_sexmats
rm(raj_sexmats)
raj_sf <- st_as_sf(raj, coords = c("ShootLong", "ShootLat"), 
                   crs = 4326, remove = FALSE)
raj_utm <- st_transform(raj_sf, 2157)

## Read in substrate file
shp_path <- "C:/Users/astroh/Desktop/Chapter 2/Environmental data/Folk"
s <- st_read(paste0(shp_path, "/", 
                    "Multiscale - folk 5/seabed_substrate_250k.shp"))
s <- s[, c("folk_5cl", "folk_5cl_t", "geometry")] # stick with 5-folk
s_val <- st_make_valid(s)
all(st_is_valid(s_val))
s_utm <- st_transform(s_val, 2157)
plot(s)

## Clip points with existing habitat data
s_utm$folk_5cl_t <- gsub("6. No data at this level", NA, s_utm$folk_5cl_t) 
table(is.na(s_utm$folk_5cl_t))

s_nomissing <- dplyr::filter(s_utm, ! is.na(s_utm$folk_5cl_t))

whg_subst <- st_intersection(whg_utm, s_nomissing)
raj_subst <- st_intersection(raj_utm, s_nomissing)

#### ---------------------------------------------
## Format and properly attach habitat data to fish data 
#### ---------------------------------------------

# whiting
unique(whg_subst$folk_5cl_t)
whg_subst$folk_5 <- gsub("[0-9]\\. ", "", whg_subst$folk_5cl_t, perl = TRUE)

ggplot() + 
  geom_sf(data = whg_subst, aes(colour = folk_5)) +
  facet_wrap(~ year)

# ray
unique(raj_subst$folk_5cl_t)
raj_subst$folk_5 <- gsub("[0-9]\\. ", "", raj_subst$folk_5cl_t, perl = TRUE)

ggplot() + 
  geom_sf(data = raj_subst, aes(colour = folk_5)) +
  facet_wrap(~ Year)

### For my record: how many records are being removed? 
whg_removed <- nrow(whg) - nrow(whg_subst) # 3150
whg_removed_perc <- (whg_removed/nrow(whg))*100 #18.4
100 - whg_removed_perc # 81.6
removed_rows <- which(! rownames(whg) %in% rownames(whg_subst))
removed <- subset(whg, rownames(whg) %in% removed_rows)
removed |>
  group_by(age) |>
  summarise(n())
gt_tib <- gt::gt(tib)
gt::gtsave(gt_tib, "age_removed.png", 
           path = "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG")
 
raj_removed <- nrow(raj) - nrow(raj_subst) 
raj_removed_perc <- (raj_removed/nrow(raj))*100 # 8.754462
100 - raj_removed_perc # 91.24554
removed_rows <- which(! rownames(raj) %in% rownames(raj_subst))
removed <- subset(raj, rownames(raj) %in% removed_rows)
tib <- removed |> # how many positive catches removed in each category?
  group_by(SexMat, HLCpueSexMat > 0) |>
  summarise(n())
gt_tib <- gt::gt(tib)
gt::gtsave(gt_tib, "sexmats_removed.png", 
       path = "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/RAJ")

#### ---------------------------------------------
## Get and attach bathymetry-based depth
#### ---------------------------------------------

load("C:/Users/astroh/Desktop/Chapter 2/Environmental data/BathyDepth_HH_IGFS_NIGFS.RData")
depth <- depth_df[, c("HaulID", "Survey", "Year", "ShootLong", "ShootLat",
                      "MidDepth", "MeanDepth", "MedianDepth")]
colnames(depth) <- tolower(colnames(depth))

# whiting
whg_substr_depth <- merge(whg_subst, depth)
ggplot() +
  geom_sf(data = whg_substr_depth[whg_substr_depth$nage > 0,], 
          aes(colour = folk_5, size = nage)) +
  facet_wrap(~ age)

# ray
colnames(raj_subst) <- tolower(colnames(raj_subst))
raj_substr_depth <- merge(raj_subst, depth)
ggplot() +
  geom_sf(data = raj_substr_depth[raj_substr_depth$hlcpuesexmat > 0,], 
          aes(colour = folk_5, size = hlcpuesexmat)) +
  facet_wrap(~ sexmat)

#### ---------------------------------------------
## Get and attach strata names (RAJ only)
#### ---------------------------------------------

## Load survey strata
igfs <- st_read("C:/Users/astroh/Desktop/Chapter 1/Plotting canvases/IGFS_Strata_final.shp")
names(igfs)
nigfs <- st_read("C:/Users/astroh/Desktop/Chapter 2/Environmental data/NIGFS/Datras_NIGFS_-_-_-.shp")
names(nigfs)

ggplot() +
  geom_sf(data = igfs, fill = "blue", colour = "white") +
  geom_sf(data = nigfs, fill = "green", colour = "white")

## Subset and bind survey strata 
igfs_sub <- igfs[, c("Primary", "Secondary", "geometry")]
colnames(igfs_sub)[1:2] <- c("Strata", "Survey")

nigfs_sub <- nigfs[, c("AreaName", "DatasetVer", "geometry")]
colnames(nigfs_sub)[1:2] <- c("Strata", "Survey")

survey_strata <- rbind(igfs_sub, nigfs_sub)
st_is_valid(survey_strata)
sf_use_s2(FALSE)

# Calculate stratum area for later use 
survey_strata$Stratum_Area <- units::set_units(st_area(survey_strata), "km2")

area_calc <- survey_strata[, c("Strata", "Stratum_Area")]
colnames(area_calc)[1] <- "Primary"
area_calc <- st_drop_geometry(area_calc)
test <- merge(igfs, area_calc)
with(test, plot(AreaKm2, units::drop_units(Stratum_Area))) # looks fine

## Retrieve strata names for each haul 
survey_strata_utm <- st_transform(survey_strata, 2157)
fish_strata <- st_intersection(raj_substr_depth, survey_strata_utm)
fish_strata2 <- st_intersection(whg_substr_depth, survey_strata_utm)

ggplot() + 
  geom_sf(data = fish_strata2, aes(colour = Survey)) + # IGFS in NI in 2003-2005
  geom_sf(data = survey_strata_utm, fill = NA, colour = "black") +
  scale_colour_viridis_d()

names(fish_strata)
colnames(fish_strata)[27:28] <- c("Survey_sf", "Stratum_AreaKm2")
colnames(fish_strata) <- tolower(colnames(fish_strata))
colnames(fish_strata2)[23:24] <- c("Survey_sf", "Stratum_AreaKm2")
colnames(fish_strata2) <- tolower(colnames(fish_strata2))

## There seems to be a small data loss
ggplot() + 
  geom_sf(data = fish_strata, aes(colour = survey)) + # IGFS in NI in 2003-2005
  geom_sf(data = survey_strata, fill = NA, colour = "black") +
  scale_colour_viridis_d()
# as expected some hauls fall outside survey domain
#nrow(fish_sf) - nrow(fish_strata) # 52
52/4 # 13 hauls loss


#### ---------------------------------------------
## Format and attach gear information (IGFS gears for whiting) 
#### ---------------------------------------------
load("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/RAJ/raw_HH.RData")
colnames(all_HH) <- tolower(colnames(all_HH))
all_HH$survey[all_HH$survey == "IE-IGFS"] <- "IGFS"
all_HH$haulid <- with(all_HH, paste0(survey, year, quarter, haulno))
HH_sub <- all_HH[, c("haulid", "gearex")]

fish_strata2_gear <- merge(fish_strata2, HH_sub)
fish_strata2_gear$gearex[fish_strata2_gear$gearex == "S"] <- "GOV"
fish_strata2_gear$gearex[fish_strata2_gear$gearex == "I2"] <- "GOVmod"

#### ---------------------------------------------
## Save data 
#### ---------------------------------------------

## Final formatting
whg_final <- st_drop_geometry(whg_substr_depth)
raj_final <- st_drop_geometry(fish_strata)

#colnames(whg_final) <- tolower(colnames(whg_final))
colnames(raj_final) <- tolower(colnames(raj_final))


## Write out data
save(whg_final, 
     file = paste0(whg_path, "/", "whg_cnaa_base.RData"))

save(raj_final, 
     file = paste0(raj_path, "/", "raj_cpue_base.RData"))




## Make stratum plot for my documentation
survey_strata$Mid <- st_centroid(survey_strata$geometry)
ggplot() + 
  geom_sf(data = survey_strata, aes(fill = Strata)) +
  geom_sf_label(data = survey_strata, aes(label = Strata)) +
  scale_fill_viridis_d()
ggsave(filename = "IGFS_NIGFS_strata_plot.png", plot = last_plot(), dpi = 300)


