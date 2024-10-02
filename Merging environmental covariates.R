#####
## Merging environmental covariates: temps, depth, habitat
## Date: 02/10/2024
## 
#####

### Load and merge data 

temps_depth <- read.csv("whiting_sst_sbt_depth.csv")
habitat <- read.csv("whiting_bbht.csv")

t_d_h <- merge(temps_depth, habitat)
head(t_d_h)

### Export data

write.csv(t_d_h, file = "environmental_covariates.csv", row.names = FALSE)
