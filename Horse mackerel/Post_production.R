#############
## HOM fits
## diagnostics, plots, index
#############

library(sdmTMB)
library(sdmTMBextra)
library(visreg)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(DHARMa)
library(patchwork)
library(purrr)

#path <- 
#plot_path <- 

# load functions
source("funs_hom.R")

### Check AIC ------

load(paste0(path, "fits/aic.RData"))
aic
min(aic$AIC) # m_dlogn_dxsst

sink(file = paste0(path, "fits/", "aic.txt"))
aic
sink(file = NULL)

mods <- readRDS(paste0(path, "fits/final_model_object.rds"))
names(mods)

### Read data ------
mods$m_dlogn_dxsst[[1]]

sanity(mods$m_dlogn_spattemp[[1]])
sanity(mods$m_dlogn_mldxsst[[1]])

mod_spattemp <- mods$m_dlogn_spattemp[[1]]
mod_mldxsst <- mods$m_dlogn_mldxsst[[1]] # best

sink(file = paste0(path, "fits/","base_model.txt"))
mod_spattemp

print("SDREPORT")
mod_spattemp$sd_report
sink(file = NULL)

sink(file = paste0(path, "fits/","cov_model.txt"))
mod_mldxsst

print("SDREPORT")
mod_mldxsst$sd_report
sink(file = NULL)

########## -------------------------
### GET OUTPUTS AND PRODUCTS
########## -------------------------
load("prediction_grid.RData")
head(pred_grid)
pred_grid <- add_utm_columns(pred_grid, c("lon", "lat"), units = "km")
grid_yrs <- data.frame(X = pred_grid$X,
                       Y = pred_grid$Y,
                       sst_scaled = pred_grid$sst_scaled, 
                       mld_scaled = pred_grid$mld_scaled,
                       fyear = as.factor(pred_grid$fyear),
                       fgear = as.factor(pred_grid$gear),
                       year = as.numeric(pred_grid$year))
grid_yrs <- distinct(grid_yrs)

rm(list = c("mods", "pred_grid"))

base_index <- get_index_split(mod_spattemp, 
                              newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                              bias_correct = TRUE, 
                              nsplit = 6)
save(base_index, file = paste0(path, "fits/","base_index.RData"))

rm(list = c( "base_index", "mod_spattemp"))
cov_index <- get_index_split(mod_mldxsst, 
                             newdata = grid_yrs, bias_correct = TRUE, 
                             nsplit = 12)
save(cov_index, file = paste0(path, "fits/","cov_index.RData"))

load(paste0(path, "fits/","base_index.RData"))
base_index$model <- "base"
cov_index$model <- "base+mld*sst"
all_indices <- rbind(base_index, cov_index)

### Visualise indices
index_hom <- ggplot(all_indices, aes(year, est, colour = model), size=2) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  #scale_color_manual(values = c("red", "blue", "green")) +
  #scale_fill_manual(values = c("red", "blue", "green")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(fill = "Model", colour = "Model") +
  xlab('Year') + ylab('Relative Index of Abundance') +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15)) +
  ggtitle("Abundance index by model approach: Age 0 western horse mackerel")
ggsave(filename = paste0(path, "fits/","indices.jpg"), plot = index_hom)
saveRDS(index_hom, file = paste0(plot_path, "hom_index_ribbon.rds"))

index_hom2 <- ggplot(all_indices, aes(year, est, colour = model, linetype = model), size=2) + 
  geom_line(linewidth = 2) +
  scale_color_manual(values = c("grey", "black")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('Relative Index of Abundance') +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = c(.1,.95),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "fits/","indices2.jpg"), plot = index_hom2)
saveRDS(index_hom2, file = paste0(plot_path, "hom_index.rds"))

all_indices$CI_width <- with(all_indices, upr - lwr)
head(all_indices)

index_hom_error <- ggplot(all_indices, aes(year, CI_width, colour = model, linetype = model), size=2) + 
  geom_line(linewidth = 2) +
  #scale_color_manual(values = c("red", "blue", "green")) +
  scale_color_manual(values = c("grey", "black")) +
  #scale_color_manual(values = c("red", "blue")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('CI width of Index') +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = c(.1,.95),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "fits/","hom_index_CI.jpg"), plot = index_hom_error)
saveRDS(index_hom_error, file = paste0(plot_path, "hom_index_CI.rds"))

index_hom_se <- ggplot(all_indices, aes(year, log(se), colour = model, linetype = model), size=2) + 
  geom_line(linewidth = 2) +
  scale_color_manual(values = c("grey", "black")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('Log-SE of the Index') +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = c(.25,.95),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "fits/","indices_se.jpg"), plot = index_hom_se)
saveRDS(index_hom_se, file = paste0(plot_path, "hom_index_se.rds"))

#save(all_indices, file = paste0(path, "fits/", "allindices.RData"))
load(paste0(path, "fits/", "allindices.RData"))

### Make plot of error ------
base_rf_enc <- tidy(mod_spattemp, model = 1, effects = "ran_pars")
base_rf_pos <- tidy(mod_spattemp, model = 2, effects = "ran_pars")
t_base <- rbind(base_rf_enc, base_rf_pos)
colnames(t_base)[1] <- "model_comp"
t_base$model_comp[t_base$model_comp == 1] <- "Binomial component"
t_base$model_comp[t_base$model_comp == 2] <- "Lognormal component"
t_base$model <- "base"
head(t_base)

cov_rf_enc <- tidy(mod_mldxsst, model = 1, effects = "ran_pars")
cov_rf_pos <- tidy(mod_mldxsst, model = 2, effects = "ran_pars")
t_cov <- rbind(cov_rf_enc, cov_rf_pos)
colnames(t_cov)[1] <- "model_comp"
t_cov$model_comp[t_cov$model_comp == 1] <- "Binomial component"
t_cov$model_comp[t_cov$model_comp == 2] <- "Lognormal component"
t_cov$model <- "base+mld*sst"
head(t_cov)

all_rfs <- rbind(t_base, t_cov)

all_rfs_sub <- subset(all_rfs, term %in% c("sigma_O", "sigma_E"))
all_rfs_sub$term[all_rfs_sub$term == "sigma_O"] <- "spatial"
all_rfs_sub$term[all_rfs_sub$term == "sigma_E"] <- "spatiotemporal"
save(all_rfs_sub, file = paste0(path, "fits/", "rf_data.RData"))

rfs_plot2D <- ggplot(all_rfs_sub, aes(term, estimate, fill = model)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width=.2, position=position_dodge(.9), 
                colour = "grey40", linewidth = 2) +
  scale_fill_manual(values = c("grey", "black")) +
  labs(x = "Random effect", y = "SD Estimate", fill = "Model") +
  facet_wrap(~ model_comp) +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = "top",
        #legend.position = c(.8,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "fits/","rf_2D.jpg"), plot = rfs_plot2D)
saveRDS(rfs_plot2D, file = paste0(path, "fits/", "hom_rf.rds"))

### MARGINAL model effects -----
### Base model
# year
y_base2a <- visreg_delta(mod_spattemp, xvar = "fyear", 
                        model = 1, scale = "response", gg = TRUE, 
                        partial = FALSE, rug = FALSE)
y_base2b <- visreg_delta(mod_spattemp, xvar = "fyear", 
                        model = 2, trans = exp, gg = TRUE, 
                        partial = FALSE, rug = FALSE)
y_base2adat <- visreg_delta(mod_spattemp, xvar = "fyear", 
                         model = 1, scale = "response", gg = TRUE, 
                         partial = FALSE, rug = FALSE, plot = FALSE)
y_base2bdat <- visreg_delta(mod_spattemp, xvar = "fyear", 
                         model = 2, trans = exp, gg = TRUE, 
                         partial = FALSE, rug = FALSE, plot = FALSE)
y_p_base <- format_mareff_plots(mareff1 = y_base2a, mareff2 = y_base2b)
ggsave(filename = paste0(path, "fits/","year_effect_base.jpg"), plot = y_p_base)

# gear
g_base2a <- visreg_delta(mod_spattemp, xvar = "fgear", 
                         model = 1, scale = "response", gg = TRUE, 
                         partial = FALSE, rug = FALSE)
g_base2b <- visreg_delta(mod_spattemp, xvar = "fgear", 
                         model = 2, trans = exp, gg = TRUE, 
                         partial = FALSE, rug = FALSE)
g_p_base <- format_mareff_plots(mareff1 = g_base2a, mareff2 = g_base2b)
ggsave(filename = paste0(path, "fits/","gear_effect_base.jpg"), plot = g_p_base)


# Covariate model
# year
y_cova <- visreg_delta(mod_mldxsst, xvar = "fyear", model = 1, 
                      gg = TRUE)
y_covb <- visreg_delta(mod_mldxsst, xvar = "fyear", model = 2, 
                       gg = TRUE)
y_cov2a <- visreg_delta(mod_mldxsst, xvar = "fyear", 
                        model = 1, scale = "response", gg = TRUE, 
                        partial = FALSE, rug = FALSE)
y_cov2b <- visreg_delta(mod_mldxsst, xvar = "fyear", 
                        model = 2, trans = exp, gg = TRUE, 
                        partial = FALSE, rug = FALSE)
y_cov2adat <- visreg_delta(mod_mldxsst, xvar = "fyear", 
                            model = 1, scale = "response", gg = TRUE, 
                            partial = FALSE, rug = FALSE, plot = FALSE)
y_cov2bdat <- visreg_delta(mod_mldxsst, xvar = "fyear", 
                            model = 2, trans = exp, gg = TRUE, 
                            partial = FALSE, rug = FALSE, plot = FALSE)
y_p <- format_mareff_plots(mareff1 = y_cov2a, mareff2 = y_cov2b)
ggsave(filename = paste0(path, "fits/","year_effect_cov.jpg"), plot = y_p)

# gear
g_cova <- visreg_delta(mod_mldxsst, xvar = "fgear", model = 1,
                      gg = TRUE)
g_covb <- visreg_delta(mod_mldxsst, xvar = "fgear", model = 2, 
                      gg = TRUE)
g_cov2a <- visreg_delta(mod_mldxsst, xvar = "fgear", model = 1,
                        scale = "response", gg = TRUE, partial = FALSE, rug = FALSE)
g_cov2b <- visreg_delta(mod_mldxsst, xvar = "fgear", model = 2,
                       trans = exp, gg = TRUE, partial = FALSE, rug = FALSE)
g_p <- format_mareff_plots(mareff1 = g_cov2a, mareff2 = g_cov2b)
ggsave(filename = paste0(path, "fits/","gear_effect.jpg"), plot = g_p)

# mld x sst
d_mld_cova <- visreg2d_delta(mod_mldxsst, model = 1,
                            xvar = "mld_scaled", yvar = "sst_scaled", 
                            plot.type="persp")
d_mld_covb <- visreg2d_delta(mod_mldxsst, model = 2,
                             xvar = "mld_scaled", yvar = "sst_scaled", 
                            plot.type="persp")
d_mld_cov2a <- visreg2d_delta(mod_mldxsst, model = 1,
                              xvar = "mld_scaled", yvar = "sst_scaled",
                              scale = "response", plot.type="persp")
d_mld_cov2b <- visreg2d_delta(mod_mldxsst, model = 2,
                              xvar = "mld_scaled", yvar = "sst_scaled", 
                             trans = exp, plot.type="persp")

png(paste0(path, "fits/", "mld_sst_encounter.png"))
plot(d_mld_cov2a, main = "Encounter probability", 
     cex.main = 1.5, cex.lab = 1.3)
dev.off()
png(paste0(path, "fits/", "mld_sst_biomass.png"))
plot(d_mld_cov2b, main = "Biomass", 
     cex.main = 1.5, cex.lab = 1.3)
dev.off()

# save all
mar_eff_list <- list("fyear_m1" = y_cova, "fyear_m2" = y_covb,
                     "fyear2_m1" = y_cov2a, "fyear2_m2" = y_cov2b,
                     "fgear_m1" = g_cova, "fgear_m2" = g_covb, 
                     "fgear2_m1" = g_cov2a, "fgear2_m2" = g_cov2b,
                     "mld_sst_m1" = d_mld_cova, "mld_sst_m2" = d_mld_covb,
                     "mld_sst2_m1" = d_mld_cov2a, "mld_sst2_m2" = d_mld_cov2b)
saveRDS(mar_eff_list, paste0(path, "fits/", "allmarginaleffects.rds"))

year_eff_list <- list("fyear2_base_m1_dat" = y_base2adat, 
                      "fyear2_base_m2_dat" = y_base2bdat,
                      "fyear2_cov_m1_dat" = y_cov2adat, 
                      "fyear2_cov_m2_dat" = y_cov2bdat)
saveRDS(year_eff_list, paste0(path, "fits/", "year_marginaleffects.rds"))


### Partial deviance explained (PDE) -----
mod_sst <- mods$m_dlogn_sst[[1]]
#mod_depthxmld <- mods$m_dlogn_dxmld[[1]]
#mod_depth <- mods$m_dlogn_d[[1]]
mod_mld <- mods$m_dlogn_mld[[1]]


null_model <- sdmTMB(
  data = mod_spattemp$data,
  list(numkm ~ 1,
       numkm ~ 1),
  mesh = mod_spattemp$spde,
  family = delta_lognormal(),
  offset = log(mod_spattemp$data$sweptareakm2)
)

sink(file = paste0(path, "fits/", "pde_report.txt"))
deviance_report(full_cov_model = mod_mldxsst, null_model = null_model, 
                base_model = mod_spattemp, 
                reduced_models = list(mod_sst, mod_mld))
sink(file = NULL)

# Check significance of fixed covariates -----
# https://github.com/pbs-assess/sdmTMB/discussions/360
# https://github.com/pbs-assess/sdmTMB/discussions/228
# base - year effects
SDrepFixed_base <- as.data.frame(summary(mod_spattemp$sd_report, 
                                         p.value=TRUE, select = "fixed"))
SDrepFixed_base$Variable <- row.names(SDrepFixed_base)
SDrepFixed_base$Variable <- gsub("\\.[0-9]+", "", SDrepFixed_base$Variable)
row.names(SDrepFixed_base) <- NULL
table(SDrepFixed_base$Variable)

terms_base <- rbind(tidy(mod_spattemp, model = 1), 
                    tidy(mod_spattemp, model = 2))
print(terms_base, n = nrow(terms_base))
SDrepFixedFac_base <- SDrepFixed_base[SDrepFixed_base$Variable %in% c("b_j", "b_j2"),]
SDrepFixedFac_base$ModelComp <- ifelse(SDrepFixedFac_base$Variable == "b_j",
                                       1, 2)

nrow(SDrepFixedFac_base) == nrow(terms_base)
all_terms_base <- cbind(terms_base[, "term"], SDrepFixedFac_base[, c(1:4, 6)], 
                        terms_base[, c("conf.low", "conf.high")])
all_terms_base$sig001 <- ifelse(all_terms_base$`Pr(>|z^2|)` < 0.001,
                                "sig", "insig")
all_terms_base$sig05 <- ifelse(all_terms_base$`Pr(>|z^2|)` < 0.05,
                               "sig", "insig")
write.csv(all_terms_base, 
          file = paste0(path, "fits/", "base_factor_level_pvalues.csv"))

# cov - year + covariates
# factor effects
SDrepFixed_cov <- as.data.frame(summary(mod_mldxsst$sd_report, 
                                        p.value=TRUE, select = "fixed"))
SDrepFixed_cov$Variable <- row.names(SDrepFixed_cov)
SDrepFixed_cov$Variable <- gsub("\\.[0-9]+", "", SDrepFixed_cov$Variable)
row.names(SDrepFixed_cov) <- NULL
table(SDrepFixed_cov$Variable)
SDrepFixedFac_cov <- SDrepFixed_cov[SDrepFixed_cov$Variable %in% c("b_j", "b_j2"),]
SDrepFixedFac_cov$ModelComp <- ifelse(SDrepFixedFac_cov$Variable == "b_j", 
                                      1, 2)

terms_cov <- rbind(tidy(mod_mldxsst, model = 1), 
                   tidy(mod_mldxsst, model = 2))
terms_cov <- terms_cov[-c(25:26, 51:52),]
print(terms_cov, n = nrow(terms_cov))

nrow(SDrepFixedFac_cov) == nrow(terms_cov)
all_terms_cov <- cbind(terms_cov[, "term"], SDrepFixedFac_cov[, c(1:4, 6)], 
                       terms_cov[, c("conf.low", "conf.high")])
all_terms_cov$sig001 <- ifelse(all_terms_cov$`Pr(>|z^2|)` < 0.001,
                               "sig", "insig")
all_terms_cov$sig05 <- ifelse(all_terms_cov$`Pr(>|z^2|)` < 0.05,
                              "sig", "insig")
write.csv(all_terms_cov, 
          file = paste0(path, "fits/", "cov_factor_level_pvalues.csv"))

########## -------------------------
### TESTING MODEL DESIGN AND FIT
########## -------------------------

### CLEAR SESSION

path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/HOM/"
mods <- readRDS(paste0(path, "fits/final_model_object.rds"))
mod_spattemp <- mods$m_dlogn_spattemp[[1]]
mod_mldxsst <- mods$m_dlogn_mldxsst[[1]]
rm(mods)
gc()


### Crossvalidation ------

# separate data in test and training folds
set.seed(5)
seeds <- sample(1:1500, 20, replace = TRUE)
partitions <- holdout(model = mod_spattemp, 
                      training_cutoff = 0.7, seeds = seeds)

# run CVs and calculate evaluation scores
baseCVs <- map_dfr(partitions, 
                   rep_cv, formula = "base", 
                   .id = "randomPartition_nr")
head(baseCVs)
colnames(baseCVs)[2:3] <- c("RMSE_base", "MAE_base")

covCVs <- map_dfr(partitions, 
                  rep_cv, formula = "cov", 
                  .id = "randomPartition_nr")
head(covCVs)
colnames(covCVs)[2:3] <- c("RMSE_cov", "MAE_cov")

CV_results <- merge(baseCVs, covCVs)
CV_results$RMSE_dif <- CV_results$RMSE_cov - CV_results$RMSE_base
CV_results$MAE_dif <- CV_results$MAE_cov - CV_results$MAE_base
CV_results$RMSE_dif_mean <- mean(CV_results$RMSE_dif)
CV_results$MAE_dif_mean <- mean(CV_results$MAE_dif)

write.table(CV_results, 
            file = paste0(path, "fits/", "crossvalidation_results.txt"), 
            append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)

### Check residuals ------
#Looking at tails
# randomised quantile residuals - normal
mod_base_res1 <- residuals(mod_spattemp, model = 1, type = "mle-mvn")
mod_base_res1 <- mod_base_res1[!is.na(mod_base_res1) & !is.infinite(mod_base_res1)]
mod_base_res2 <- residuals(mod_spattemp, model = 2, type = "mle-mvn")
mod_base_res2 <- mod_base_res2[!is.na(mod_base_res2) & !is.infinite(mod_base_res2)]

mod_cov_res1 <- residuals(mod_mldxsst, model = 1, type = "mle-mvn")
mod_cov_res1 <- mod_cov_res1[!is.na(mod_cov_res1) & !is.infinite(mod_cov_res1)]
mod_cov_res2 <- residuals(mod_mldxsst, model = 2, type = "mle-mvn")
mod_cov_res2 <- mod_cov_res2[!is.na(mod_cov_res2) & !is.infinite(mod_cov_res2)]


png(paste0(path, "fits/", "randomised_quantile_base_bestfit.png"))
par(mfrow = c(2,2))
qqnorm(mod_base_res1, main = "Base model - binomial: Normal Q-Q Plot");abline(0, 1)
qqnorm(mod_base_res2, main = "Base model - lognormal: Normal Q-Q Plot");abline(0, 1)
qqnorm(mod_cov_res1, main = "Covariate model - binomial: Normal Q-Q Plot");abline(0, 1)
qqnorm(mod_cov_res2, main = "Covariate model - lognormal: Normal Q-Q Plot");abline(0, 1)
dev.off()

# simulation-based - uniform
s_mod_base <- make_dharma(model = mod_spattemp, model_comp = c(1, 2))
s_mod_cov <- make_dharma(model = mod_mldxsst, model_comp = c(1, 2))

# Uniformity, dispersion, outliers
png(paste0(path, "fits/", "dharma_covbestfit_test_m1.png"))
sink(file = paste0(path, "fits/", "dharma_standardtests_covmod_m1.txt"))
DHARMa::testResiduals(s_mod_cov$m1)
sink(file = NULL)
dev.off()
png(paste0(path, "fits/", "dharma_covbestfit_test_m2.png"))
sink(file = paste0(path, "fits/", "dharma_standardtests_covmod_m2.txt"))
DHARMa::testResiduals(s_mod_cov$m2)
sink(file = NULL)
dev.off()
png(paste0(path, "fits/", "dharma_base_test_m1.png"))
sink(file = paste0(path, "fits/", "dharma_standardtests_basemod_m1.txt"))
DHARMa::testResiduals(s_mod_base$m1)
sink(file = NULL)
dev.off()
png(paste0(path, "fits/", "dharma_base_test_m2.png"))
sink(file = paste0(path, "fits/", "dharma_standardtests_basemod_m2.txt"))
DHARMa::testResiduals(s_mod_base$m2)
sink(file = NULL)
dev.off()

# Res vs preds
png(paste0(path, "fits/", "dharma_mods_vs_preds.png"))
par(mfrow = c(2,2))
DHARMa::plotResiduals(s_mod_base$m1)
DHARMa::plotResiduals(s_mod_base$m2)
DHARMa::plotResiduals(s_mod_cov$m1)
DHARMa::plotResiduals(s_mod_cov$m2)
title(main = list("Base model (upper), Covariate model (lower)", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()

# Residuals vs covariates
mod_plot_data <- subset(mod_mldxsst$data, numkm > 0)
png(paste0(path, "fits/", "dharma_bestfit_res_vs_covs.png"))
par(mfrow = c(3,2))
DHARMa::plotResiduals(s_mod_cov$m1, mod_mldxsst$data$mld_scaled)
DHARMa::plotResiduals(s_mod_cov$m2, mod_plot_data$mld_scaled)
DHARMa::plotResiduals(s_mod_cov$m1, mod_mldxsst$data$sst_scaled)
DHARMa::plotResiduals(s_mod_cov$m2, mod_plot_data$sst_scaled)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()

# Zeroinflation
png(paste0(path, "fits/", "zeroinflation_mods.png"))
par(mfrow = c(2,2))
DHARMa::testZeroInflation(s_mod_base$m1)
DHARMa::testZeroInflation(s_mod_cov$m1)
DHARMa::testZeroInflation(s_mod_base$m2)
DHARMa::testZeroInflation(s_mod_cov$m2)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/", "dharma_zeroinf_basemod_m1.txt"))
DHARMa::testZeroInflation(s_mod_base$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_zeroinf_basemod_m2.txt"))
DHARMa::testZeroInflation(s_mod_base$m2)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_zeroinf_covmod_m1.txt"))
DHARMa::testZeroInflation(s_mod_cov$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_zeroinf_covmod_m2.txt"))
DHARMa::testZeroInflation(s_mod_cov$m2)
sink(file = NULL)

# Uniformity
png(paste0(path, "fits/", "uniformity_mods.png"))
par(mfrow = c(2,2))
DHARMa::testUniformity(s_mod_base$m1)
DHARMa::testUniformity(s_mod_cov$m1)
DHARMa::testUniformity(s_mod_base$m2)
DHARMa::testUniformity(s_mod_cov$m2)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/", "dharma_uniformity_basemod_m1.txt"))
DHARMa::testUniformity(s_mod_base$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_uniformity_basemod_m2.txt"))
DHARMa::testUniformity(s_mod_base$m2)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_uniformity_covmod.txt"))
DHARMa::testUniformity(s_mod_cov$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_uniformity_covmod.txt"))
DHARMa::testUniformity(s_mod_cov$m2)
sink(file = NULL)


# Dispersion
png(paste0(path, "fits/", "dispersion_mods.png"))
par(mfrow = c(2,2))
DHARMa::testDispersion(s_mod_base$m1)
DHARMa::testDispersion(s_mod_cov$m1)
DHARMa::testDispersion(s_mod_base$m2)
DHARMa::testDispersion(s_mod_cov$m2)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/", "dharma_dispersion_basemod_m1.txt"))
DHARMa::testDispersion(s_mod_base$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_dispersion_basemod_m2.txt"))
DHARMa::testDispersion(s_mod_base$m2)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_dispersion_covmod.txt"))
DHARMa::testDispersion(s_mod_cov$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_dispersion_covmod.txt"))
DHARMa::testDispersion(s_mod_cov$m2)
sink(file = NULL)

# Temporal autocorrelation
grouped_basemod1 <- recalculateResiduals(s_mod_base$m1, group = unique(mod_mldxsst$data$fyear))
grouped_basemod2 <- recalculateResiduals(s_mod_base$m2, group = unique(mod_plot_data$fyear))

grouped_covmod1<- recalculateResiduals(s_mod_cov$m1, group = unique(mod_mldxsst$data$fyear))
grouped_covmod2 <- recalculateResiduals(s_mod_cov$m2, group = unique(mod_plot_data$fyear))

png(paste0(path, "fits/", "tempautocorr_basemod_m1.png"))
DHARMa::testTemporalAutocorrelation(grouped_basemod1, time = unique(mod_mldxsst$data$fyear))
title(main = list("Base model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
png(paste0(path, "fits/", "tempautocorr_basemod_m2.png"))
DHARMa::testTemporalAutocorrelation(grouped_basemod2, time = unique(mod_plot_data$fyear))
title(main = list("Base model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
png(paste0(path, "fits/", "tempautocorr_covmod_m1.png"))
DHARMa::testTemporalAutocorrelation(grouped_covmod1, time = unique(mod_mldxsst$data$fyear))
title(main = list("Covariate model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
png(paste0(path, "fits/", "tempautocorr_covmod_m2.png"))
DHARMa::testTemporalAutocorrelation(grouped_covmod2, time = unique(mod_plot_data$fyear))
title(main = list("Covariate model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/", "dharma_dispersion_basemod_m1.txt"))
testTemporalAutocorrelation(grouped_basemod1, time = unique(mod_mldxsst$data$fyear))
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_dispersion_basemod_m2.txt"))
testTemporalAutocorrelation(grouped_basemod2, time = unique(mod_plot_data$fyear))
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_dispersion_covmod_m1.txt"))
testTemporalAutocorrelation(grouped_covmod1, time = unique(mod_mldxsst$data$fyear))
sink(file = NULL)
sink(file = paste0(path, "fits/", "dharma_dispersion_covmod_m2.txt"))
testTemporalAutocorrelation(grouped_covmod2, time = unique(mod_plot_data$fyear))
sink(file = NULL)

## Spatial autocorrelation
mod_spattemp$data$resids <- residuals(mod_spattemp, model = 1, type = "mle-mvn")
mod_spattemp$data$resids2 <- residuals(mod_spattemp, model = 2, type = "mle-mvn")
mod_mldxsst$data$resids <- residuals(mod_mldxsst, model = 1, type = "mle-mvn")
mod_mldxsst$data$resids2 <- residuals(mod_mldxsst, model = 2, type = "mle-mvn")

spatial_basemod1 <- ggplot(mod_spattemp$data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/","spatial_res_basemod_m1.png"), plot = spatial_basemod1)
spatial_basemod2 <- ggplot(mod_spattemp$data, aes(X, Y, col = resids2)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/","spatial_res_basemod_m2.png"), plot = spatial_basemod2)
spatial_covmod1 <- ggplot(mod_mldxsst$data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/","spatial_res_covmod_m1.png"), plot = spatial_covmod1)
spatial_covmod2 <- ggplot(mod_mldxsst$data, aes(X, Y, col = resids2)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/","spatial_res_covmod_m2.png"), plot = spatial_covmod2)

#coords_mod_spattemp <- mod_spattemp$data |> 
#  select(year, fquarter, lat, lon, X, Y, resids) 
#head(coords_mod_spattemp)

#ggplot(coords_mod_spattemp, aes(lon, lat, shape = fquarter)) +
#  geom_point() +
#  facet_wrap(~ year)
#ggsave(filename = paste0(path, "fits/","survey_stations_map.png"), plot = last_plot())

#coords_mod_cov <- mod_mldxsst$data |> 
#  select(year, lat, lon, X, Y, resids)
#head(coords_mod_cov)

## Global Moran's I
# base
knw_spattemp_binom <- make_kn_dist_obj_binom(coords = mod_spattemp$data)
knw_spattemp_pos <- make_kn_dist_obj_pos(coords = mod_spattemp$data)

morans_spattemp_m1 <- purrr::map(knw_spattemp_binom, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, nsim = 999)
  #moran.test(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw)
))
morans_spattemp_m2 <- purrr::map(knw_spattemp_pos, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, nsim = 999)
  #moran.test(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw)
))

# we also need to record the sample size  
n_spattemp_m1 <- imap_dfr(knw_spattemp_binom, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})
n_spattemp_m2 <- imap_dfr(knw_spattemp_pos, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})

sink(file = paste0(path, "fits/","global_morans_results_base_binom.txt"))
cat("Sample size in each year \n")
n_spattemp_m1
cat("\n")
cat("Moran's I in each year \n")
morans_spattemp_m1
sink(file = NULL)
sink(file = paste0(path, "fits/","global_morans_results_base_logn.txt"))
cat("Sample size in each year \n")
n_spattemp_m2
cat("\n")
cat("Moran's I in each year \n")
morans_spattemp_m2
sink(file = NULL)

# finally we extract plots
morans_spattemp_m1_plots <- purrr::map(knw_spattemp_binom, ~ 
                                         moran.plot(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, return_df = TRUE)
)
morans_spattemp_m1_gg <- purrr::map(morans_spattemp_m1_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 1
))
globmorans_spattemp_m1_grid <- cowplot::plot_grid(plotlist = morans_spattemp_m1_gg,
                                                  align = "hv", axis = "tb")
buffered_spattemp_m1_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_spattemp_m1_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/", "globm_plots_base_binom.jpg"), 
       buffered_spattemp_m1_grid)

morans_spattemp_m2_plots <- purrr::map(knw_spattemp_pos, ~ 
                                         moran.plot(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, return_df = TRUE)
)
morans_spattemp_m2_gg <- purrr::map(morans_spattemp_m2_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 2
))
globmorans_spattemp_m2_grid <- cowplot::plot_grid(plotlist = morans_spattemp_m2_gg,
                                                  align = "hv", axis = "tb")
buffered_spattemp_m2_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_spattemp_m2_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/", "globm_plots_base_logn.jpg"), 
       buffered_spattemp_m2_grid)


# covariate
knw_cov_binom <- make_kn_dist_obj_binom(coords = mod_mldxsst$data)
knw_cov_pos <- make_kn_dist_obj_pos(coords = mod_mldxsst$data)

morans_cov_m1 <- purrr::map(knw_cov_binom, ~ list(
  #moran.test(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw)
  moran.mc(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, nsim = 999)
))
morans_cov_m2 <- purrr::map(knw_cov_pos, ~ list(
  #moran.test(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw)
  moran.mc(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, nsim = 999)
))

# we also need to record the sample size  
n_cov_m1 <- imap_dfr(knw_cov_binom, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})
n_cov_m2 <- imap_dfr(knw_cov_pos, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})

sink(file = paste0(path, "fits/","global_morans_results_cov_binom.txt"))
cat("Sample size in each year \n")
n_cov_m1
cat("\n")
cat("Moran's I in each year \n")
morans_cov_m1
sink(file = NULL)
sink(file = paste0(path, "fits/","global_morans_results_cov_logn.txt"))
cat("Sample size in each year \n")
n_cov_m2
cat("\n")
cat("Moran's I in each year \n")
morans_cov_m2
sink(file = NULL)

# finally we extract plots
morans_cov_m1_plots <- purrr::map(knw_cov_binom, ~ 
                                    moran.plot(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, return_df = TRUE)
)
morans_cov_m1_gg <- purrr::map(morans_cov_m1_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 1
))
globmorans_cov_m1_grid <- cowplot::plot_grid(plotlist = morans_cov_m1_gg,
                                             align = "hv", axis = "tb")
buffered_cov_m1_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_cov_m1_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/", "globm_plots_cov_binom.jpg"), 
       buffered_cov_m1_grid)

morans_cov_m2_plots <- purrr::map(knw_cov_pos, ~ 
                                    moran.plot(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, return_df = TRUE)
)
morans_cov_m2_gg <- purrr::map(morans_cov_m2_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 2
))
globmorans_cov_m2_grid <- cowplot::plot_grid(plotlist = morans_cov_m2_gg,
                                             align = "hv", axis = "tb")
buffered_cov_m2_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_cov_m2_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/", "globm_plots_cov_logn.jpg"), 
       buffered_cov_m2_grid)


## Local Moran's I
# apply functions to base and covariate models
library(parallel)
invisible(spdep::set.coresOption(max(detectCores()-1L, 1L)))

#base
locm_spattemp_binom <- purrr::map(knw_spattemp_binom, make_local_morans_obj_binom)
locm_spattemp_binom_df <- do.call(rbind, locm_spattemp_binom)
save(locm_spattemp_binom_df, 
     file = paste0(path, "fits/","local_morans_psignificance_base_binom.RData"))
locm_spattemp_binom_sig <- purrr::map(locm_spattemp_binom, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/","local_morans_psignificance_base_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_spattemp_binom_sig)
sink(file = NULL)

locm_spattemp_pos <- purrr::map(knw_spattemp_pos, make_local_morans_obj_pos)
locm_spattemp_pos_df <- do.call(rbind, locm_spattemp_pos)
save(locm_spattemp_pos_df, 
     file = paste0(path, "fits/","local_morans_psignificance_base_pos.RData"))
locm_spattemp_pos_sig <- purrr::map(locm_spattemp_pos, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/","local_morans_psignificance_base_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_spattemp_pos_sig)
sink(file = NULL)

# covariate
locm_cov_binom <- purrr::map(knw_cov_binom, make_local_morans_obj_binom)
locm_cov_binom_df <- do.call(rbind, locm_cov_binom)
save(locm_cov_binom_df, 
     file = paste0(path, "fits/","local_morans_psignificance_cov_binom.RData"))
locm_cov_binom_sig <- purrr::map(locm_cov_binom, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/","local_morans_psignificance_cov_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_cov_binom_sig)
sink(file = NULL)

locm_cov_pos <- purrr::map(knw_cov_pos, make_local_morans_obj_pos)
locm_cov_pos_df <- do.call(rbind, locm_cov_pos)
save(locm_cov_pos_df, 
     file = paste0(path, "fits/","local_morans_psignificance_cov_pos.RData"))
locm_cov_pos_sig <- purrr::map(locm_cov_pos, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/","local_morans_psignificance_cov_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_cov_pos_sig)
sink(file = NULL)

## Getis-Ord G
# base
getisord_spattemp_binom <- purrr::map(knw_spattemp_binom, make_getis_ord_obj_binom)
getisord_spattemp_binom_df <- do.call(rbind, getisord_spattemp_binom)
save(getisord_spattemp_binom_df, 
     file = paste0(path, "fits/","getis_ord_psignificance_spattemp_binom.RData"))
getisord_spattemp_binom_sig <- purrr::map(getisord_spattemp_binom, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/","getis_ord_psignificance_base_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_spattemp_binom_sig)
sink(file = NULL)

getisord_spattemp_pos <- purrr::map(knw_spattemp_pos, make_getis_ord_obj_pos)
getisord_spattemp_pos_df <- do.call(rbind, getisord_spattemp_pos)
save(getisord_spattemp_pos_df, 
     file = paste0(path, "fits/","getis_ord_psignificance_spattemp_pos.RData"))
getisord_spattemp_pos_sig <- purrr::map(getisord_spattemp_pos, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/","getis_ord_psignificance_base_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_spattemp_pos_sig)
sink(file = NULL)

# covariate
getisord_cov_binom <- purrr::map(knw_cov_binom, make_getis_ord_obj_binom)
getisord_cov_binom_df <- do.call(rbind, getisord_cov_binom)
save(getisord_cov_binom_df, 
     file = paste0(path, "fits/","getis_ord_psignificance_cov_binom.RData"))
getisord_cov_binom_sig <- purrr::map(getisord_cov_binom, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/","getis_ord_psignificance_cov_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_cov_binom_sig)
sink(file = NULL)

getisord_cov_pos <- purrr::map(knw_cov_pos, make_getis_ord_obj_pos)
getisord_cov_pos_df <- do.call(rbind, getisord_cov_pos)
save(getisord_cov_pos_df, 
     file = paste0(path, "fits/","getis_ord_psignificance_cov_pos.RData"))
getisord_cov_pos_sig <- purrr::map(getisord_cov_pos, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/","getis_ord_psignificance_cov_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_cov_pos_sig)
sink(file = NULL)

# Make a final dataframe encompassing all local autocorrelation measures
# Base model
# binomial
locm_spattemp_df_binom <- locm_spattemp_binom_df
rownames(locm_spattemp_df_binom) <- NULL
locg_spattemp_df_binom <- getisord_spattemp_binom_df
rownames(locg_spattemp_df_binom) <- NULL
nrow(locm_spattemp_df_binom) == nrow(locg_spattemp_df_binom)
local_autocor_spattemp_binom <- cbind(locm_spattemp_df_binom, 
                                      locg_spattemp_df_binom[, c("Gi", "E.Gi",
                                                                 "Var.Gi", "StdDev.Gi",
                                                                 "Pr(z != E(Gi))",
                                                                 "Pr(z != E(Gi)) Sim",
                                                                 "locg_p_pv", 
                                                                 "locg_significance",
                                                                 "locg_p_clust")])
# lognormal
locm_spattemp_df_logn <- locm_spattemp_pos_df
rownames(locm_spattemp_df_logn) <- NULL
locg_spattemp_df_logn <- getisord_spattemp_pos_df
rownames(locg_spattemp_df_logn) <- NULL
nrow(locm_spattemp_df_logn) == nrow(locg_spattemp_df_logn)
local_autocor_spattemp_logn <- cbind(locm_spattemp_df_logn, 
                                     locg_spattemp_df_logn[, c("Gi", "E.Gi",
                                                               "Var.Gi", "StdDev.Gi",
                                                               "Pr(z != E(Gi))",
                                                               "Pr(z != E(Gi)) Sim",
                                                               "locg_p_pv", 
                                                               "locg_significance",
                                                               "locg_p_clust")])

# Covariate model
# binomial
locm_cov_df_binom <- locm_cov_binom_df
rownames(locm_cov_df_binom) <- NULL
locg_cov_df_binom <- getisord_cov_binom_df
rownames(locg_cov_df_binom) <- NULL
nrow(locm_cov_df_binom) == nrow(locg_cov_df_binom)
local_autocor_cov_binom <- cbind(locm_cov_df_binom, 
                                 locg_cov_df_binom[, c("Gi", "E.Gi", "Var.Gi", 
                                                       "StdDev.Gi", "Pr(z != E(Gi))",
                                                       "Pr(z != E(Gi)) Sim",
                                                       "locg_p_pv", 
                                                       "locg_significance",
                                                       "locg_p_clust")])
# lognormal
locm_cov_df_logn <- locm_cov_pos_df
rownames(locm_cov_df_logn) <- NULL
locg_cov_df_logn <- getisord_cov_pos_df
rownames(locg_cov_df_logn) <- NULL
nrow(locm_cov_df_logn) == nrow(locg_cov_df_logn)
local_autocor_cov_logn <- cbind(locm_cov_df_logn, 
                                locg_cov_df_logn[, c("Gi", "E.Gi", "Var.Gi", 
                                                     "StdDev.Gi", "Pr(z != E(Gi))",
                                                     "Pr(z != E(Gi)) Sim",
                                                     "locg_p_pv", 
                                                     "locg_significance",
                                                     "locg_p_clust")])

# Save all outputs 
spat_autocor_res <- list("base-binom" = local_autocor_spattemp_binom,
                         "base-lognormal" = local_autocor_spattemp_logn,
                         "cov-binom" = local_autocor_cov_binom,
                         "cov-lognormal" = local_autocor_cov_logn)
saveRDS(spat_autocor_res, file = paste0(path, "fits/", "spatialAutocorResults.rds"))

# Plot out
# Make the Voronoi tesselation and hull
neatl <- st_read("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/HOM/plotting aids/neatl_canvas.shp")

h_v_spattemp_binom <- make_hull_and_vortess(local_autocor_spattemp_binom, neatl)
h_v_spattemp_logn <- make_hull_and_vortess(local_autocor_spattemp_logn, neatl)
h_v_cov_binom <- make_hull_and_vortess(local_autocor_cov_binom, neatl)
h_v_cov_logn <- make_hull_and_vortess(local_autocor_cov_logn, neatl)

# plotting local Moran's
# Base
# binomial
all_locm_plots_base_binom <- purrr::map(h_v_spattemp_binom, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_base_binom <- cowplot::plot_grid(plotlist = all_locm_plots_base_binom)
cowplot::save_plot(paste0(path, "fits/", "locm_plots_base_binom.jpg"), 
                   locm_plots_base_binom, base_asp = 1.6)

# lognormal
all_locm_plots_base_logn <- purrr::map(h_v_spattemp_logn, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_base_logn <- cowplot::plot_grid(plotlist = all_locm_plots_base_logn)
cowplot::save_plot(paste0(path, "fits/", "locm_plots_base_logn.jpg"), 
                   locm_plots_base_logn, base_asp = 1.6)


# Covariate
# binomial
all_locm_plots_cov_binom <- purrr::map(h_v_cov_binom, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_cov_binom <- cowplot::plot_grid(plotlist = all_locm_plots_cov_binom)
cowplot::save_plot(paste0(path, "fits/", "locm_plots_cov_binom.jpg"), 
                   locm_plots_cov_binom, base_asp = 1.6)

# lognormal
all_locm_plots_cov_logn <- purrr::map(h_v_cov_logn, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_cov_logn <- cowplot::plot_grid(plotlist = all_locm_plots_cov_logn)
cowplot::save_plot(paste0(path, "fits/", "locm_plots_cov_logn.jpg"), 
                   locm_plots_cov_logn, base_asp = 1.6)


# plotting Getis Ord
# Base
# binomial
all_locg_plots_base_binom <- purrr::map(h_v_spattemp_binom, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_base_binom <- cowplot::plot_grid(plotlist = all_locg_plots_base_binom)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_plots_base_binom.jpg"), 
                   gord_plots_base_binom, base_asp = 1.6)

# lognormal
all_locg_plots_base_logn <- purrr::map(h_v_spattemp_logn, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_base_logn <- cowplot::plot_grid(plotlist = all_locg_plots_base_logn)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_plots_base_logn.jpg"), 
                   gord_plots_base_logn, base_asp = 1.6)

# Covariate
# binomial
all_locg_plots_cov_binom <- purrr::map(h_v_cov_binom, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_cov_binom <- cowplot::plot_grid(plotlist = all_locg_plots_cov_binom)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_plots_cov_binom.jpg"), 
                   gord_plots_cov_binom, base_asp = 1.6)

# lognormal
all_locg_plots_cov_logn <- purrr::map(h_v_cov_logn, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_cov_logn <- cowplot::plot_grid(plotlist = all_locg_plots_cov_logn)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_plots_cov_logn.jpg"), 
                   gord_plots_cov_logn, base_asp = 1.6)

# plotting Getis Ord hotspots
# Base
# binomial
all_locg_hotspot_plots_base_binom <- purrr::map(h_v_spattemp_binom, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_base_binom <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_base_binom)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_hotspot_plots_base_binom.jpg"), 
                   hotspot_plots_base_binom, base_asp = 1.6)

# lognormal
all_locg_hotspot_plots_base_logn <- purrr::map(h_v_spattemp_logn, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_base_logn <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_base_logn)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_hotspot_plots_base_logn.jpg"), 
                   hotspot_plots_base_logn, base_asp = 1.6)

# Covariate
# binomial
all_locg_hotspot_plots_cov_binom <- purrr::map(h_v_cov_binom, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_cov_binom <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_cov_binom)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_hotspot_plots_cov_binom.jpg"), 
                   hotspot_plots_cov_binom, base_asp = 1.6)

# lognormal
all_locg_hotspot_plots_cov_logn <- purrr::map(h_v_cov_logn, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_cov_logn <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_cov_logn)
cowplot::save_plot(paste0(path, "fits/", "getis_ord_hotspot_plots_cov_logn.jpg"), 
                   hotspot_plots_cov_logn, base_asp = 1.6)


### Map spatial and spatiotemporal variability -----
### Make predictions ------
load("prediction_grid.RData")
head(pred_grid)
pred_grid <- add_utm_columns(pred_grid, c("lon", "lat"), units = "km")
grid_yrs <- data.frame(X = pred_grid$X,
                       Y = pred_grid$Y,
                       sst_scaled = pred_grid$sst_scaled, 
                       mld_scaled = pred_grid$mld_scaled,
                       fyear = as.factor(pred_grid$fyear),
                       fgear = as.factor(pred_grid$gear),
                       year = as.numeric(pred_grid$year))
grid_yrs <- distinct(grid_yrs)

rm("pred_grid")

path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/HOM/"
mods <- readRDS(paste0(path, "fits/final_model_object.rds"))
mod_spattemp <- mods$m_dlogn_spattemp[[1]]
mod_mldxsst <- mods$m_dlogn_mldxsst[[1]]
rm(mods)
gc()

neatl <- st_read("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/HOM/plotting aids/neatl_canvas.shp")
neatl_utm <- st_transform(neatl, st_crs(32629)) # get shapefile for plotting

# load functions
source("funs_hom.R")

# With spatial random effects and spatiotemporal rf in both predictors
predictions_spattemp <- predict(mod_spattemp, newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                                offset = rep(0, nrow(grid_yrs)), 
                                return_tmb_object = TRUE)
# Including all covariates
predictions_mldsst <- predict(mod_mldxsst, newdata = grid_yrs, 
                             offset = rep(0, nrow(grid_yrs)),
                             return_tmb_object = TRUE)

# Save predictions
density_base <- predict(mod_spattemp, newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                     offset = rep(0, nrow(grid_yrs)), type = "response")
predictions_spattemp$data$density <- density_base$est
saveRDS(predictions_spattemp, file = paste0(path, "fits/","fullobj_base_model_predictions.rds"))
write.csv(predictions_spattemp$data, file = paste0(path, "fits/","base_model_predictions.csv"), row.names = FALSE)

density_cov <- predict(mod_mldxsst, newdata = grid_yrs, 
                    offset = rep(0, nrow(grid_yrs)), type = "response")
predictions_mldsst$data$density <- density_cov$est
saveRDS(predictions_mldsst, file = paste0(path, "fits/","fullobj_cov_model_predictions.rds"))
write.csv(predictions_mldsst$data, file = paste0(path, "fits/","cov_model_predictions.csv"), row.names = FALSE)

### Map outputs ------
base_preds <- readRDS(paste0(path, "fits/","fullobj_base_model_predictions.rds"))
cov_preds <- readRDS(paste0(path, "fits/","fullobj_cov_model_predictions.rds"))

neatl_pol <- map_data("world", region = c("UK", "Ireland", 
                                                   "France", "Spain", 
                                                   "Portugal", "Norway"))
neatl_pol <- add_utm_columns(neatl_pol, c("long", "lat"), units = "km")

# DENSITY
base_density <- plot_map(base_preds$data, density) + 
  ggtitle("Base model: Density") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  #geom_polygon(data = neatl_pol, aes(long, lat, group = group), fill = "grey67") +
  theme_map()
ggsave(filename = paste0(path, "fits/","density_base.jpg"), plot = base_density, dpi = 300)

cov_density <- plot_map(predictions_mldsst$data, density) + 
  ggtitle("Covariate model: Density") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  theme_map()
ggsave(filename = paste0(path, "fits/","density_cov.jpg"), plot = cov_density, dpi = 300)

# Average density
ggplot(neatl_pol, aes(long, lat, group = group), fill = "grey67") +
  geom_polygon() +
  coord_quickmap(xlim = c(-10, 10), ylim = c(40, 65))
ggsave(filename = paste0(path, "fits/","neatl_check.jpg"), 
       plot = last_plot(), dpi = 300)

mean_base_density <- base_preds$data |> 
  group_by(X, Y) |>
  summarise(mean_density = mean(density))

base_ave_density <- plot_map2(mean_base_density, mean_density) + 
  ggtitle("Base model: Density") +
  scale_fill_viridis_c(option = "G", direction = -1) #+
  #geom_polygon(data = neatl_pol, aes(X, Y, group = group), fill = "grey67") +
  #coord_quickmap(xlim = c(-15, 10), ylim = c(40, 65))
ggsave(filename = paste0(path, "fits/","average_density_base.jpg"), 
       plot = base_ave_density, dpi = 300)




# Uncertainty around model estimates
predictions_base_est1 <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")],
                                model = 1, nsim = 500, 
                                sims_var = "est")
grid_yrs$base_est1_se <- apply(predictions_base_est1, 1, sd) # uncertainty of spatial effect
rm(predictions_base_est1)

predictions_base_est2 <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")],
                                model = 2, nsim = 500, 
                                sims_var = "est")
grid_yrs$base_est2_se <- apply(predictions_base_est2, 1, sd)
rm(predictions_base_est2)

est_uncertainty_base1 <- plot_map(grid_yrs, base_est1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "fits/","est_se_base1.jpg"), plot = est_uncertainty_base1)

est_uncertainty_base2 <- plot_map(grid_yrs, base_est2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "fits/","est_se_base2.jpg"), plot = est_uncertainty_base2)


# covariate model
predictions_cov_est1 <- predict(mod_mldxsst, 
                               newdata = grid_yrs,
                               model = 1, nsim = 500, 
                               sims_var = "est")
grid_yrs$cov_est1_se <- apply(predictions_cov_est1, 1, sd) # uncertainty of spatial effect
rm(predictions_cov_est1)

predictions_cov_est2 <- predict(mod_mldxsst, 
                               newdata = grid_yrs, 
                               model = 2, nsim = 500, 
                               sims_var = "est")
grid_yrs$cov_est2_se <- apply(predictions_cov_est2, 1, sd)
rm(predictions_cov_est2)

est_uncertainty_cov1 <- plot_map(grid_yrs, cov_est1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "fits/","est_se_cov1.jpg"), plot = est_uncertainty_cov1)

est_uncertainty_cov2 <- plot_map(grid_yrs, cov_est2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "fits/","est_se_cov2.jpg"), plot = est_uncertainty_cov2)

# FIXED AN RANDOM EFFECTS
base_est1 <- plot_map(predictions_spattemp$data, plogis(est1)) + 
  ggtitle("Base model: Binomial model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Encounter\nprobability") +
  theme_map()
ggsave(filename = paste0(path, "fits/","enc_prob_base.jpg"), plot = base_est1, dpi = 300)

base_est2 <- plot_map(predictions_spattemp$data, exp(est2)) + 
  ggtitle("Base model: Lognormal model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Biomass\n(kg)") +
  theme_map()
ggsave(filename = paste0(path, "fits/","pos_catch_base.jpg"), plot = base_est2, dpi = 300)


cov_est1 <- plot_map(predictions_mldsst$data, plogis(est1)) + 
  ggtitle("Covariate model: Binomial model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Encounter\nprobability") +
  theme_map()
ggsave(filename = paste0(path, "fits/","enc_prob_cov.jpg"), plot = cov_est1)

cov_est2 <- plot_map(predictions_mldsst$data, exp(est2)) + 
  ggtitle("Covariate model: Lognormal model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Biomass\n(kg)") +
  theme_map()
ggsave(filename = paste0(path, "fits/","pos_catch_cov.jpg"), plot = cov_est2, dpi = 300)


# FIXED EFFECTS
base_ff1 <- plot_map(predictions_spattemp$data, est_non_rf1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component: Fixed effects only", 
       subtitle = "Fixed effects: fyear + fgear") +
  theme_map()
ggsave(filename = paste0(path, "fits/","ff_effects_base1.jpg"), plot = base_ff1)

base_ff2 <- plot_map(predictions_spattemp$data, est_non_rf2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component: Fixed effects only", 
       subtitle = "Fixed effects: fyear + fgear") +
  theme_map()
ggsave(filename = paste0(path, "fits/","ff_effects_base2.jpg"), plot = base_ff2)

cov_ff1 <- plot_map(predictions_mldsst$data, est_non_rf1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component: Fixed effects only", 
       subtitle = "Fixed effects: fyear + fgear +depth + MLD + SST") +
  theme_map()
ggsave(filename = paste0(path, "fits/","ff_effects_cov1.jpg"), plot = cov_ff1)

cov_ff2 <- plot_map(predictions_mldsst$data, est_non_rf2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component: Fixed effects only", 
       subtitle = "Fixed effects: fyear + fgear +depth + MLD + SST") +
  theme_map()
ggsave(filename = paste0(path, "fits/","ff_effects_cov2.jpg"), plot = cov_ff2)

# SPATIAL EFFECTS
base_omega1 <- plot_map2(base_preds$data, omega_s1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))

base_omega2 <- plot_map2(base_preds$data, omega_s2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))

cov_omega1 <- plot_map2(cov_preds$data, omega_s1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")) 

cov_omega2 <- plot_map2(cov_preds$data, omega_s2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))

omega <- (base_omega1 + base_omega2) / (cov_omega1 + cov_omega2) + 
  plot_annotation(title = "Spatial random effects only",
                  theme = theme(plot.title = element_text(hjust = 0.45, size = 16, face = "bold")))
ggsave(filename = paste0(path, "fits/","all_omega_base_cov.jpg"), plot = omega)

# Uncertainty around spatial random effects

### CLEAR SESSION
rm(list = c("predictions_spattemp", "predictions_d_s_m"))

# base model
predictions_base_om1 <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")],
                                model = 1, nsim = 500, 
                                sims_var = "omega_s")
grid_yrs$base_om1_se <- apply(predictions_base_om1, 1, sd) # uncertainty of spatial effect
rm(predictions_base_om1)

predictions_base_om2 <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")],
                                model = 2, nsim = 500, 
                                sims_var = "omega_s")
grid_yrs$base_om2_se <- apply(predictions_base_om2, 1, sd)
rm(predictions_base_om2)

omega_uncertainty_base1 <- plot_map2(grid_yrs, base_om1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "fits/","omega_sd_base1.png"), plot = omega_uncertainty_base1)

omega_uncertainty_base2 <- plot_map2(grid_yrs, base_om2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "fits/","omega_sd_base.png"), plot = omega_uncertainty_base)


# covariate model
predictions_cov_om1 <- predict(mod_mldxsst, 
                               newdata = grid_yrs,
                               model = 1, nsim = 500, 
                               sims_var = "omega_s")
grid_yrs$cov_om1_se <- apply(predictions_cov_om1, 1, sd) # uncertainty of spatial effect
rm(predictions_cov_om1)

predictions_cov_om2 <- predict(mod_mldxsst, 
                               newdata = grid_yrs, 
                               model = 2, nsim = 500, 
                               sims_var = "omega_s")
grid_yrs$cov_om2_se <- apply(predictions_cov_om2, 1, sd)
rm(predictions_cov_om2)

omega_uncertainty_cov1 <- plot_map2(grid_yrs, cov_om1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "fits/","omega_sd_base.png"), plot = omega_uncertainty_base)

omega_uncertainty_cov2 <- plot_map2(grid_yrs, cov_om2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "fits/","omega_sd_cov.png"), plot = omega_uncertainty_cov)

omega_se <- (omega_uncertainty_base1 + omega_uncertainty_base2) / 
  (omega_uncertainty_cov1 + omega_uncertainty_cov2) + 
  plot_annotation(title = "Point-wise uncertainty in spatial predictions",
                  theme = theme(plot.title = element_text(hjust = 0.45, size = 16, face = "bold")))
ggsave(filename = paste0(path, "fits/","all_omega_se_base_cov.jpg"), plot = omega_se)

# SPATIOTEMPORAL EFFECTS
base_preds <- read.csv(file = paste0(path, "fits/","base_model_predictions.csv"))
base_eps1 <- plot_map(base_preds, epsilon_st1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_effects_base1.jpg"), plot = base_eps1, dpi = 300)

base_eps2 <- plot_map(base_preds, epsilon_st2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_effects_base2.jpg"), plot = base_eps2, dpi = 300)

cov_preds <- read.csv(file = paste0(path, "fits/","cov_model_predictions.csv"))
cov_eps1 <- plot_map(cov_preds, epsilon_st1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_effects_cov1.jpg"), plot = cov_eps1, dpi = 300)

cov_eps2 <- plot_map(cov_preds, epsilon_st2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_effects_cov2.jpg"), plot = cov_eps2, dpi = 300)

# Uncertainty of spatiotemporal effects
# base model
predictions_base_eps1 <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                                model = 1, nsim = 500, 
                                sims_var = "epsilon_st")
grid_yrs$base_eps1_se <- apply(predictions_base_eps1, 1, sd) # uncertainty of spatial effect
rm(predictions_base_eps1)

predictions_base_eps2 <- predict(mod_spattemp, 
                                 newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                                 model = 1, nsim = 500, 
                                 sims_var = "epsilon_st")
grid_yrs$base_eps2_se <- apply(predictions_base_eps2, 1, sd)
rm(predictions_base_eps2)

eps_uncertainty_base1 <- plot_map(grid_yrs, base_eps1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binom comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_se_base1.jpg"), plot = eps_uncertainty_base1, dpi = 300)

eps_uncertainty_base2 <- plot_map(grid_yrs, base_eps2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, logn comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_se_base2.jpg"), plot = eps_uncertainty_base2, dpi = 300)


# covariate model
predictions_cov_eps1 <- predict(mod_mldxsst, 
                                newdata = grid_yrs,
                                model = 1, nsim = 500, 
                                sims_var = "epsilon_st")
grid_yrs$cov_eps1_se <- apply(predictions_cov_eps1, 1, sd) # uncertainty of spatial effect
rm(predictions_cov_eps1)

predictions_cov_eps2 <- predict(mod_mldxsst, 
                                newdata = grid_yrs,
                                model = 2, nsim = 500, 
                                sims_var = "epsilon_st")
grid_yrs$cov_eps2_se <- apply(predictions_cov_eps2, 1, sd)
rm(predictions_cov_eps2)

eps_uncertainty_cov1 <- plot_map(grid_yrs, cov_eps1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binom comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_se_cov1.jpg"), plot = eps_uncertainty_cov1)

eps_uncertainty_cov2 <- plot_map(grid_yrs, cov_eps2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, logn comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map()
ggsave(filename = paste0(path, "fits/","eps_se_cov2.jpg"), plot = eps_uncertainty_cov2)

# Save all predictions
save(grid_yrs, file = paste0(path, "fits/", "uncertainty_mapping_data.RData"))

